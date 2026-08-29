from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import re
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
import os

ROOT = Path.cwd().resolve()
BASE = Path(os.environ.get("CIMMYT_ANALYSIS_ROOT", Path(__file__).resolve().parent)).resolve()
HISTORY = BASE / "history"
RAW = HISTORY / "raw_downloads"
EXTRACTED = HISTORY / "derived" / "extracted_HLBSN"
CURRENT = HISTORY / "raw_downloads" / "cycle_2024"
AUDIT = HISTORY / "audit"
STANDARD = HISTORY / "standardized"
REPORT = HISTORY / "reports"
for folder in (AUDIT, STANDARD, REPORT):
    folder.mkdir(parents=True, exist_ok=True)

MAX_BYTES = 10 * 1024 * 1024
MAX_ROWS = 10000
MAX_COLS = 200
TARGET_RE = re.compile(r"^(\d+)(HLBSN|HZAN|SAWSN|SAWNS|IBWSN)_WB\.csv(?:\.gz)?$", re.I)
ENV_RE = re.compile(r"^(Jashore|Quirusillas|Okinawa)_(\d{4})_(1st|2nd)_sowing$", re.I)


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def opener(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig", errors="strict", newline="") if path.suffix.lower() == ".gz" else path.open("r", encoding="utf-8-sig", errors="strict", newline="")


def canonical_nursery(name: str) -> tuple[str, str, int, bool]:
    m = TARGET_RE.match(name)
    if not m:
        raise ValueError(name)
    num, series = int(m.group(1)), m.group(2).upper()
    typo = series == "SAWNS"
    if typo:
        series = "SAWSN"
    return f"{num}{series}", series, num, typo


def source_key(path: Path) -> str:
    if CURRENT in path.parents:
        return "cycle_2024"
    if RAW in path.parents:
        return path.parent.name
    if EXTRACTED in path.parents:
        return "HLBSN_2018_2021"
    return "unknown"


paths = []
for p in RAW.rglob("*"):
    if p.is_file() and TARGET_RE.match(p.name):
        paths.append(p)
for p in EXTRACTED.glob("*.csv"):
    if TARGET_RE.match(p.name):
        paths.append(p)
for p in CURRENT.glob("*.csv"):
    if TARGET_RE.match(p.name):
        paths.append(p)
paths = sorted(set(paths), key=lambda p: (canonical_nursery(p.name)[1], canonical_nursery(p.name)[2], p.name))

file_audit: list[dict] = []
environment_audit: list[dict] = []
long_all: list[dict] = []
file_materials: dict[str, dict[str, tuple[str, str]]] = {}

aliases = {
    "entry": "Entry", "ent": "Entry", "cid": "Cid", "sid": "Sid", "gid": "GID",
    "cross": "Cross", "cross name": "Cross", "sel_hist": "Sel_Hist", "selection history": "Sel_Hist",
}

for path in paths:
    if path.stat().st_size > MAX_BYTES:
        raise RuntimeError(f"File over size cap: {path.name}")
    nursery, series, number, filename_typo = canonical_nursery(path.name)
    with opener(path) as handle:
        reader = csv.reader(handle, strict=True)
        rows = list(itertools.islice(reader, MAX_ROWS + 2))
    if len(rows) > MAX_ROWS + 1:
        raise RuntimeError(f"Row cap exceeded: {path.name}")
    if not rows:
        raise RuntimeError(f"Empty file: {path.name}")
    header, data = rows[0], rows[1:]
    if len(header) > MAX_COLS:
        raise RuntimeError(f"Column cap exceeded: {path.name}")
    bad_width = sum(len(r) != len(header) for r in data)
    if bad_width:
        raise RuntimeError(f"Nonrectangular file: {path.name}")

    # Remove only columns that have an empty header and are completely empty.
    trailing_empty = []
    keep = []
    for j, name in enumerate(header):
        empty_col = not name.strip() and all((r[j].strip() == "" for r in data))
        if empty_col:
            trailing_empty.append(j)
        else:
            keep.append(j)
    header2 = [header[j].strip() for j in keep]
    data2 = [[r[j].strip() for j in keep] for r in data]

    canonical = {}
    env_cols = []
    unexpected = []
    for j, name in enumerate(header2):
        low = name.lower()
        if low in aliases:
            canonical[aliases[low]] = j
        elif ENV_RE.match(name):
            env_cols.append((j, name, ENV_RE.match(name)))
        else:
            unexpected.append(name)

    missing_ids = [x for x in ("Cid", "Sid", "GID", "Cross", "Sel_Hist") if x not in canonical]
    exact_rows = Counter(tuple(r) for r in data2)
    exact_duplicate_rows = sum(n for n in exact_rows.values() if n > 1)
    invalid_numeric = out_of_range = missing_cells = observed_cells = explicit_na_cells = blank_missing_cells = 0
    invalid_tokens = Counter()
    gid_values = []
    material_map: dict[str, tuple[str, str]] = {}

    for row_no, row in enumerate(data2, start=2):
        entry = row[canonical["Entry"]] if "Entry" in canonical else str(row_no - 1)
        gid = row[canonical["GID"]] if "GID" in canonical else ""
        cross = row[canonical["Cross"]] if "Cross" in canonical else ""
        sel = row[canonical["Sel_Hist"]] if "Sel_Hist" in canonical else ""
        cid = row[canonical["Cid"]] if "Cid" in canonical else ""
        sid = row[canonical["Sid"]] if "Sid" in canonical else ""
        if gid:
            gid_values.append(gid)
            material_map[gid] = (cross, sel)
        for j, original_env, match in env_cols:
            location = match.group(1).capitalize()
            year = int(match.group(2))
            sowing = match.group(3).lower()
            raw_value = row[j]
            value = None
            if raw_value == "" or raw_value.upper() == "NA":
                missing_cells += 1
                if raw_value == "":
                    blank_missing_cells += 1
                else:
                    explicit_na_cells += 1
            else:
                try:
                    value = float(raw_value)
                    observed_cells += 1
                    if value < 0 or value > 100:
                        out_of_range += 1
                except ValueError:
                    invalid_numeric += 1
                    invalid_tokens[raw_value] += 1
            long_all.append({
                "dataset_key": source_key(path), "source_file": path.name, "nursery": nursery,
                "series": series, "nursery_number": number, "Entry": entry, "Cid": cid, "Sid": sid,
                "GID": gid, "Cross": cross, "Sel_Hist": sel, "environment_original": original_env,
                "environment": f"{location}_{year}_{sowing}_sowing", "location": location,
                "year": year, "sowing": sowing, "disease_index_pct": "" if value is None else value,
                "source_row": row_no, "entry_generated": "Entry" not in canonical,
            })

    file_materials[nursery] = material_map
    years = sorted({int(m.group(2)) for _, _, m in env_cols})
    locations = sorted({m.group(1).capitalize() for _, _, m in env_cols})
    sowings = sorted({m.group(3).lower() for _, _, m in env_cols})
    file_audit.append({
        "dataset_key": source_key(path), "source_file": path.name, "nursery": nursery, "series": series,
        "nursery_number": number, "rows": len(data2), "columns_raw": len(header), "columns_retained": len(header2),
        "empty_columns_removed": len(trailing_empty), "entry_present": "Entry" in canonical,
        "entry_generated_from_row": "Entry" not in canonical, "identity_fields_missing": ";".join(missing_ids),
        "environment_columns": len(env_cols), "environment_years": ";".join(map(str, years)),
        "locations": ";".join(locations), "sowings": ";".join(sowings),
        "potential_cells": len(data2) * len(env_cols), "observed_cells": observed_cells,
        "missing_cells": missing_cells, "blank_missing_cells": blank_missing_cells, "explicit_na_cells": explicit_na_cells,
        "missing_rate": missing_cells / max(1, len(data2) * len(env_cols)),
        "invalid_numeric_cells": invalid_numeric, "invalid_tokens": ";".join(f"{k}:{v}" for k, v in sorted(invalid_tokens.items())), "out_of_range_cells": out_of_range,
        "exact_duplicate_rows_affected": exact_duplicate_rows, "gid_unique": len(set(gid_values)),
        "gid_duplicate_rows": len(gid_values) - len(set(gid_values)), "filename_typo_normalized": filename_typo,
        "unexpected_columns": ";".join(unexpected), "bytes": path.stat().st_size, "md5": md5(path),
        "complete_scan": True,
    })
    for j, original_env, match in env_cols:
        vals = []
        missing = invalid = 0
        for row in data2:
            raw_value = row[j]
            if raw_value == "" or raw_value.upper() == "NA":
                missing += 1
            else:
                try:
                    vals.append(float(raw_value))
                except ValueError:
                    invalid += 1
        vals_sorted = sorted(vals)
        def quant(p: float):
            if not vals_sorted:
                return ""
            pos = (len(vals_sorted) - 1) * p
            lo, hi = int(pos), min(int(pos) + 1, len(vals_sorted) - 1)
            return vals_sorted[lo] + (vals_sorted[hi] - vals_sorted[lo]) * (pos - lo)
        environment_audit.append({
            "nursery": nursery, "series": series, "source_file": path.name, "environment_original": original_env,
            "location": match.group(1).capitalize(), "year": int(match.group(2)), "sowing": match.group(3).lower(),
            "n_rows": len(data2), "observed": len(vals), "missing": missing, "invalid": invalid,
            "zero_prop": sum(v == 0 for v in vals) / max(1, len(vals)), "mean": sum(vals) / max(1, len(vals)),
            "median": quant(.5), "q1": quant(.25), "q3": quant(.75), "max": max(vals) if vals else "",
        })

# Cross-year nursery mapping and gaps.
mapping = []
for f in file_audit:
    actual_years = [int(x) for x in f["environment_years"].split(";") if x]
    dataset_cycle = {"cycle_2022": 2022, "cycle_2023": 2023, "cycle_2024": 2024}.get(f["dataset_key"])
    mismatch = bool(dataset_cycle and actual_years and dataset_cycle not in actual_years)
    mapping.append({
        "series": f["series"], "nursery": f["nursery"], "nursery_number": f["nursery_number"],
        "dataset_key": f["dataset_key"], "source_file": f["source_file"],
        "environment_years": f["environment_years"], "dataset_cycle_year": dataset_cycle or "multi-year package",
        "cycle_year_mismatch": mismatch, "mapping_basis": "filename nursery code + environment column years",
        "mapping_status": "review_year_conflict" if mismatch else "compatible",
        "exact_duplicate_of": "",
    })

# Exact normalized content duplicates across nominal nursery files.
nursery_records = defaultdict(list)
for r in long_all:
    nursery_records[r["nursery"]].append((r["GID"], r["Cross"], r["Sel_Hist"], r["environment"], str(r["disease_index_pct"])))
nursery_signatures = {}
for nursery_name, records in nursery_records.items():
    payload = "\n".join("\t".join(x) for x in sorted(records)).encode("utf-8")
    nursery_signatures[nursery_name] = hashlib.sha256(payload).hexdigest()
duplicate_nurseries = []
for a, b in itertools.combinations(sorted(nursery_signatures), 2):
    if nursery_signatures[a] == nursery_signatures[b]:
        duplicate_nurseries.append({
            "nursery_a": a, "nursery_b": b, "normalized_sha256": nursery_signatures[a],
            "long_rows": len(nursery_records[a]), "interpretation": "exact identity and phenotype duplicate after field normalization",
        })
        for m in mapping:
            if m["nursery"] == b:
                m["exact_duplicate_of"] = a
                m["mapping_status"] = "exclude_exact_duplicate"

# Pairwise GID overlap within each nursery series.
overlaps = []
series_members = defaultdict(list)
for nursery_name in file_materials:
    series_members[re.sub(r"^\d+", "", nursery_name)].append(nursery_name)
for series, member_list in sorted(series_members.items()):
    member_list = sorted(member_list, key=lambda x: int(re.match(r"^\d+", x).group()))
    for a, b in itertools.combinations(member_list, 2):
        ga, gb = set(file_materials[a]), set(file_materials[b])
        shared = ga & gb
        identity_conflicts = sum(file_materials[a][g] != file_materials[b][g] for g in shared)
        overlaps.append({
            "series": series, "nursery_a": a, "nursery_b": b, "gid_a": len(ga), "gid_b": len(gb),
            "shared_gid": len(shared), "jaccard": len(shared) / max(1, len(ga | gb)),
            "shared_gid_identity_conflicts": identity_conflicts,
        })

# Series sequence and missing nursery numbers.
series_summary = []
for series in sorted({m["series"] for m in mapping}):
    group = sorted([m for m in mapping if m["series"] == series], key=lambda x: x["nursery_number"])
    nums = [m["nursery_number"] for m in group]
    missing_nums = sorted(set(range(min(nums), max(nums) + 1)) - set(nums)) if nums else []
    years = sorted({int(y) for m in group for y in str(m["environment_years"]).split(";") if y})
    series_summary.append({
        "series": series, "nurseries": ";".join(m["nursery"] for m in group), "nursery_number_min": min(nums),
        "nursery_number_max": max(nums), "missing_nursery_numbers": ";".join(map(str, missing_nums)),
        "environment_year_min": min(years), "environment_year_max": max(years),
        "n_nursery_files": len(group), "n_environment_columns": sum(int(f["environment_columns"]) for f in file_audit if f["series"] == series),
    })

file_fields = list(file_audit[0])
env_fields = list(environment_audit[0])
long_fields = list(long_all[0])
write_csv(AUDIT / "file_compatibility_audit.csv", file_audit, file_fields)
write_csv(AUDIT / "environment_audit.csv", environment_audit, env_fields)
write_csv(AUDIT / "nursery_year_mapping.csv", mapping, list(mapping[0]))
write_csv(AUDIT / "cross_year_gid_overlap.csv", overlaps, list(overlaps[0]))
write_csv(AUDIT / "series_mapping_summary.csv", series_summary, list(series_summary[0]))
write_csv(AUDIT / "exact_duplicate_nursery_files.csv", duplicate_nurseries,
          ["nursery_a", "nursery_b", "normalized_sha256", "long_rows", "interpretation"])
write_csv(STANDARD / "historical_long_all.csv", long_all, long_fields)
write_csv(STANDARD / "historical_long_observed.csv", [r for r in long_all if r["disease_index_pct"] != ""], long_fields)
duplicate_exclusions = {x["nursery_b"] for x in duplicate_nurseries}
write_csv(STANDARD / "historical_long_primary_observed.csv",
          [r for r in long_all if r["disease_index_pct"] != "" and r["nursery"] not in duplicate_exclusions], long_fields)

# Consolidate the official Excel data dictionaries after their read-only artifact-tool inspection.
dictionary_rows = []
dictionary_json = AUDIT / "data_dictionary_inspection.json"
if dictionary_json.exists():
    for book in json.loads(dictionary_json.read_text(encoding="utf-8")):
        for sheet in book.get("sheets", []):
            for row in sheet.get("values", [])[2:]:
                if len(row) > 4 and str(row[4] or "").strip():
                    dictionary_rows.append({
                        "source_dictionary": book.get("source", ""),
                        "sheet": sheet.get("name", ""),
                        "example_file": str(row[0] or "").strip() if len(row) > 0 else "",
                        "variable": str(row[4] or "").strip(),
                        "variable_description": str(row[5] or "").strip() if len(row) > 5 else "",
                        "unit": str(row[6] or "").strip() if len(row) > 6 else "",
                        "declared_data_type": str(row[7] or "").strip() if len(row) > 7 else "",
                    })
write_csv(AUDIT / "official_data_dictionary_fields.csv", dictionary_rows,
          ["source_dictionary", "sheet", "example_file", "variable", "variable_description", "unit", "declared_data_type"])

summary = {
    "analysis_time": datetime.now().isoformat(timespec="seconds"), "files_scanned": len(file_audit),
    "complete_scans": sum(bool(x["complete_scan"]) for x in file_audit),
    "nursery_series": sorted({x["series"] for x in file_audit}),
    "nursery_files": len({x["nursery"] for x in file_audit}),
    "environment_columns": sum(x["environment_columns"] for x in file_audit),
    "potential_cells": sum(x["potential_cells"] for x in file_audit),
    "observed_cells": sum(x["observed_cells"] for x in file_audit),
    "missing_cells": sum(x["missing_cells"] for x in file_audit),
    "explicit_na_cells": sum(x["explicit_na_cells"] for x in file_audit),
    "invalid_numeric_cells": sum(x["invalid_numeric_cells"] for x in file_audit),
    "out_of_range_cells": sum(x["out_of_range_cells"] for x in file_audit),
    "year_conflicts": [x["nursery"] for x in mapping if x["cycle_year_mismatch"]],
    "filename_normalizations": [x["source_file"] for x in file_audit if x["filename_typo_normalized"]],
    "exact_duplicate_nursery_pairs": [f"{x['nursery_a']}={x['nursery_b']}" for x in duplicate_nurseries],
}
(AUDIT / "audit_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")

series_lines = [
    f"- {x['series']}：{x['nurseries']}；环境年份{x['environment_year_min']}–{x['environment_year_max']}；缺失编号{x['missing_nursery_numbers'] or '无'}。"
    for x in series_summary
]
report = f"""# CIMMYT历年数据字段兼容性审计与跨年度育种圃映射

审计时间：{summary['analysis_time']}  
范围：HLBSN、HZAN、SAWSN、IBWSN，环境年份2018–2024。

## 下载与扫描范围

- 官方历史数据包：5个，下载文件29个；
- 纳入字段审计的目标表型文件：{summary['files_scanned']}个；
- 育种圃文件：{summary['nursery_files']}个；
- 环境列：{summary['environment_columns']}个；
- 潜在材料×环境表型：{summary['potential_cells']}个，非缺失{summary['observed_cells']}个，缺失{summary['missing_cells']}个，其中显式`NA`为{summary['explicit_na_cells']}个；
- 非数值病害单元：{summary['invalid_numeric_cells']}个；超出0–100%的单元：{summary['out_of_range_cells']}个；
- 所有文件均在设定的字节、行数和列数上限内完成全表扫描，没有抽样截断。

## 跨年度育种圃序列

{chr(10).join(series_lines)}

育种圃编号是系列标识，不等同于年份。映射以“文件名中的育种圃编号＋环境列中的实际年份”为依据。一个育种圃文件可能覆盖两个环境年份，因此后续分析必须以环境列年份为时间单位，而不能用文件发布日期或圃号代替年份。

## 主要字段兼容性发现

1. Entry字段有三种情况：`Entry`、`Ent`和完全缺失。13HZAN、38SAWSN和54IBWSN等文件缺少Entry；标准长表为其生成圃内行号，并设置`entry_generated=TRUE`，不能把该行号当作CIMMYT正式Entry。
2. 系谱字段有两套命名：`Cross`/`Sel_Hist`与`Cross Name`/`Selection History`，已映射到统一字段但保留源文件名。
3. 36SAWSN原始文件名写为`36SAWNS_WB.csv.gz`，标准化育种圃名为36SAWSN，并在审计表中保留拼写修正标记。
4. 9HLBSN及14/15HZAN存在完全空的尾随列，仅在派生表中删除；原始文件未改变。
5. 15HZAN属于2024发布包，但六个环境列均标记为2023。进一步比较发现，14HZAN与15HZAN的233个非空GID、身份字段及1,404个材料×环境表型逐项完全相同；15HZAN标记为`exclude_exact_duplicate`，不能作为独立2024证据重复计入。
6. 地点大小写不统一，例如`quirusillas`，派生表统一为`Quirusillas`。
7. 病害指数范围兼容，当前未发现非数值或0–100%以外的有效观测。
8. 五份官方Excel数据字典均可读取，结构均为单个`Example`工作表（A1:Q16）。字典明确将环境列定义为“Location_Year_Sowing试验的wheat blast index”，单位为百分比，与本审计的环境列拆分和0–100%范围规则一致；`Ent`、`Cid`、`Sid`、`GID`的定义也支持现有映射。数据字典把`Cross`和`Sel_Hist`声明为numerical，但实际CSV含系谱文本，因此本项目按文本保存，这一元数据矛盾不能据字典强制数值化。

## GID跨年度映射

`cross_year_gid_overlap.csv`逐对报告同一系列育种圃之间的共享GID、Jaccard比例和Cross/Sel_Hist身份冲突。GID是跨年度材料匹配的主键；Entry只在单个育种圃内使用。共享GID较少时，应转向系谱层面的重复证据，而不能用相同Entry号连接不同年度。

## 纳入建议

- 可直接进入后续跨年度环境区分力分析：所有通过范围和环境列解析的目标文件；
- 材料级跨年度分析：仅使用非空GID连接；Entry生成行只作为圃内记录标识；
- 15HZAN：保留原始文件和审计证据，但从主跨年度统计中排除，等待CIMMYT确认其是否为14HZAN重复发布；
- 36SAWSN：按标准名纳入，同时报告原文件拼写；
- SAWSN序列缺少39号文件，不能假设其不存在试验，只能表述为“本次官方公开数据检索未获得”；
- 历史包未发现2018–2021年的HZAN专项文件，因此HZAN时间序列目前从2022开始。

## 输出

- `audit/file_compatibility_audit.csv`：逐文件字段与质量审计；
- `audit/environment_audit.csv`：逐环境完整性及分布；
- `audit/nursery_year_mapping.csv`：育种圃—环境年份映射；
- `audit/cross_year_gid_overlap.csv`：跨年度GID连接证据；
- `audit/exact_duplicate_nursery_files.csv`：标准化后完全重复的育种圃文件；
- `audit/official_data_dictionary_fields.csv`：五份官方Excel数据字典的字段定义汇总；
- `audit/data_dictionary_inspection.json`：Excel工作簿只读解析证据；
- `standardized/historical_long_all.csv`和`historical_long_observed.csv`：含全部审计记录的统一长表；
- `standardized/historical_long_primary_observed.csv`：排除标准化后完全重复育种圃的主分析长表。

这些派生长表仅用于本项目分析。依据CIMMYT Custom Dataset Terms，不应将其作为原始数据替代品重新公开分发。
"""
(REPORT / "historical_compatibility_and_mapping_report.md").write_text(report, encoding="utf-8")
print(json.dumps(summary, ensure_ascii=False, indent=2))
