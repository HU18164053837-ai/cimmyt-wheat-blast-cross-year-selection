# End-to-end run order

The package is self-contained. Run from any location with PowerShell, Python, and R 4.4.1 available. The master runner resolves the project root from its own file location and sets `CIMMYT_ANALYSIS_ROOT` for every script.

```powershell
./RUN_ALL.ps1 -Python "path/to/python" -Rscript "path/to/Rscript"
```

The runner writes `运行记录/end_to_end_rerun.log`, `运行记录/run_manifest.csv`, and `运行记录/key_file_checksums.csv`, including script SHA256 values, exit status, elapsed time, locked key counts, and key-file hashes.

The manual order is:

1. `04_历年字段审计与育种圃映射.py`
2. `05_跨年度环境与候选稳健性分析.R`
3. `06_证据分层与系谱来源分析.R`
4. `07_候选复验组合与敏感性分析.R`
5. `08_共同对照锚定优势分析.R`
6. `09_时间外伪前瞻验证.R`
7. `10_环境校正材料效应与失效边界.R`
8. `11_候选阈值敏感性与基线比较.R`
9. `12_跨年度不确定性与负对照.R`
10. `13_研究流程图.R`
11. `14_编辑意见方法透明度与敏感性分析.R`

Expected locked counts after step 7 are 30,475 records, 4,088 linkable GIDs, 104 cross-year-eligible GIDs, 55 strong or moderate GIDs, 15 priority GIDs, and 14 nonredundant resistance candidates.
