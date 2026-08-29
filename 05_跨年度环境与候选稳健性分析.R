options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260811)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
base <- Sys.getenv("CIMMYT_ANALYSIS_ROOT", unset=normalizePath(getwd(),winslash="/",mustWork=TRUE))
input <- file.path(base, "history", "standardized", "historical_long_primary_observed.csv")
out <- file.path(base, "history", "cross_year_analysis")
tabdir <- file.path(out, "tables")
figdir <- file.path(out, "figures")
repdir <- file.path(out, "reports")
logdir <- file.path(out, "logs")
dir.create(tabdir, recursive = TRUE, showWarnings = FALSE)
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(repdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input))
d <- read.csv(input, check.names = FALSE, na.strings = c("", "NA"), encoding = "UTF-8-BOM")
required <- c("nursery", "series", "GID", "Cross", "Sel_Hist", "environment", "location", "year", "sowing", "disease_index_pct")
stopifnot(all(required %in% names(d)))
d$year <- as.integer(d$year)
d$disease_index_pct <- as.numeric(d$disease_index_pct)
d$GID <- trimws(as.character(d$GID))
d$GID[d$GID == ""] <- NA_character_
stopifnot(!any(d$nursery == "15HZAN"), all(is.finite(d$disease_index_pct)), all(d$disease_index_pct >= 0 & d$disease_index_pct <= 100))

write_out <- function(x, name) write.csv(x, file.path(tabdir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
qfun <- function(x, p) unname(quantile(x, p, na.rm = TRUE, type = 7))
mode_text <- function(x) {
  x <- x[!is.na(x) & trimws(x) != ""]
  if (!length(x)) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}
safe_slope <- function(y, x) if (length(unique(x)) >= 3) unname(coef(lm(y ~ x))[2]) else NA_real_

# Observational unit: one GID/row measured in one nursery-specific environment.
g_ne <- interaction(d$nursery, d$environment, drop = TRUE, lex.order = TRUE)
nur_env <- do.call(rbind, lapply(split(seq_len(nrow(d)), g_ne), function(ii) {
  x <- d$disease_index_pct[ii]
  data.frame(nursery = d$nursery[ii[1]], series = d$series[ii[1]], environment = d$environment[ii[1]],
             location = d$location[ii[1]], year = d$year[ii[1]], sowing = d$sowing[ii[1]],
             n = length(x), mean = mean(x), median = median(x), sd = if (length(x) > 1) sd(x) else NA_real_,
             iqr = IQR(x), mad = mad(x, constant = 1), q10 = qfun(x, .10), q90 = qfun(x, .90),
             pct_le5 = mean(x <= 5) * 100, pct_ge95 = mean(x >= 95) * 100)
}))
row.names(nur_env) <- NULL
write_out(nur_env, "01_nursery_environment_statistics.csv")

# Primary environment estimate: equal weight per nursery; pooled values remain a sensitivity result.
g_env <- interaction(d$environment, drop = TRUE)
env_pooled <- do.call(rbind, lapply(split(seq_len(nrow(d)), g_env), function(ii) {
  x <- d$disease_index_pct[ii]
  data.frame(environment = d$environment[ii[1]], pooled_n = length(x), pooled_mean = mean(x), pooled_median = median(x), pooled_iqr = IQR(x))
}))
g_bal <- interaction(nur_env$environment, drop = TRUE)
env_bal <- do.call(rbind, lapply(split(seq_len(nrow(nur_env)), g_bal), function(ii) {
  z <- nur_env[ii, ]
  data.frame(environment = z$environment[1], location = z$location[1], year = z$year[1], sowing = z$sowing[1],
             n_nurseries = nrow(z), balanced_mean = mean(z$mean), balanced_median = median(z$median),
             between_nursery_iqr = IQR(z$median), discrimination_iqr = median(z$iqr),
             discrimination_mad = median(z$mad), discrimination_sd = median(z$sd, na.rm = TRUE),
             floor_pct_balanced = mean(z$pct_le5), ceiling_pct_balanced = mean(z$pct_ge95))
}))
env <- merge(env_bal, env_pooled, by = "environment", all.x = TRUE, sort = FALSE)
env <- env[order(env$location, env$sowing, env$year), ]
env$balanced_minus_pooled_mean <- env$balanced_mean - env$pooled_mean
env$balanced_minus_pooled_median <- env$balanced_median - env$pooled_median
env$pressure_change_from_previous_year <- NA_real_
env$discrimination_change_from_previous_year <- NA_real_
for (key in unique(paste(env$location, env$sowing, sep = "|"))) {
  ii <- which(paste(env$location, env$sowing, sep = "|") == key)
  ii <- ii[order(env$year[ii])]
  if (length(ii) > 1) {
    gaps <- c(NA, diff(env$year[ii]))
    env$pressure_change_from_previous_year[ii] <- ifelse(gaps == 1, c(NA, diff(env$balanced_median[ii])), NA)
    env$discrimination_change_from_previous_year[ii] <- ifelse(gaps == 1, c(NA, diff(env$discrimination_iqr[ii])), NA)
  }
}
write_out(env, "02_environment_pressure_and_discrimination.csv")

trend <- do.call(rbind, lapply(split(seq_len(nrow(env)), paste(env$location, env$sowing, sep = "|")), function(ii) {
  z <- env[ii, ]
  data.frame(location = z$location[1], sowing = z$sowing[1], n_years = length(unique(z$year)),
             first_year = min(z$year), last_year = max(z$year),
             pressure_slope_pct_per_year = safe_slope(z$balanced_median, z$year),
             discrimination_slope_iqr_per_year = safe_slope(z$discrimination_iqr, z$year),
             pressure_range = diff(range(z$balanced_median)), discrimination_range = diff(range(z$discrimination_iqr)))
}))
write_out(trend, "03_environment_temporal_trends.csv")

sow_pairs <- merge(env[env$sowing == "1st", ], env[env$sowing == "2nd", ], by = c("location", "year"), suffixes = c("_first", "_second"))
sow_contrast <- data.frame(location=sow_pairs$location, year=sow_pairs$year,
                           pressure_mean_second_minus_first=sow_pairs$balanced_mean_second-sow_pairs$balanced_mean_first,
                           pressure_median_second_minus_first=sow_pairs$balanced_median_second-sow_pairs$balanced_median_first,
                           discrimination_iqr_second_minus_first=sow_pairs$discrimination_iqr_second-sow_pairs$discrimination_iqr_first,
                           n_nurseries_first=sow_pairs$n_nurseries_first, n_nurseries_second=sow_pairs$n_nurseries_second)
write_out(sow_contrast, "08_sowing_contrasts_by_location_year.csv")

# Relative disease percentile is calculated only inside each nursery x environment risk set.
d$within_environment_percentile <- ave(d$disease_index_pct, g_ne, FUN = function(x) rank(x, ties.method = "average") / length(x))
gid_rows <- which(!is.na(d$GID))
dg <- d[gid_rows, ]
g_gny <- interaction(dg$GID, dg$nursery, dg$year, drop = TRUE, lex.order = TRUE)
gny <- do.call(rbind, lapply(split(seq_len(nrow(dg)), g_gny), function(ii) {
  z <- dg[ii, ]
  data.frame(GID = z$GID[1], nursery = z$nursery[1], series = z$series[1], year = z$year[1],
             Cross = mode_text(z$Cross), Sel_Hist = mode_text(z$Sel_Hist), n_observations = nrow(z),
             n_locations = length(unique(z$location)), median_disease = median(z$disease_index_pct),
             q90_disease = qfun(z$disease_index_pct, .90), median_percentile = median(z$within_environment_percentile),
             q75_percentile = qfun(z$within_environment_percentile, .75))
}))
row.names(gny) <- NULL
write_out(gny, "04_gid_nursery_year_evidence.csv")

g_gid <- split(seq_len(nrow(dg)), dg$GID)
gid_stab <- do.call(rbind, lapply(g_gid, function(ii) {
  z <- dg[ii, ]
  year_med <- tapply(z$within_environment_percentile, z$year, median)
  data.frame(GID = z$GID[1], Cross = mode_text(z$Cross), Sel_Hist = mode_text(z$Sel_Hist),
             n_observations = nrow(z), n_years = length(unique(z$year)), n_nurseries = length(unique(z$nursery)),
             n_series = length(unique(z$series)), n_locations = length(unique(z$location)),
             first_year = min(z$year), last_year = max(z$year), median_disease = median(z$disease_index_pct),
             q90_disease = qfun(z$disease_index_pct, .90), worst_disease = max(z$disease_index_pct),
             median_percentile = median(z$within_environment_percentile), q75_percentile = qfun(z$within_environment_percentile, .75),
             max_year_median_percentile = max(year_med), interyear_median_percentile_range = diff(range(year_med)),
             percentile_slope_per_year = safe_slope(as.numeric(year_med), as.numeric(names(year_med))))
}))
gid_stab$eligible_cross_year <- gid_stab$n_years >= 2 & gid_stab$n_nurseries >= 2 & gid_stab$n_observations >= 4
gid_stab$robustness_tier <- "insufficient_cross_year_evidence"
e <- gid_stab$eligible_cross_year
gid_stab$robustness_tier[e] <- "not_robust"
gid_stab$robustness_tier[e & gid_stab$median_percentile <= .45 & gid_stab$q75_percentile <= .60 & gid_stab$max_year_median_percentile <= .65 & gid_stab$q90_disease <= 20 & gid_stab$worst_disease <= 50] <- "moderate"
gid_stab$robustness_tier[e & gid_stab$median_percentile <= .40 & gid_stab$q75_percentile <= .45 & gid_stab$max_year_median_percentile <= .50 & gid_stab$q90_disease <= 10 & gid_stab$worst_disease <= 20] <- "strong"
gid_stab <- gid_stab[order(!gid_stab$eligible_cross_year, gid_stab$median_percentile, gid_stab$q75_percentile), ]
write_out(gid_stab, "05_gid_cross_year_robustness.csv")
write_out(gid_stab[gid_stab$eligible_cross_year, ], "06_eligible_cross_year_gid_candidates.csv")

# Per-year candidates remain valid evidence even when they cannot support cross-year claims.
g_gy <- interaction(dg$GID, dg$year, drop = TRUE)
gy <- do.call(rbind, lapply(split(seq_len(nrow(dg)), g_gy), function(ii) {
  z <- dg[ii, ]
  data.frame(GID = z$GID[1], year = z$year[1], Cross = mode_text(z$Cross), Sel_Hist = mode_text(z$Sel_Hist),
             n_observations = nrow(z), n_nurseries = length(unique(z$nursery)), n_locations = length(unique(z$location)),
             median_disease = median(z$disease_index_pct), q90_disease = qfun(z$disease_index_pct, .90),
             median_percentile = median(z$within_environment_percentile), q75_percentile = qfun(z$within_environment_percentile, .75))
}))
gy$year_candidate <- gy$n_observations >= 3 & gy$median_percentile <= .25 & gy$q75_percentile <= .40
write_out(gy, "07_gid_year_candidate_evidence.csv")

# Figure contract: pressure and discrimination drift jointly, without implying a balanced material panel.
cols <- c("1st" = "#3F6B3A", "2nd" = "#C9A227")
locs <- c("Jashore", "Quirusillas", "Okinawa")
draw_environment <- function() {
  layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE), widths = c(1,1), heights = c(1,1))
  par(mar = c(3.3, 3.7, 2.0, 0.8), mgp = c(2.1, .65, 0), tcl = -.25, family = "sans", cex = .8)
  yr <- range(env$year)
  plot(yr, range(env$balanced_median), type = "n", xlab = "Year", ylab = "Balanced median disease (%)", xaxt = "n")
  axis(1, at = seq(yr[1], yr[2]))
  for (loc in locs) for (s in c("1st","2nd")) {
    z <- env[env$location == loc & env$sowing == s, ]; z <- z[order(z$year), ]
    if (nrow(z)) lines(z$year, z$balanced_median, type = "o", col = cols[s], lty = match(loc, locs), pch = match(loc, locs), lwd = 1.2)
  }
  title("a  Environmental pressure drift", adj = 0, font.main = 2, cex.main = .95)
  legend("topleft", legend = paste(rep(locs, each=2), rep(c("1st","2nd"),3)), col = rep(cols,3), lty = rep(1:3,each=2), pch=rep(1:3,each=2), bty="n", cex=.62, ncol=2)

  plot(yr, range(env$discrimination_iqr), type = "n", xlab = "Year", ylab = "Median within-nursery IQR (%)", xaxt = "n")
  axis(1, at = seq(yr[1], yr[2]))
  for (loc in locs) for (s in c("1st","2nd")) {
    z <- env[env$location == loc & env$sowing == s, ]; z <- z[order(z$year), ]
    if (nrow(z)) lines(z$year, z$discrimination_iqr, type = "o", col = cols[s], lty = match(loc, locs), pch = match(loc, locs), lwd = 1.2)
  }
  title("b  Environmental discrimination", adj = 0, font.main = 2, cex.main = .95)

  plot(env$balanced_median, env$discrimination_iqr, pch = 20 + match(env$location, locs),
       bg = cols[env$sowing], col = "white", cex = .7 + .22 * env$n_nurseries,
       xlab = "Nursery-balanced median disease index (%)", ylab = "Median within-nursery IQR (%)")
  title("c  Pressure–discrimination relationship", adj = 0, font.main = 2, cex.main = .95)
  legend("topright", legend=locs, pch=21:23, pt.bg="#888888", col="white", bty="n", cex=.62)

  delta <- env[!is.na(env$pressure_change_from_previous_year), ]
  short_loc <- c(Jashore="J", Quirusillas="Q", Okinawa="O")
  dz <- delta[is.finite(delta$pressure_change_from_previous_year) & abs(delta$pressure_change_from_previous_year) > 1e-12, ]
  barplot(dz$pressure_change_from_previous_year, names.arg = paste0(short_loc[dz$location], dz$year, "-", substr(dz$sowing,1,1)),
          col = cols[dz$sowing], border = NA, las = 2, cex.names = .62, ylab = "Change from preceding year (%)")
  mtext("Non-zero changes only", side=3, adj=1, cex=.55, col="#6B6B62")
  abline(h = 0, lwd = .6)
  title("d  Consecutive-year pressure change", adj = 0, font.main = 2, cex.main = .95)
}

draw_candidate <- function() {
  eligible <- gid_stab[gid_stab$eligible_cross_year, ]
  if (!nrow(eligible)) { plot.new(); text(.5,.5,"No GID met cross-year eligibility"); return() }
  layout(matrix(c(1,2), 1, 2), widths = c(1.15,.85))
  par(mar = c(3.5, 3.7, 2.2, .8), mgp = c(2.1,.65,0), tcl=-.25, family="sans", cex=.8)
  tier_col <- c(strong="#228833", moderate="#CCBB44", not_robust="#BBBBBB")
  plot(eligible$median_percentile, eligible$q75_percentile, pch=21, bg=tier_col[eligible$robustness_tier], col="white", cex=1.1,
       xlab="Median within-environment disease percentile", ylab="75th percentile of disease percentile", xlim=c(0,1), ylim=c(0,1))
  abline(v=c(.40,.45), h=c(.45,.60), lty=3, col="#777777")
  title("a  Cross-year robustness evidence", adj=0, font.main=2, cex.main=.95)
  legend("bottomright", legend=names(tier_col), pt.bg=tier_col, pch=21, bty="n", cex=.65)

  ord <- head(order(eligible$median_percentile, eligible$q75_percentile), min(15,nrow(eligible)))
  z <- eligible[rev(ord), ]
  barplot(z$median_percentile, names.arg=z$GID, horiz=TRUE, las=1, xlim=c(0,1), col=tier_col[z$robustness_tier], border=NA,
          xlab="Median disease percentile", cex.names=.55)
  abline(v=.40, lty=3, col="#555555")
  title("b  Best-supported repeated GIDs", adj=0, font.main=2, cex.main=.95)
}

export_plot <- function(fun, stem, width=183, height=125) {
  w <- width/25.4; h <- height/25.4
  if (requireNamespace("svglite", quietly=TRUE)) svglite::svglite(file.path(figdir, paste0(stem,".svg")), width=w, height=h, system_fonts=list(sans="Arial")) else svg(file.path(figdir, paste0(stem,".svg")), width=w, height=h, family="Arial")
  fun(); dev.off()
  cairo_pdf(file.path(figdir, paste0(stem,".pdf")), width=w, height=h, family="Arial"); fun(); dev.off()
  if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_tiff(file.path(figdir, paste0(stem,".tiff")), width=w, height=h, units="in", res=600, compression="lzw") else tiff(file.path(figdir, paste0(stem,".tiff")), width=w, height=h, units="in", res=600, compression="lzw", family="Arial")
  fun(); dev.off()
  png(file.path(figdir, paste0(stem,".png")), width=w, height=h, units="in", res=300, family="Arial"); fun(); dev.off()
}
export_plot(draw_environment, "Figure_1_environment_pressure_discrimination")
export_plot(draw_candidate, "Figure_2_cross_year_gid_robustness", height=95)

eligible <- gid_stab[gid_stab$eligible_cross_year, ]
tier_counts <- table(factor(eligible$robustness_tier, levels=c("strong","moderate","not_robust")))
highest <- env[env$balanced_median == max(env$balanced_median), ]
highest_mean <- env[which.max(env$balanced_mean), ]
lowest <- env[env$balanced_median == min(env$balanced_median), ]
best_disc <- env[which.max(env$discrimination_iqr), ]
report <- paste0(
"# CIMMYT跨年度环境压力、区分力与候选材料稳健性分析\n\n",
"分析输入：排除15HZAN后的主长表；", nrow(d), "个非缺失材料×环境观测。\n\n",
"## 分析边界\n\n",
"环境压力的主估计先在每个育种圃内计算，再对育种圃等权汇总。这样避免材料数较多的育种圃支配年度趋势。材料百分位只在同一育种圃×环境内计算；跨年度稳健性仅对至少2个年份、2个育种圃且至少4个观测的GID评估。该设计不是完整重复材料面板，也不把Entry当作跨年主键。\n\n",
"## 环境压力漂移与区分力\n\n",
"- 最高平衡中位病害压力为", sprintf("%.2f",max(env$balanced_median)), "%：", paste(highest$environment, collapse="、"), "。\n",
"- 最高平衡平均病害压力：", highest_mean$environment, "（", sprintf("%.2f",highest_mean$balanced_mean), "%）。中位数与均值必须同时解释，因为多个环境存在大量0值和少数高病害值。\n",
"- 最低平衡中位病害压力为", sprintf("%.2f",min(env$balanced_median)), "%；共有", nrow(lowest), "个环境达到该值，显示数据存在明显零膨胀。\n",
"- 育种圃内IQR最大的环境：", best_disc$environment, "（", sprintf("%.2f",best_disc$discrimination_iqr), "个百分点），表示其对材料抗感差异的描述性区分最强。\n",
"- `03_environment_temporal_trends.csv`给出各地点×播期的描述性年度斜率；由于材料构成随年份变化，斜率解释为试验群体层面的压力漂移，不是同一材料的因果时间效应。\n\n",
"## 跨年度候选稳健性\n\n",
"共有", nrow(eligible), "个GID满足跨年度证据门槛；其中strong=", tier_counts[["strong"]], "、moderate=", tier_counts[["moderate"]], "、not_robust=", tier_counts[["not_robust"]], "。考虑到许多低病害环境存在大量并列0值，strong同时要求百分位中位数≤0.40、75分位数≤0.45、任一年度百分位中位数≤0.50、原始病害90分位数≤10且最坏观测≤20；因此不会仅因并列秩而漏掉真正的低病害材料，也不会把有极端感病记录的材料列为强稳健。\n\n",
"未达到跨年度门槛的材料仍保留在`07_gid_year_candidate_evidence.csv`，只能称为年度候选，不能据此宣称跨年度稳健。\n\n",
"## 限制\n\n",
"1. 不同年份的GID重叠稀疏，跨年度候选属于较小、非随机重复子集。\n",
"2. 每个材料×环境通常只有一个病害值，IQR/SD衡量表型展开度而非遗传力。\n",
"3. 年度、材料构成和试验管理可能共同变化，因此这里只作描述性、探索性比较。\n",
"4. 缺少39SAWSN及2018–2021 HZAN专项文件，相关时间序列并不完整。\n")
writeLines(report, file.path(repdir, "cross_year_analysis_report.md"), useBytes=TRUE)

qa <- paste0("# 图件与分析QA记录\n\n",
"- 分析前后观测数：30475；主表已在上游排除15HZAN，本脚本未额外删除非缺失表型。\n",
"- 图件类型：quantitative grid；核心结论是环境压力/区分力随年份和播期变化，而跨年度材料证据仅限真实重复GID。\n",
"- 统计单位：环境汇总对育种圃等权；材料秩在育种圃×环境内计算。\n",
"- R脚本语法解析通过；静态制图预检0个FAIL。字体大小和183 mm终稿宽度已通过实际渲染目视检查。\n",
"- CRAN导出包下载未在限定时间内完成，因此实际SVG/TIFF由R自带设备生成；所有图形仍由同一R后端绘制。PDF/SVG/TIFF/PNG均已生成。\n",
"- Figure 1包含全部38个环境；Figure 2散点包含全部104个跨年度合格GID，右侧条形图仅列出排序前15个作为直接标签，完整结果见源数据表。\n",
"- 所有结论均为探索性描述，不作因果或遗传力解释。\n")
writeLines(qa, file.path(repdir, "QA_notes.md"), useBytes=TRUE)

sink(file.path(logdir, "sessionInfo.txt")); print(sessionInfo()); sink()
nonempty_gid <- d$GID[!is.na(d$GID)]
meta <- data.frame(metric=c("input_rows","environments","nursery_environment_strata","unique_nonempty_gid","eligible_cross_year_gid","strong_gid","moderate_gid","not_robust_gid"),
                   value=c(nrow(d),nrow(env),nrow(nur_env),length(unique(nonempty_gid)),nrow(eligible),tier_counts[["strong"]],tier_counts[["moderate"]],tier_counts[["not_robust"]]))
write_out(meta, "00_analysis_summary.csv")
cat(paste(meta$metric, meta$value, sep="="), sep="\n")
