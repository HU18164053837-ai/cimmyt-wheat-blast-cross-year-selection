options(stringsAsFactors=FALSE, warn=1)
set.seed(20260812)

root <- normalizePath(getwd(), winslash="/", mustWork=TRUE)
base <- Sys.getenv("CIMMYT_ANALYSIS_ROOT", unset=normalizePath(getwd(),winslash="/",mustWork=TRUE))
hist <- file.path(base, "history")
out <- file.path(hist, "candidate_sensitivity_baseline_analysis")
tabdir <- file.path(out, "tables")
figdir <- file.path(out, "figures")
repdir <- file.path(out, "reports")
logdir <- file.path(out, "logs")
for (p in c(tabdir, figdir, repdir, logdir)) dir.create(p, recursive=TRUE, showWarnings=FALSE)
write_out <- function(x, n) write.csv(x, file.path(tabdir, n), row.names=FALSE, na="", fileEncoding="UTF-8")

portfolio <- read.csv(file.path(hist, "validation_portfolio_analysis", "tables", "01_nonredundant_validation_portfolio.csv"), check.names=FALSE)
candidates <- read.csv(file.path(hist, "evidence_pedigree_analysis", "tables", "08_candidates_with_pedigree_clusters.csv"), check.names=FALSE)
pred <- read.csv(file.path(hist, "temporal_validation_analysis", "tables", "02_all_temporal_transition_predictions.csv"), check.names=FALSE)
temporal <- read.csv(file.path(hist, "temporal_validation_analysis", "tables", "06_portfolio_temporal_validation_summary.csv"), check.names=FALSE)
adjusted <- read.csv(file.path(hist, "environment_adjusted_response_analysis", "tables", "03_environment_adjusted_candidate_summary.csv"), check.names=FALSE)
checkadv <- read.csv(file.path(hist, "reference_anchored_analysis", "tables", "05_local_check_candidate_advantage_summary.csv"), check.names=FALSE)
loyo <- read.csv(file.path(hist, "validation_portfolio_analysis", "tables", "07_leave_one_year_out_summary.csv"), check.names=FALSE)
profiles <- read.csv(file.path(hist, "validation_portfolio_analysis", "tables", "04_portfolio_location_sowing_profiles.csv"), check.names=FALSE)

stopifnot(nrow(candidates)==104, sum(portfolio$selected_for_resistance_validation)==14)
as_logical <- function(x) if (is.logical(x)) x else tolower(as.character(x))=="true"
portfolio$selected_for_resistance_validation <- as_logical(portfolio$selected_for_resistance_validation)
pred$eligible_train <- as_logical(pred$eligible_train)
pred$selected_train <- as_logical(pred$selected_train)
pred$passed_test <- as_logical(pred$passed_test)

# Exact transition-wise comparisons at the same selection size as the reference rule.
eligible <- pred[pred$eligible_train,]
transitions <- split(eligible, paste(eligible$cutoff_year, eligible$test_year, sep="->"))
baseline_rows <- lapply(transitions, function(z) {
  k <- sum(z$selected_train)
  n <- nrow(z)
  pass_n <- sum(z$passed_test)
  if (k==0) return(NULL)
  ref_pass <- sum(z$passed_test[z$selected_train])
  ord_pct <- order(z$train_median_percentile, z$train_q75_percentile, z$train_q90_disease, na.last=TRUE)
  ord_raw <- order(z$train_q90_disease, z$train_median_percentile, z$train_q75_percentile, na.last=TRUE)
  pct_pass <- sum(z$passed_test[ord_pct[seq_len(k)]])
  raw_pass <- sum(z$passed_test[ord_raw[seq_len(k)]])
  rnd_low <- qhyper(.025, pass_n, n-pass_n, k) / k
  rnd_high <- qhyper(.975, pass_n, n-pass_n, k) / k
  p_enrich <- phyper(ref_pass-1, pass_n, n-pass_n, k, lower.tail=FALSE)
  data.frame(
    cutoff_year=z$cutoff_year[1], test_year=z$test_year[1], eligible_gid=n,
    selected_gid=k, eligible_pass_rate=pass_n/n,
    reference_pass=ref_pass, reference_precision=ref_pass/k,
    simple_percentile_pass=pct_pass, simple_percentile_precision=pct_pass/k,
    simple_raw_tail_pass=raw_pass, simple_raw_tail_precision=raw_pass/k,
    random_expected_precision=pass_n/n, random_precision_low=rnd_low,
    random_precision_high=rnd_high, exact_enrichment_p=p_enrich
  )
})
baseline <- do.call(rbind, baseline_rows)
write_out(baseline, "01_transition_baseline_comparison.csv")

pooled <- data.frame(
  method=c("reference_multicriterion", "simple_historical_percentile", "simple_raw_disease_tail", "random_same_size_expected"),
  selected_events=sum(baseline$selected_gid),
  passed_events=c(sum(baseline$reference_pass), sum(baseline$simple_percentile_pass), sum(baseline$simple_raw_tail_pass), NA),
  pooled_precision=c(
    sum(baseline$reference_pass)/sum(baseline$selected_gid),
    sum(baseline$simple_percentile_pass)/sum(baseline$selected_gid),
    sum(baseline$simple_raw_tail_pass)/sum(baseline$selected_gid),
    weighted.mean(baseline$random_expected_precision, baseline$selected_gid)
  )
)
write_out(pooled, "02_pooled_baseline_comparison.csv")

# Full prespecified threshold grid. This is a sensitivity surface, not model tuning.
grid <- expand.grid(
  median_pct=c(.35,.40,.45), q75_pct=c(.40,.45,.50,.55),
  max_year_pct=c(.45,.50,.55,.60), q90_disease=c(5,10,15,20),
  worst_disease=c(10,20,30,50), KEEP.OUT.ATTRS=FALSE
)
ref_ids <- as.character(candidates$GID[
  candidates$median_percentile<=.40 & candidates$q75_percentile<=.45 &
  candidates$max_year_median_percentile<=.50 & candidates$q90_disease<=10 &
  candidates$worst_disease<=20
])
pass_mat <- matrix(FALSE, nrow(candidates), nrow(grid))
grid_rows <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  s <- grid[i,]
  pass <- candidates$median_percentile<=s$median_pct & candidates$q75_percentile<=s$q75_pct &
    candidates$max_year_median_percentile<=s$max_year_pct & candidates$q90_disease<=s$q90_disease &
    candidates$worst_disease<=s$worst_disease
  pass[is.na(pass)] <- FALSE
  pass_mat[,i] <- pass
  ids <- as.character(candidates$GID[pass])
  union_n <- length(union(ids, ref_ids))
  grid_rows[[i]] <- data.frame(
    grid[i,], selected_gid=sum(pass), selected_two_year=sum(pass & candidates$n_years==2),
    selected_three_plus=sum(pass & candidates$n_years>=3),
    reference_overlap=length(intersect(ids, ref_ids)),
    jaccard_with_reference=if (union_n) length(intersect(ids, ref_ids))/union_n else NA
  )
}
grid_out <- do.call(rbind, grid_rows)
write_out(grid_out, "03_candidate_threshold_grid.csv")

stability <- data.frame(
  GID=candidates$GID, Cross=candidates$Cross, n_years=candidates$n_years,
  robustness_tier=candidates$robustness_tier,
  threshold_selection_frequency=rowMeans(pass_mat),
  selected_reference=as.character(candidates$GID)%in%ref_ids,
  selected_conservative=rowSums(pass_mat[, grid$median_pct==.35 & grid$q75_pct==.40 & grid$max_year_pct==.45 & grid$q90_disease==5 & grid$worst_disease==10, drop=FALSE])>0,
  selected_liberal=rowSums(pass_mat[, grid$median_pct==.45 & grid$q75_pct==.55 & grid$max_year_pct==.60 & grid$q90_disease==20 & grid$worst_disease==50, drop=FALSE])>0
)
stability <- stability[order(-stability$threshold_selection_frequency, -stability$n_years),]
write_out(stability, "04_candidate_threshold_stability.csv")

# End-to-end portfolio multiverse. Each strong-rule setting is paired with a
# moderate rule that preserves the reference strong-to-moderate relaxation.
# Evidence-duration gates and pedigree de-duplication are then applied exactly
# as in the locked workflow. Three pedigree cut heights propagate redundancy
# choice through final membership rather than stopping at the robustness tier.
pairwise <- read.csv(file.path(hist, "evidence_pedigree_analysis", "tables", "07_pairwise_pedigree_similarity.csv"), check.names=FALSE)
gid_chr <- as.character(candidates$GID)
sim <- diag(1, nrow(candidates)); rownames(sim) <- colnames(sim) <- gid_chr
for (j in seq_len(nrow(pairwise))) {
  a <- match(as.character(pairwise$GID1[j]), gid_chr); b <- match(as.character(pairwise$GID2[j]), gid_chr)
  if (!is.na(a) && !is.na(b)) sim[a,b] <- sim[b,a] <- pairwise$jaccard_similarity[j]
}
hc <- hclust(as.dist(1-sim), method="average")
cluster_settings <- data.frame(nominal_similarity=c(.20,.30,.40), cut_height=c(.80,.70,.60))
full_rows <- list(); membership <- matrix(FALSE, nrow(candidates), nrow(grid)*nrow(cluster_settings)); col_i <- 1L
for (i in seq_len(nrow(grid))) {
  s <- grid[i,]
  strong_pass <- pass_mat[,i]
  moderate_pass <- candidates$median_percentile<=min(.45,s$median_pct+.05) &
    candidates$q75_percentile<=min(.60,s$q75_pct+.15) &
    candidates$max_year_median_percentile<=min(.65,s$max_year_pct+.15) &
    candidates$q90_disease<=min(20,s$q90_disease+10) &
    candidates$worst_disease<=min(50,s$worst_disease+30)
  moderate_pass[is.na(moderate_pass)] <- FALSE
  priority <- (candidates$n_years==2 & strong_pass) | (candidates$n_years>=3 & (strong_pass | moderate_pass))
  for (cc in seq_len(nrow(cluster_settings))) {
    cl <- cutree(hc, h=cluster_settings$cut_height[cc])[gid_chr]
    keep <- priority & candidates$n_years>=3
    short <- which(priority & candidates$n_years==2 & strong_pass)
    if (length(short)) {
      ord <- short[order(cl[short], candidates$q90_disease[short], candidates$worst_disease[short], candidates$median_percentile[short])]
      reps <- ord[!duplicated(cl[ord])]
      keep[reps] <- TRUE
    }
    membership[,col_i] <- keep
    full_rows[[col_i]] <- data.frame(grid_id=i, grid[i,], nominal_pedigree_similarity=cluster_settings$nominal_similarity[cc],
      strong_gid=sum(strong_pass), strong_or_moderate_gid=sum(strong_pass|moderate_pass), priority_gid=sum(priority), final_portfolio_gid=sum(keep))
    col_i <- col_i + 1L
  }
}
full_grid <- do.call(rbind, full_rows)
write_out(full_grid, "07_full_funnel_multiverse.csv")
ref_port_ids <- as.character(portfolio$GID[as_logical(portfolio$selected_for_resistance_validation)])
portfolio_stability <- data.frame(GID=candidates$GID, Cross=candidates$Cross, n_years=candidates$n_years,
  reference_portfolio=gid_chr%in%ref_port_ids, full_funnel_selection_frequency=rowMeans(membership))
portfolio_stability$stability_class <- ifelse(portfolio_stability$full_funnel_selection_frequency>=.80,"stable_core",
  ifelse(portfolio_stability$full_funnel_selection_frequency>=.50,"moderately_stable","choice_sensitive"))
portfolio_stability <- portfolio_stability[order(!portfolio_stability$reference_portfolio,-portfolio_stability$full_funnel_selection_frequency,-portfolio_stability$n_years),]
write_out(portfolio_stability, "08_full_funnel_candidate_stability.csv")

# Candidate-level breeding-use table for the 14-material prospective portfolio.
p14 <- portfolio[portfolio$selected_for_resistance_validation,]
keep <- c("GID","n_environments","adjusted_effect","adjusted_effect_ci_low","adjusted_effect_ci_high","pressure_sensitivity","sensitivity_ci_low","sensitivity_ci_high","failure_fraction","worst_environment","sensitivity_class")
p14 <- merge(p14, adjusted[,keep], by="GID", all.x=TRUE)
p14 <- merge(p14, checkadv[,c("GID","n_strata","median_delta","median_delta_ci_low","median_delta_ci_high","pct_strict_advantage")], by="GID", all.x=TRUE)
p14 <- merge(p14, temporal[,c("GID","later_years","later_pass_fraction","later_median_percentile","later_q90_disease_max")], by="GID", all.x=TRUE)
p14 <- merge(p14, loyo[,c("GID","at_least_moderate_fraction")], by="GID", all.x=TRUE)

prof_risk <- do.call(rbind, lapply(split(profiles, profiles$GID), function(z) {
  z <- z[order(-z$q75_percentile, -z$q90_disease),]
  data.frame(GID=z$GID[1], observed_location_sowing_strata=nrow(z),
             vulnerable_strata=sum(z$vulnerability_flag),
             worst_location_sowing=paste(z$location[1], z$sowing[1], sep=" | "))
}))
p14 <- merge(p14, prof_risk, by="GID", all.x=TRUE)
p14$recommended_use <- ifelse(
  p14$portfolio_role=="long_term_anchor", "temporal calibration anchor",
  ifelse(p14$portfolio_role=="intermediate_candidate", "multi-environment resistance candidate", "confirmation candidate with short evidence duration")
)
p14$risk_flag <- ifelse(
  (!is.na(p14$failure_fraction) & p14$failure_fraction>.15) |
    (!is.na(p14$at_least_moderate_fraction) & p14$at_least_moderate_fraction<.80) |
    p14$sensitivity_class=="advantage_erodes_under_pressure",
  "targeted stress testing recommended", "no major historical risk flag"
)
p14$Cross <- ifelse(is.na(p14$Cross) | p14$Cross=="", "not reported", p14$Cross)
cols <- c("GID","Cross","Sel_Hist","portfolio_role","recommended_use","risk_flag",
          "n_observations","n_years","n_locations","first_year","last_year","median_disease","q90_disease","worst_disease",
          "median_percentile","q75_percentile","n_strata","median_delta","median_delta_ci_low","median_delta_ci_high","pct_strict_advantage",
          "later_years","later_pass_fraction","later_median_percentile","later_q90_disease_max","at_least_moderate_fraction",
          "n_environments","adjusted_effect","adjusted_effect_ci_low","adjusted_effect_ci_high","pressure_sensitivity","sensitivity_ci_low","sensitivity_ci_high",
          "sensitivity_class","failure_fraction","vulnerable_strata","observed_location_sowing_strata","worst_location_sowing","worst_environment")
p14 <- p14[, cols]
p14 <- p14[order(factor(p14$portfolio_role, levels=c("long_term_anchor","intermediate_candidate","provisional_strong_unique_pedigree")), p14$median_percentile),]
write_out(p14, "05_candidate_breeding_application_table.csv")

# Compact main-text version.
compact <- p14[,c("GID","Cross","portfolio_role","n_years","n_locations","median_disease","q90_disease","median_delta","later_pass_fraction","adjusted_effect","failure_fraction","risk_flag","recommended_use")]
names(compact) <- c("GID","Name_or_cross","Evidence_role","Years","Locations","Median_disease_pct","Disease_q90_pct","Median_difference_from_LOCAL_CHECK","Later_year_pass_fraction","Background_adjusted_rank","Failure_fraction","Risk_flag","Recommended_use")
write_out(compact, "06_candidate_breeding_application_compact.csv")

# Submission-grade supplementary figure.
publication_font_family <- "Arial"
draw <- function() {
  layout(matrix(c(1,2,3),1,3))
  par(mar=c(7,4.2,2.5,.8), family="sans", cex=.78)
  bp <- barplot(pooled$pooled_precision, names.arg=c("Multicriterion","Percentile only","Raw tail only","Random expected"),
                las=2, ylim=c(0,1), col=c("#3F6B3A","#C9A227","#8B5E3C","#B8B2A7"), border=NA,
                ylab="Pooled conditional pass proportion")
  text(bp, pooled$pooled_precision, labels=sprintf("%.1f%%",100*pooled$pooled_precision), pos=3, cex=.7)
  title("a  Temporal baselines", adj=0, font.main=2, cex.main=.9)
  par(mar=c(4.2,4.2,2.5,.8))
  plot(grid_out$selected_gid, grid_out$jaccard_with_reference, pch=21, bg="#3F6B3A55", col="#3F6B3A",
       xlab="GIDs selected", ylab="Jaccard overlap with reference rule", xlim=c(0,max(grid_out$selected_gid)+2), ylim=c(0,1))
  title("b  Robustness-rule sensitivity", adj=0, font.main=2, cex.main=.9)
  z <- portfolio_stability[portfolio_stability$reference_portfolio,]
  z <- z[order(z$full_funnel_selection_frequency),]
  par(mar=c(4.2,5.2,2.5,.8))
  barplot(z$full_funnel_selection_frequency,names.arg=z$GID,horiz=TRUE,las=1,xlim=c(0,1),
          col=ifelse(z$stability_class=="stable_core","#3F6B3A",ifelse(z$stability_class=="moderately_stable","#C9A227","#8B5E3C")),
          border=NA,xlab="Selection frequency",cex.names=.58)
  abline(v=c(.5,.8),lty=3,col="#777777")
  title("c  Portfolio stability",adj=0,font.main=2,cex.main=.9)
}
export_plot <- function(fun, stem, width=183, height=92) {
  w <- width/25.4; h <- height/25.4
  if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_tiff(file.path(figdir,paste0(stem,".tiff")), width=w,height=h,units="in",res=600,compression="lzw") else tiff(file.path(figdir,paste0(stem,".tiff")),width=w,height=h,units="in",res=600,compression="lzw",family="sans")
  fun(); dev.off()
  png(file.path(figdir,paste0(stem,".png")),width=w,height=h,units="in",res=300,family="sans"); fun(); dev.off()
  if (requireNamespace("svglite", quietly=TRUE)) svglite::svglite(file.path(figdir,paste0(stem,".svg")),width=w,height=h,system_fonts=list(sans=publication_font_family)) else svg(file.path(figdir,paste0(stem,".svg")),width=w,height=h,family="sans",pointsize=8)
  fun(); dev.off()
  cairo_pdf(file.path(figdir,paste0(stem,".pdf")),width=w,height=h,family="sans"); fun(); dev.off()
}
export_plot(draw, "Figure_S1_threshold_and_baseline_sensitivity", height=110)

report <- paste0(
  "# Candidate threshold sensitivity and temporal baselines\n\n",
  "The reference multicriterion rule selected ", sum(baseline$selected_gid), " material-transition events and ", sum(baseline$reference_pass),
  " passed the next-year rule (", sprintf("%.1f",100*pooled$pooled_precision[1]), "%). Same-size percentile-only and raw-tail baselines reached ",
  sprintf("%.1f",100*pooled$pooled_precision[2]), "% and ", sprintf("%.1f",100*pooled$pooled_precision[3]), "%, respectively. The weighted random same-size expectation was ",
  sprintf("%.1f",100*pooled$pooled_precision[4]), "%. Exact transition-wise hypergeometric enrichment tests are reported without treating transition events as independent biological replicates.\n\n",
  "The threshold surface contains ", nrow(grid_out), " fixed combinations. Combined with three pedigree cut heights, the full-funnel multiverse contains ", nrow(full_grid), " complete decision paths that propagate robustness thresholds, evidence-duration gates, and pedigree de-duplication to final portfolio membership. It is used to show choice sensitivity, not to choose a post hoc optimal threshold.\n\n",
  "The breeding application table integrates evidence duration, disease tail, LOCAL CHECK contrast, temporal outcomes, background-adjusted rank, observed failure frequency, and location-by-sowing risk for all 14 prospective candidates.\n"
)
writeLines(report, file.path(repdir, "candidate_sensitivity_baseline_report.md"), useBytes=TRUE)

qa <- paste0(
  "# QA notes\n\n",
  "- Candidate threshold analysis used the locked 104 repeated-GID table.\n",
  "- Temporal baselines used the same eligible material-transition events and the same selection size within every transition.\n",
  "- Random expectations and intervals were computed exactly from the hypergeometric distribution.\n",
  "- No future-year information entered any training ranking.\n",
  "- The threshold grid was reported as a sensitivity surface and was not optimized against next-year outcomes.\n",
  "- The end-to-end multiverse propagated 768 paired strong/moderate threshold settings and three pedigree cut heights through evidence-duration gating and two-year pedigree de-duplication.\n",
  "- The candidate application table contains only the fixed 14-member resistance portfolio.\n"
)
writeLines(qa, file.path(repdir, "QA_notes.md"), useBytes=TRUE)
sink(file.path(logdir,"sessionInfo.txt")); print(sessionInfo()); sink()

summary <- data.frame(metric=c("threshold_grid_combinations","full_funnel_paths","portfolio_candidates","stable_core_reference_candidates","temporal_selected_events","temporal_passed_events","reference_conditional_pass_proportion","percentile_baseline_conditional_pass_proportion","raw_tail_baseline_conditional_pass_proportion","random_expected_conditional_pass_proportion"),
                      value=c(nrow(grid_out),nrow(full_grid),nrow(p14),sum(portfolio_stability$reference_portfolio & portfolio_stability$stability_class=="stable_core"),sum(baseline$selected_gid),sum(baseline$reference_pass),pooled$pooled_precision))
write_out(summary, "00_analysis_summary.csv")
print(summary)
