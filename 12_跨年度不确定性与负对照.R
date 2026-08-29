options(stringsAsFactors = FALSE)
set.seed(20260812)

base <- Sys.getenv("CIMMYT_ANALYSIS_ROOT", unset=normalizePath(getwd(),winslash="/",mustWork=TRUE))
infile <- file.path(base,"history","temporal_validation_analysis","tables","02_all_temporal_transition_predictions.csv")
out <- file.path(base,"history","temporal_uncertainty_negative_control")
dir.create(file.path(out, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "reports"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "logs"), recursive = TRUE, showWarnings = FALSE)

d <- read.csv(infile, check.names = FALSE)
d$transition <- paste0(d$cutoff_year, "-", d$test_year)
e <- d[d$eligible_train %in% TRUE, ]
s <- e[e$selected_train %in% TRUE, ]
stopifnot(nrow(s) > 0, all(!is.na(s$passed_test)))

tr <- aggregate(cbind(selected = as.integer(e$selected_train),
                      passed = as.integer(e$selected_train & e$passed_test),
                      eligible = rep(1L, nrow(e))),
                list(transition = e$transition), sum)
tr$precision <- ifelse(tr$selected > 0, tr$passed / tr$selected, NA_real_)
tr_pos <- tr[tr$selected > 0, ]
observed <- sum(tr_pos$passed) / sum(tr_pos$selected)
macro <- mean(tr_pos$precision)

B <- 10000L
# Transition-cluster bootstrap: resample entire selected transitions.
boot_transition <- replicate(B, {
  z <- tr_pos[sample.int(nrow(tr_pos), nrow(tr_pos), replace = TRUE), ]
  sum(z$passed) / sum(z$selected)
})
# GID-cluster bootstrap: resample GIDs and retain all selected events for each draw.
gid_split <- split(s$passed_test, s$GID)
gids <- names(gid_split)
boot_gid <- replicate(B, {
  z <- sample(gids, length(gids), replace = TRUE)
  mean(unlist(gid_split[z], use.names = FALSE))
})

boot_summary <- data.frame(
  estimand = c("Event-weighted precision", "Event-weighted precision"),
  cluster = c("Transition", "GID"), B = B, estimate = observed,
  lower_95 = c(quantile(boot_transition, .025), quantile(boot_gid, .025)),
  upper_95 = c(quantile(boot_transition, .975), quantile(boot_gid, .975))
)

loto <- do.call(rbind, lapply(tr_pos$transition, function(x) {
  z <- tr_pos[tr_pos$transition != x, ]
  data.frame(omitted_transition = x, selected = sum(z$selected),
             passed = sum(z$passed), precision = sum(z$passed) / sum(z$selected))
}))

# Negative control: permute next-year pass labels within each transition,
# preserving transition sizes, pass prevalence, and the historical selection set.
perm <- numeric(B)
for (b in seq_len(B)) {
  num <- 0L; den <- 0L
  for (tt in unique(e$transition)) {
    z <- e[e$transition == tt, ]
    y <- sample(z$passed_test, nrow(z), replace = FALSE)
    pick <- z$selected_train %in% TRUE
    num <- num + sum(y[pick])
    den <- den + sum(pick)
  }
  perm[b] <- num / den
}
perm_p <- (1 + sum(perm >= observed)) / (B + 1)
perm_summary <- data.frame(B = B, observed_precision = observed,
                           null_mean = mean(perm), null_sd = sd(perm),
                           null_lower_95 = quantile(perm, .025),
                           null_upper_95 = quantile(perm, .975),
                           one_sided_p = perm_p)

summary_tbl <- data.frame(
  selected_events = nrow(s), passed_events = sum(s$passed_test),
  unique_selected_GIDs = length(unique(s$GID)),
  informative_transitions = nrow(tr_pos),
  event_weighted_precision = observed,
  transition_macro_precision = macro,
  loto_min = min(loto$precision), loto_max = max(loto$precision),
  transition_cluster_ci_low = quantile(boot_transition, .025),
  transition_cluster_ci_high = quantile(boot_transition, .975),
  gid_cluster_ci_low = quantile(boot_gid, .025),
  gid_cluster_ci_high = quantile(boot_gid, .975),
  permutation_null_mean = mean(perm), permutation_p = perm_p
)

write.csv(summary_tbl, file.path(out,"tables","00_analysis_summary.csv"), row.names=FALSE)
write.csv(tr, file.path(out,"tables","01_transition_precision.csv"), row.names=FALSE)
write.csv(boot_summary, file.path(out,"tables","02_cluster_bootstrap_summary.csv"), row.names=FALSE)
write.csv(loto, file.path(out,"tables","03_leave_one_transition_out.csv"), row.names=FALSE)
write.csv(perm_summary, file.path(out,"tables","04_permutation_summary.csv"), row.names=FALSE)
write.csv(data.frame(draw=seq_len(B), precision=perm), file.path(out,"tables","05_permutation_draws.csv"), row.names=FALSE)

cols <- c(navy="#4F718C", blue="#3F6B3A", orange="#C9A227", gray="#B8B2A7", red="#8B5E3C")
publication_font_family <- "Arial"
draw_fig <- function() {
  par(mfrow=c(1,3), mar=c(4.2,4.2,2.1,0.8), oma=c(0,0,0.5,0), family="sans", las=1)
  x <- seq_len(nrow(tr_pos))
  point_cex <- 1 + 2.4*sqrt(tr_pos$selected/max(tr_pos$selected))
  point_fill <- ifelse(tr_pos$selected>=100,cols["blue"],cols["orange"])
  plot(x, tr_pos$precision, ylim=c(0,1.03), xaxt="n", pch=21, bg=point_fill, col="white",
       cex=point_cex, xlab="Year transition", ylab="Conditional next-year pass proportion")
  axis(1, x, tr_pos$transition, las=2, cex.axis=.72)
  abline(h=observed, col=cols["orange"], lwd=2, lty=2)
  label_y <- ifelse(tr_pos$precision >= .93, tr_pos$precision-.04, tr_pos$precision+.08)
  text(x, label_y, paste0("n=",tr_pos$selected), cex=.68)
  mtext("a", side=3, adj=0, font=2, line=.3)

  y <- rev(seq_len(nrow(loto)))
  plot(loto$precision, y, xlim=c(0.68,0.98), yaxt="n", pch=21, bg=cols["navy"], col="white",
       cex=1.4, xlab="Conditional pass proportion after omission", ylab="")
  axis(2, y, loto$omitted_transition, las=1, cex.axis=.72)
  abline(v=observed, col=cols["orange"], lwd=2, lty=2)
  mtext("b", side=3, adj=0, font=2, line=.3)

  hist(perm, breaks=35, freq=FALSE, col="#DCE7D5", border="white",
       xlim=c(min(perm)-.005, observed+.008), xlab="Pooled conditional pass proportion", ylab="Density", main="")
  abline(v=observed, col=cols["red"], lwd=2.4)
  abline(v=mean(perm), col=cols["navy"], lwd=2, lty=2)
  legend("topleft", c(sprintf("Observed %.3f",observed),sprintf("Null mean %.3f",mean(perm))),
         col=c(cols["red"],cols["navy"]), lty=c(1,2), lwd=2, bty="n", cex=.72)
  mtext("c", side=3, adj=0, font=2, line=.3)
}
base <- file.path(out,"figures","Figure_temporal_uncertainty_negative_control")
png(paste0(base,".png"), width=3600, height=1350, res=300); draw_fig(); dev.off()
tiff(paste0(base,".tiff"), width=7.2, height=2.7, units="in", res=600, compression="lzw"); draw_fig(); dev.off()
if (requireNamespace("svglite", quietly=TRUE)) svglite::svglite(paste0(base,".svg"),width=7.2,height=2.7,system_fonts=list(sans="Arial")) else svg(paste0(base,".svg"),width=7.2,height=2.7,family="sans",pointsize=8); draw_fig(); dev.off()
cairo_pdf(paste0(base,".pdf"), width=7.2, height=2.7, family="sans", pointsize=8); draw_fig(); dev.off()

writeLines(c(
  "# Temporal uncertainty and negative-control analysis",
  sprintf("Observed event-weighted precision: %.3f (%d/%d).", observed, sum(s$passed_test), nrow(s)),
  sprintf("Transition macro-average: %.3f across %d informative transitions.", macro, nrow(tr_pos)),
  sprintf("Transition-cluster bootstrap 95%% interval: %.3f–%.3f.", quantile(boot_transition,.025), quantile(boot_transition,.975)),
  sprintf("GID-cluster bootstrap 95%% interval: %.3f–%.3f.", quantile(boot_gid,.025), quantile(boot_gid,.975)),
  sprintf("Leave-one-transition-out range: %.3f–%.3f.", min(loto$precision), max(loto$precision)),
  sprintf("Within-transition permutation null mean: %.3f; one-sided p = %.5f.", mean(perm), perm_p),
  "Intervals are descriptive because only five transitions produced selected candidates; they do not establish independent external validation."
), file.path(out,"reports","analysis_report.md"))
capture.output(sessionInfo(), file=file.path(out,"logs","sessionInfo.txt"))
