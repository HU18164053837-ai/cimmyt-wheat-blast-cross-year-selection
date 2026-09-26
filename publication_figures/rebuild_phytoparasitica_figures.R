options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
base_override <- Sys.getenv("PHYTO_CIMMYT_BASE", unset = "")
base <- if (nzchar(base_override)) base_override else root
out <- file.path(root, "publication_figures", "outputs")
src <- file.path(out, "Source_Data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(src, recursive = TRUE, showWarnings = FALSE)

font_family <- "sans" # Windows maps the sans family to TT Arial.
green <- "#315B3A"
green_light <- "#8FA98B"
ochre <- "#B8860B"
ochre_light <- "#D8BE69"
blue <- "#496D89"
grey <- "#A9A59D"
grey_dark <- "#5E625E"
ink <- "#20231F"
paper <- "#FFFFFF"

pub_par <- function(mar = c(3.2, 3.5, 0.8, 0.5), cex = 0.78) {
  par(mar = mar, mgp = c(2.0, 0.58, 0), tcl = -0.22, family = font_family,
      cex = cex, las = 1, bty = "l", lend = "round")
}

export_figure <- function(draw, stem, width_mm = 174, height_mm = 120) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  # Some Windows graphics devices cannot open Unicode output paths. Render to
  # the ASCII-safe session temp directory and then copy into the repository.
  finish_export <- function(tmp, ext) {
    dev.off()
    destination <- file.path(out, paste0(stem, ".", ext))
    if (!file.copy(tmp, destination, overwrite = TRUE)) {
      stop("Could not copy figure export to ", destination)
    }
    unlink(tmp)
  }
  tmp <- tempfile(fileext = ".svg")
  svglite::svglite(tmp, width = w, height = h, bg = paper, pointsize = 8)
  draw(); finish_export(tmp, "svg")
  tmp <- tempfile(fileext = ".pdf")
  cairo_pdf(tmp, width = w, height = h, family = font_family, pointsize = 8)
  draw(); finish_export(tmp, "pdf")
  tmp <- tempfile(fileext = ".tiff")
  tiff(tmp, width = w, height = h, units = "in", res = 600,
       compression = "lzw", family = font_family, pointsize = 8)
  draw(); finish_export(tmp, "tiff")
  tmp <- tempfile(fileext = ".png")
  png(tmp, width = w, height = h, units = "in", res = 300,
      family = font_family, pointsize = 8, type = "cairo-png", bg = paper)
  draw(); finish_export(tmp, "png")
}

panel_letter <- function(letter) {
  mtext(letter, side = 3, adj = 0, line = 0.05, font = 2, cex = 1.05)
}

# Fig. 1: workflow. The diagram is descriptive; it does not encode quantitative data.
draw_fig1 <- function() {
  par(mar = c(0, 0, 0, 0), family = font_family)
  plot.new(); plot.window(xlim = c(0, 18), ylim = c(0, 10))
  edge <- "#385044"
  fills <- c("#E8EFE5", "#D9E5D3", "#F3E5B4", "#E8D8C8", "#DAE5E8",
             "#E8EFE5", "#D9E5D3", "#F3E5B4")
  draw_box <- function(x, y, w, h, title, subtitle, fill) {
    rect(x, y, x + w, y + h, col = fill, border = edge, lwd = 1.4)
    text(x + w / 2, y + h * 0.66, paste(title, collapse = "\n"),
         cex = 0.80, font = 2, col = ink)
    text(x + w / 2, y + h * 0.23, paste(subtitle, collapse = "\n"),
         cex = 0.55, col = grey_dark)
  }
  arr <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, length = 0.08,
                                           lwd = 1.35, col = edge)
  xs <- c(0.8, 4.2, 7.6, 11.0, 14.4)
  titles <- list(c("Historical field", "nursery records"), "Data curation",
                 c("Within-nursery", "disease comparison"), "Repeated GIDs",
                 c("Candidate", "prioritization"))
  subtitles <- list("2018-2024; 24 phenotype files",
                    c("harmonized fields", "duplicate file removed"),
                    "structural absences not imputed",
                    c("same GID observed", "across years"),
                    c("14 GIDs retained", "for prospective field tests"))
  for (i in seq_along(xs)) {
    draw_box(xs[i], 6.8, 2.6, 2.0, titles[[i]], subtitles[[i]], fills[i])
    if (i > 1) arr(xs[i] - 0.35, 7.8, xs[i], 7.8)
  }
  arr(15.7, 6.8, 15.7, 5.0); arr(15.7, 5.0, 4.4, 5.0); arr(4.4, 5.0, 4.4, 3.9)
  text(9.4, 5.35, "Historical evidence for candidate prioritization",
       cex = 0.82, col = edge)
  bx <- c(2.2, 7.0, 11.8)
  bt <- list(c("Comparison with", "LOCAL CHECK"),
             c("Historical next-year", "signal"),
             c("Prospective replicated", "field trials"))
  bs <- list("shared nursery-environment strata",
             c("reobserved GIDs only", "not external validation"),
             "multi-site test with common standards")
  for (i in seq_along(bx)) {
    draw_box(bx[i], 1.5, 3.5, 2.0, bt[[i]], bs[[i]], fills[5 + i])
    if (i > 1) arr(bx[i] - 0.4, 2.5, bx[i], 2.5)
  }
  text(9, 0.55,
       "Scope: historical candidate prioritization; resistance confirmation requires prospective replicated trials",
       cex = 0.70, col = "#81582C", font = 3)
}
export_figure(draw_fig1, "Fig1_workflow", 174, 92)

# Fig. 2: environment-level disease pressure and within-nursery separation.
env_file <- file.path(base, "history", "cross_year_analysis", "tables",
                      "02_environment_pressure_and_discrimination.csv")
env <- read.csv(env_file, check.names = FALSE)
env <- env[order(env$location, env$sowing, env$year), ]
write.csv(env, file.path(src, "Fig2_environment_source_data.csv"), row.names = FALSE)
locs <- c("Jashore", "Quirusillas", "Okinawa")
sow_col <- c("1st" = green, "2nd" = ochre)
loc_pch <- c(Jashore = 21, Quirusillas = 24, Okinawa = 23)
loc_lty <- c(Jashore = 1, Quirusillas = 2, Okinawa = 3)

draw_fig2 <- function() {
  layout(matrix(1:4, 2, 2, byrow = TRUE))
  yr <- range(env$year)
  pub_par(c(3.0, 3.6, 0.9, 0.4), 0.72)
  plot(yr, range(env$balanced_median), type = "n", xaxt = "n", xlab = "Year",
       ylab = "Nursery-balanced median disease severity (%)")
  axis(1, at = seq(yr[1], yr[2]))
  for (loc in locs) for (s in c("1st", "2nd")) {
    z <- env[env$location == loc & env$sowing == s, ]; z <- z[order(z$year), ]
    if (nrow(z)) lines(z$year, z$balanced_median, type = "o", col = sow_col[s],
                       lty = loc_lty[loc], pch = loc_pch[loc], lwd = 1.25,
                       cex = 0.72, bg = paper)
  }
  panel_letter("a")
  legend("topleft", legend = c("First sowing", "Second sowing"), col = sow_col,
         lty = 1, pch = 21, pt.bg = paper, bty = "n", cex = 0.62)

  pub_par(c(3.0, 3.6, 0.9, 0.4), 0.72)
  plot(yr, range(env$discrimination_iqr), type = "n", xaxt = "n", xlab = "Year",
       ylab = "Median within-nursery IQR (%)")
  axis(1, at = seq(yr[1], yr[2]))
  for (loc in locs) for (s in c("1st", "2nd")) {
    z <- env[env$location == loc & env$sowing == s, ]; z <- z[order(z$year), ]
    if (nrow(z)) lines(z$year, z$discrimination_iqr, type = "o", col = sow_col[s],
                       lty = loc_lty[loc], pch = loc_pch[loc], lwd = 1.25,
                       cex = 0.72, bg = paper)
  }
  panel_letter("b")
  legend("topright", legend = locs, lty = loc_lty, pch = loc_pch,
         col = grey_dark, pt.bg = paper, bty = "n", cex = 0.60)

  pub_par(c(3.2, 3.6, 0.9, 0.4), 0.72)
  plot(env$balanced_median, env$discrimination_iqr,
       pch = loc_pch[env$location], bg = sow_col[env$sowing], col = ink,
       cex = 0.65 + 0.15 * env$n_nurseries,
       xlab = "Nursery-balanced median disease severity (%)",
       ylab = "Median within-nursery IQR (%)")
  panel_letter("c")
  legend("topright", legend = locs, pch = loc_pch, pt.bg = grey,
         col = ink, bty = "n", cex = 0.60)

  pub_par(c(4.5, 3.6, 0.9, 0.4), 0.72)
  dz <- env[is.finite(env$pressure_change_from_previous_year) &
              abs(env$pressure_change_from_previous_year) > 1e-12, ]
  labels <- paste0(substr(dz$location, 1, 3), "\n", dz$year, " ",
                   ifelse(dz$sowing == "1st", "S1", "S2"))
  barplot(dz$pressure_change_from_previous_year, names.arg = labels,
          col = sow_col[dz$sowing], border = NA, las = 2, cex.names = 0.55,
          ylab = "Change from preceding year (%)")
  abline(h = 0, lwd = 0.65, col = ink)
  panel_letter("d")
}
export_figure(draw_fig2, "Fig2_environment_disease_profiles", 174, 122)

# Fig. 3: evidence duration and historical response classes.
evidence_file <- file.path(base, "history", "evidence_pedigree_analysis", "tables",
                           "01_evidence_stratified_candidates.csv")
ev <- read.csv(evidence_file, check.names = FALSE)
ev$response_class <- factor(ev$robustness_tier,
                            levels = c("strong", "moderate", "not_robust"),
                            labels = c("Favourable", "Intermediate", "Inconsistent"))
ev$duration_group <- factor(ifelse(ev$n_years == 2, "2 years",
                            ifelse(ev$n_years <= 4, "3-4 years", "5-7 years")),
                            levels = c("2 years", "3-4 years", "5-7 years"))
write.csv(ev, file.path(src, "Fig3_candidate_evidence_source_data.csv"), row.names = FALSE)

draw_fig3 <- function() {
  layout(matrix(c(1, 2), 1, 2), widths = c(1.35, 0.85))
  pub_par(c(3.3, 3.8, 0.9, 0.5), 0.76)
  class_col <- c(Favourable = green, Intermediate = ochre, Inconsistent = grey)
  set.seed(20260926)
  xj <- ev$n_years + runif(nrow(ev), -0.055, 0.055)
  plot(xj, ev$q90_disease, pch = 21, bg = class_col[ev$response_class],
       col = ink, cex = pmin(1.35, 0.62 + 0.055 * ev$n_nurseries),
       xlim = c(1.7, 7.3), ylim = c(0, 102), xaxt = "n",
       xlab = "Years with observations",
       ylab = "90th-percentile disease severity (%)")
  axis(1, at = 2:7)
  abline(h = c(10, 20), lty = 3, col = grey_dark, lwd = 0.7)
  panel_letter("a")
  legend("top", legend = names(class_col), pt.bg = class_col, pch = 21,
         col = ink, bty = "n", cex = 0.63, horiz = TRUE)

  pub_par(c(3.3, 3.5, 0.9, 0.5), 0.76)
  tab <- table(ev$response_class, ev$duration_group)
  barplot(tab, col = class_col[rownames(tab)], border = NA,
          ylab = "Number of repeated GIDs", xlab = "Evidence duration",
          names.arg = c("2 years", "3-4 years", "5-7 years"),
          cex.names = 0.70)
  panel_letter("b")
}
export_figure(draw_fig3, "Fig3_candidate_evidence_duration", 174, 95)

# Fig. 4: matched disease difference against the local susceptible check.
anchor_file <- file.path(base, "history", "reference_anchored_analysis", "tables",
                         "05_local_check_candidate_advantage_summary.csv")
anc <- read.csv(anchor_file, check.names = FALSE)
anc$check_minus_candidate <- -anc$median_delta
anc$ci_low <- -anc$median_delta_ci_high
anc$ci_high <- -anc$median_delta_ci_low
anc <- anc[order(anc$check_minus_candidate), ]
write.csv(anc, file.path(src, "Fig4_local_check_source_data.csv"), row.names = FALSE)

draw_fig4 <- function() {
  pub_par(c(3.5, 4.3, 0.9, 0.5), 0.78)
  yy <- seq_len(nrow(anc))
  lim <- range(c(0, anc$ci_low, anc$ci_high), na.rm = TRUE)
  plot(lim, range(yy), type = "n", yaxt = "n",
       xlab = "LOCAL CHECK severity minus candidate severity (%)", ylab = "")
  segments(anc$ci_low, yy, anc$ci_high, yy, col = grey_dark, lwd = 1.15)
  points(anc$check_minus_candidate, yy, pch = 21, bg = green, col = ink, cex = 1.05)
  axis(2, at = yy, labels = anc$GID, las = 1, cex.axis = 0.70)
  abline(v = 0, lty = 3, col = grey_dark, lwd = 0.75)
}
export_figure(draw_fig4, "Fig4_local_check_differences", 174, 112)

# Fig. 5: transition-specific and leave-one-transition-out historical signal.
tr_file <- file.path(base, "history", "temporal_uncertainty_negative_control", "tables",
                     "01_transition_precision.csv")
loto_file <- file.path(base, "history", "temporal_uncertainty_negative_control", "tables",
                       "03_leave_one_transition_out.csv")
tr <- read.csv(tr_file, check.names = FALSE)
tr_pos <- tr[tr$selected > 0, ]
loto <- read.csv(loto_file, check.names = FALSE)
observed <- sum(tr_pos$passed) / sum(tr_pos$selected)
write.csv(tr_pos, file.path(src, "Fig5_transition_source_data.csv"), row.names = FALSE)
write.csv(loto, file.path(src, "Fig5_leave_one_transition_out_source_data.csv"), row.names = FALSE)

draw_fig5 <- function() {
  layout(matrix(c(1, 2), 1, 2), widths = c(1.0, 1.05))
  pub_par(c(6.0, 3.7, 0.9, 0.5), 0.76)
  x <- seq_len(nrow(tr_pos))
  plot(x, tr_pos$precision, ylim = c(0, 1.03), xaxt = "n", pch = 21,
       bg = ifelse(tr_pos$selected >= 10, green, ochre), col = ink, cex = 1.05,
       xlab = "", ylab = "Conditional next-year pass proportion")
  axis(1, at = x, labels = tr_pos$transition, las = 2, cex.axis = 0.62)
  mtext("Year transition", side = 1, line = 4.7)
  abline(h = observed, col = grey_dark, lwd = 1.1, lty = 2)
  label_y <- ifelse(tr_pos$precision >= 0.98,
                    tr_pos$precision - 0.08,
                    ifelse(tr_pos$precision >= 0.93,
                           tr_pos$precision - 0.13,
                           tr_pos$precision + 0.08))
  text(x, label_y, paste0("n = ", tr_pos$selected), cex = 0.58)
  panel_letter("a")

  pub_par(c(3.4, 4.3, 0.9, 0.5), 0.76)
  yy <- rev(seq_len(nrow(loto)))
  plot(loto$precision, yy, xlim = c(0.68, 0.97), yaxt = "n", pch = 21,
       bg = blue, col = ink, cex = 1.05,
       xlab = "Pass proportion after omitting one transition", ylab = "")
  axis(2, at = yy, labels = loto$omitted_transition, las = 1, cex.axis = 0.66)
  abline(v = observed, col = grey_dark, lwd = 1.1, lty = 2)
  panel_letter("b")
}
export_figure(draw_fig5, "Fig5_historical_next_year_signal", 174, 78)

# Supplementary Fig. S5: permutation reference formerly shown in Fig. 5c.
perm_file <- file.path(base, "history", "temporal_uncertainty_negative_control", "tables",
                       "05_permutation_draws.csv")
perm <- read.csv(perm_file, check.names = FALSE)
write.csv(perm, file.path(src, "FigS5_permutation_source_data.csv"), row.names = FALSE)
draw_figs1 <- function() {
  pub_par(c(3.3, 3.5, 0.8, 0.5), 0.82)
  hist(perm$precision, breaks = 35, freq = FALSE, col = "#DDE8D8", border = paper,
       xlim = range(c(perm$precision, observed)), xlab = "Pooled pass proportion",
       ylab = "Density", main = "")
  abline(v = observed, col = "#81582C", lwd = 1.7)
  abline(v = mean(perm$precision), col = blue, lwd = 1.4, lty = 2)
  legend("topleft", c(sprintf("Observed %.3f", observed),
                      sprintf("Permutation mean %.3f", mean(perm$precision))),
         col = c("#81582C", blue), lty = c(1, 2), lwd = 1.4,
         bty = "n", cex = 0.66)
}
export_figure(draw_figs1, "FigS5_transition_matched_permutation", 90, 72)

# Supplementary Fig. S6: background-relative candidate ranks.
adj_file <- file.path(base, "history", "environment_adjusted_response_analysis", "tables",
                      "03_environment_adjusted_candidate_summary.csv")
adj <- read.csv(adj_file, check.names = FALSE)
adj$background_minus_candidate <- -adj$adjusted_effect
adj$ci_low <- -adj$adjusted_effect_ci_high
adj$ci_high <- -adj$adjusted_effect_ci_low
adj <- adj[order(adj$background_minus_candidate), ]
write.csv(adj, file.path(src, "FigS6_background_adjusted_source_data.csv"), row.names = FALSE)
draw_figs2 <- function() {
  pub_par(c(3.5, 4.3, 0.8, 0.5), 0.78)
  yy <- seq_len(nrow(adj))
  plot(range(c(0, adj$ci_low, adj$ci_high), na.rm = TRUE), range(yy),
       type = "n", yaxt = "n",
       xlab = "Nursery background rank minus candidate rank", ylab = "")
  segments(adj$ci_low, yy, adj$ci_high, yy, col = grey_dark, lwd = 1.1)
  points(adj$background_minus_candidate, yy, pch = 21, bg = green, col = ink, cex = 1.0)
  axis(2, at = yy, labels = adj$GID, las = 1, cex.axis = 0.70)
  abline(v = 0, lty = 3, col = grey_dark, lwd = 0.75)
}
export_figure(draw_figs2, "FigS6_background_adjusted_candidates", 174, 110)

qa <- c(
  "Phytoparasitica figure revision QA",
  "Backend: R only.",
  "Main figures: five; supplementary figures: two.",
  "No analytical thresholds, observations, confidence intervals, or candidate membership were changed.",
  "Fig. 2 contains all 38 location-year-sowing environments.",
  "Fig. 3 contains all 104 cross-year-eligible GIDs.",
  "Fig. 4 contains all 14 candidate GIDs; exact shared-stratum counts remain in the source-data table.",
  "Fig. 5 contains all five informative transitions and all leave-one-transition-out results.",
  "Primary exports: editable SVG and PDF; submission raster: 600-dpi TIFF; PNG is for manuscript embedding and QA.",
  "Palette uses dark green, ochre, blue and neutral grey with shape/line redundancy."
)
writeLines(qa, file.path(out, "Figure_QA_Notes.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(out, "R_sessionInfo.txt"))
