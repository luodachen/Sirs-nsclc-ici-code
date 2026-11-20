############################################################
## 02_Step1b_strategyA_CSI_robustness.R
## Role: exposure-level descriptive analysis downstream of Step1, all based on [Strategy A: CSI uses published parameters]
##       tau = 26 years, delta = 0.7 years (Leffondre 2006 / ICARE 2013 pooled male+female values)
## Depends on: tidyverse
## Input: step1_output/smoking_index_clinical.csv (produced by 01_Step1_exposure_quantification.R)
## Outputs (all under step1_output/):
##   step1b_CSI_tertile_baseline.csv      baseline characteristics by CSI tertile (paper Table 1 structure)
##   step1b_exposure_correlation.csv      Spearman/Kendall correlation among the 5 exposure metrics
##   step1b_CSI_parameter_robustness.csv  agreement between Strategy A main parameters and other published parameters
##   step1b_tertile_counts.csv            tertile sample sizes per cohort (for Step2 group comparisons)
##   step1b_paper_snippets.md             paper methods/results phrasings (actual n filled in)
##   figs/step1b_*.pdf                    distribution / scatter / sample-size plots / correlation heatmap
## Author: see manuscript author list   Date: 2026-10-08
############################################################

## Prerequisite: first run 01_Step1_exposure_quantification.R to generate smoking_index_clinical.csv
##    This script does NOT re-download data; it only consumes the Step1 outputs and finishes in seconds.

suppressPackageStartupMessages({ library(tidyverse) })

## ---------- 0. Paths --------------------------------------------------------
proj_dir <- "D:/Projects/smoking-index-ici"          # keep consistent with script 01
out_dir  <- file.path(proj_dir, "step1_output")
fig_dir  <- file.path(out_dir, "figs")
invisible(lapply(c(out_dir, fig_dir), function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

clin_path <- file.path(out_dir, "smoking_index_clinical.csv")
if (!file.exists(clin_path)) stop("Not found: ", clin_path, "\nPlease run 01_Step1_exposure_quantification.R first")
clin <- readr::read_csv(clin_path, show_col_types = FALSE)

## ---------- 1. CSI parameters (Strategy A) and formula --------------------------------------
## Identical formula and parameters as Section 3.4 of script 01, ensuring a consistent definition (recomputed here = a self-check)
CSI_MAIN <- c(tau = 26, delta = 0.7)                 # main analysis: published parameters
CSI_SENS <- list(d10 = c(tau = 26, delta = 1.0),     # Leffondre male value
                 t20 = c(tau = 20, delta = 0.7),
                 t35 = c(tau = 35, delta = 0.7))
csi_leffondre <- function(int, dur, tsc,
                          tau = CSI_MAIN[["tau"]], delta = CSI_MAIN[["delta"]]) {
  tsc      <- ifelse(is.na(tsc), 0, tsc)             # current smokers: tsc = 0
  tsc_star <- pmax(tsc - delta, 0)
  dur_star <- pmax(dur + tsc - delta, 0) - tsc_star
  (1 - 0.5^(dur_star / tau)) * 0.5^(tsc_star / tau) * log(int + 1)
}

## Recompute the 4 parameter sets from raw fields (independent of whether Step1 exported sensitivity columns)
dat <- clin %>%
  mutate(
    csi_main = csi_leffondre(cpd, duration, years_since_quit),
    csi_d10  = csi_leffondre(cpd, duration, years_since_quit, CSI_SENS$d10[["tau"]], CSI_SENS$d10[["delta"]]),
    csi_t20  = csi_leffondre(cpd, duration, years_since_quit, CSI_SENS$t20[["tau"]], CSI_SENS$t20[["delta"]]),
    csi_t35  = csi_leffondre(cpd, duration, years_since_quit, CSI_SENS$t35[["tau"]], CSI_SENS$t35[["delta"]])
  )

## 1.1 Self-check: agreement with csi_lc exported by Step1 (reproducibility evidence for Strategy A)
if ("csi_lc" %in% names(dat)) {
  d <- abs(dat$csi_main - dat$csi_lc)
  cat("Self-check max |csi_main - Step1.csi_lc| =",
      signif(max(d, na.rm = TRUE), 6),
      " (should be 0 or on the order of 1e-15; nonzero means the two scripts differ in parameters/formula)\n")
}

## ---------- 2. Build the ever-smoker analysis set --------------------------------------
ever <- dat %>%
  filter(smoking_grp %in% c("current", "former"),
         is.finite(csi_main), !is.na(duration), !is.na(cpd))
cat("\n== Analysis set ==\n")
cat("Total ever-smokers: ", nrow(ever), " (full cohort ", nrow(dat), ")\n")
print(table(ever$cohort, ever$smoking_grp))

## ---------- 3. CSI tertile grouping (main analysis parameters) ---------------------------------
## Tertiles are cut on the pooled ever-smoker set (a common threshold across cohorts, facilitating cross-cancer comparison);
## if a reviewer requests within-cohort cutting, simply add group_by(cohort) back.
tert3 <- function(x, ref = NULL) {
  out <- rep(NA_character_, length(x))
  ok  <- is.finite(x)
  brk <- if (is.null(ref)) x[ok] else ref
  q   <- quantile(brk, probs = c(1/3, 2/3), na.rm = TRUE)
  out[ok] <- ifelse(x[ok] <= q[1], "T1_Low",
                    ifelse(x[ok] <= q[2], "T2_Mid", "T3_High"))
  out
}
ever <- ever %>% mutate(csi_t3 = tert3(csi_main))
ever$csi_t3 <- factor(ever$csi_t3, levels = c("T1_Low", "T2_Mid", "T3_High"))

## Tertile cut-offs and sample sizes
cuts <- quantile(ever$csi_main, probs = c(1/3, 2/3), na.rm = TRUE)
cat("\nCSI tertile cut-offs (all ever-smokers): ", round(cuts, 4), "\n")
tc <- ever %>% count(cohort, csi_t3, .drop = FALSE) %>% pivot_wider(names_from = csi_t3, values_from = n)
print(as.data.frame(tc))
write_csv(tc, file.path(out_dir, "step1b_tertile_counts.csv"))

## ---------- 4. Baseline characteristics table (CSI tertiles) --------------------------------------
## Paper Table 1 structure: continuous variables reported as mean(SD)/median(IQR), categorical as n(%)
stage_simpl <- function(x) {
  s <- toupper(as.character(x))
  dplyr::case_when(
    grepl("^STAGE ?IV|^IV", s)            ~ "IV",
    grepl("^STAGE ?III|^III", s)          ~ "III",
    grepl("^STAGE ?II(?!I)|^II", s, perl = TRUE) ~ "II",
    grepl("^STAGE ?I(?!I)|^I", s, perl = TRUE)   ~ "I",
    TRUE ~ NA_character_)
}
ever <- ever %>%
  mutate(stage_s = stage_simpl(stage),
         age_num = suppressWarnings(as.numeric(age_y)))

cont_stat <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("NA")
  sprintf("%.2f (%.2f) / %.2f [%.2f-%.2f]", mean(v), sd(v),
          median(v), quantile(v, .25), quantile(v, .75))
}
med_iqr <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("NA")
  sprintf("%.1f [%.1f-%.1f]", median(v), quantile(v, .25), quantile(v, .75))
}
pct <- function(lgl) {
  lgl <- lgl[!is.na(lgl)]
  if (!length(lgl)) return("NA")
  sprintf("%.1f", mean(lgl) * 100)
}
tab1 <- map_dfr(levels(ever$csi_t3), function(g) {
  d <- filter(ever, csi_t3 == g)
  tibble(
    csi_tertile = g, n = nrow(d),
    csi_mean_sd = cont_stat(d$csi_main),
    cpd_mean_sd = cont_stat(d$cpd),
    dur_mean_sd = cont_stat(d$duration),
    tsc_mean_sd = cont_stat(d$years_since_quit),
    cigyears_median = med_iqr(d$cig_years),
    py_median       = med_iqr(d$py_gdc),
    age_mean_sd = cont_stat(d$age_num),
    male_pct     = pct(tolower(d$gender) == "male"),
    LUAD_pct     = pct(d$cohort == "TCGA-LUAD"),
    current_pct  = pct(d$smoking_grp == "current"),
    stageIII_IV_pct = pct(d$stage_s %in% c("III", "IV"))
  )
})
cat("\n== Baseline characteristics by CSI tertile (Table 1 structure) ==\n"); print(as.data.frame(tab1))
write_csv(tab1, file.path(out_dir, "step1b_CSI_tertile_baseline.csv"))

## ---------- 5. Correlation among exposure metrics ---------------------------------------------
## Purpose: show that CSI is both correlated with and adds information beyond existing metrics (pack-years/smoking index/logcig-years)
met <- c(csi_main = "CSI (tau=26, delta=0.7)", csi_d10 = "CSI (26/1.0)",
         csi_t20 = "CSI (20/0.7)", csi_t35 = "CSI (35/0.7)",
         cig_years = "Cigarette-years", log_cig_years = "logcig-years",
         py_gdc = "Pack-years (GDC)", cpd = "CPD", duration = "Duration")
mtx <- ever %>% select(any_of(names(met)))
sp <- cor(mtx, method = "spearman", use = "pairwise.complete.obs")
ke <- cor(mtx, method = "kendall",  use = "pairwise.complete.obs")
cor_long <- as.data.frame(as.table(sp)) %>%
  rename(var1 = Var1, var2 = Var2, spearman = Freq) %>%
  mutate(kendall = as.vector(ke))
## Variable name -> readable label (lookup table instead of recode, to avoid !!! splicing differences across dplyr versions)
lab_of <- function(x) { out <- unname(met[as.character(x)]); ifelse(is.na(out), as.character(x), out) }
write_csv(cor_long, file.path(out_dir, "step1b_exposure_correlation.csv"))
cat("\n== Rank correlation between main-analysis CSI and other metrics ==\n")
print(cor_long %>% filter(var1 == "csi_main", var2 != "csi_main") %>%
        mutate(var2 = lab_of(var2)) %>% select(-var1) %>% as.data.frame())

## ---------- 6. Parameter robustness (key supporting material for Strategy A) ----------------------------
## No parameter re-estimation: only demonstrate that switching to other published parameters leaves the exposure ranking and tertile membership essentially unchanged
sens <- c("csi_d10", "csi_t20", "csi_t35")
sens_lab <- c(csi_d10 = "tau=26, delta=1.0 (Leffondre male)",
              csi_t20 = "tau=20, delta=0.7", csi_t35 = "tau=35, delta=0.7")
ever <- ever %>% mutate(across(all_of(sens), list(g = ~tert3(.x, ref = ever$csi_main)), .names = "{.col}_g"))
rob <- tibble(
  main_reference = paste0("CSI tau=", CSI_MAIN[["tau"]], ", delta=", CSI_MAIN[["delta"]], " (published)"),
  sensitivity    = unname(sens_lab[sens]),
  n_pair         = vapply(sens, function(p) sum(is.finite(ever$csi_main) & is.finite(ever[[p]])), integer(1)),
  spearman_rho   = vapply(sens, function(p) round(cor(ever$csi_main, ever[[p]], method = "spearman", use = "complete.obs"), 4), numeric(1)),
  kendall_tau    = vapply(sens, function(p) round(cor(ever$csi_main, ever[[p]], method = "kendall",  use = "complete.obs"), 4), numeric(1)),
  tertile_agree  = vapply(sens, function(p)
                          round(mean(as.character(ever[[paste0(p, "_g")]]) == as.character(ever$csi_t3),
                                     na.rm = TRUE) * 100, 1), numeric(1))
)
cat("\n== CSI parameter robustness (Strategy A) ==\n"); print(as.data.frame(rob))
write_csv(rob, file.path(out_dir, "step1b_CSI_parameter_robustness.csv"))
cat("Interpretation: rho >= 0.95 and tertile_agree >= 80% -> the conclusion is robust to parameter choice; may be written as\n",
    "      \"results were robust to alternative published (tau, delta) values\".\n")

## ---------- 7. Plots --------------------------------------------------------
th <- theme_bw(base_size = 11)
## 7.1 CSI distribution (cohort x smoking status)
p1 <- ggplot(filter(ever, is.finite(csi_main)), aes(smoking_grp, csi_main, fill = smoking_grp)) +
  geom_boxplot(outlier.size = .3, width = .6) +
  geom_jitter(width = .15, size = .25, alpha = .25) +
  facet_wrap(~cohort, scales = "free_y") +
  labs(title = "CSI (tau=26, delta=0.7) by smoking status",
       x = NULL, y = "Comprehensive smoking index") +
  th + theme(legend.position = "none")
ggsave(file.path(fig_dir, "step1b_CSI_by_status.pdf"), p1, width = 7.5, height = 4)

## 7.2 CSI vs pack-years (to see whether the information agrees / where it deviates)
p2 <- ggplot(filter(ever, is.finite(csi_main), is.finite(py_gdc)), aes(py_gdc, csi_main)) +
  geom_point(size = .4, alpha = .3) +
  geom_smooth(method = "loess", se = TRUE, linewidth = .8, color = "firebrick") +
  facet_wrap(~cohort, scales = "free") +
  labs(title = "CSI vs pack-years", x = "Pack-years (GDC)", y = "CSI (tau=26, delta=0.7)") + th
ggsave(file.path(fig_dir, "step1b_CSI_vs_packyears.pdf"), p2, width = 7.5, height = 4)

## 7.3 Tertile sample sizes
p3 <- ggplot(tc %>% pivot_longer(-cohort, names_to = "tertile", values_to = "n") %>%
               mutate(tertile = factor(tertile, levels = c("T1_Low", "T2_Mid", "T3_High"))),
             aes(tertile, n, fill = tertile)) +
  geom_col(width = .7) + geom_text(aes(label = n), vjust = -0.3, size = 3) +
  facet_wrap(~cohort, scales = "free_y") +
  labs(title = "CSI tertile sample sizes", x = NULL, y = "n") +
  th + theme(legend.position = "none")
ggsave(file.path(fig_dir, "step1b_tertile_counts.pdf"), p3, width = 7.5, height = 4)

## 7.4 Correlation heatmap of exposure metrics
p4 <- cor_long %>%
  mutate(v1 = lab_of(var1), v2 = lab_of(var2)) %>%
  ggplot(aes(v1, v2, fill = spearman)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", spearman)), size = 2.6) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(title = "Spearman correlation among exposure metrics", x = NULL, y = NULL) +
  th + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(fig_dir, "step1b_exposure_correlation.pdf"), p4, width = 7.5, height = 6)

## ---------- 8. Paper snippets (actual n filled in automatically) -----------------------------------
## Helper to fetch correlation values (avoids sprintf receiving length 0 when a column is missing and erroring)
get_rho <- function(a, b) {
  v <- cor_long$spearman[cor_long$var1 == a & cor_long$var2 == b]
  if (!length(v) || !is.finite(v[1])) "NA" else sprintf("%.3f", v[1])
}
rho_py  <- get_rho("csi_main", "py_gdc")
rho_cy  <- get_rho("csi_main", "cig_years")
rho_min <- if (all(!is.finite(rob$spearman_rho))) "NA" else sprintf("%.2f", min(rob$spearman_rho, na.rm = TRUE))
agr_min <- if (all(!is.finite(rob$tertile_agree)))  "NA" else sprintf("%.0f", min(rob$tertile_agree, na.rm = TRUE))

snip <- c(
  "# Strategy A (CSI using published parameters) -- paper snippet draft",
  "",
  "> Auto-generated by 02_Step1b_strategyA_CSI_robustness.R (with actual sample sizes); may be adapted directly into Methods/Results.",
  "",
  "## Methods phrasing",
  "",
  sprintf(paste0("The comprehensive smoking index (CSI) was calculated as `(1 - 0.5^(dur*/tau)) * 0.5^(tsc*/tau) * ln(int+1)`, ",
                 "with `tsc* = max(tsc - delta, 0)` and `dur* = max(dur + tsc - delta, 0) - tsc*`, where int, dur and tsc denote ",
                 "average daily cigarette number, smoking duration and time since cessation. As TCGA-LUAD/LUSC comprises lung cancer ",
                 "cases only, the outcome-specific parameters could not be re-estimated in cohort; we therefore adopted the published ",
                 "lung-cancer values (tau = %s years, delta = %s years). Analyses were restricted to ever-smokers (n = %d). ",
                 "In sensitivity analyses, CSI was recomputed using alternative published parameter values (tau/delta = 26/1.0, 20/0.7, 35/0.7)."),
          CSI_MAIN[["tau"]], CSI_MAIN[["delta"]], nrow(ever)),
  "",
  "## Results phrasing",
  "",
  sprintf("CSI tertiles (cut-offs at %.3f and %.3f) comprised %s ever-smokers.",
          cuts[1], cuts[2],
          paste(sapply(levels(ever$csi_t3), function(g) sprintf("%s: n=%d", g, sum(ever$csi_t3 == g, na.rm = TRUE))),
                collapse = ", ")),
  "",
  sprintf(paste0("CSI was strongly correlated with pack-years (Spearman rho = %s) and cigarette-years (rho = %s), ",
                 "while the alternative published parameter sets yielded nearly identical ranking ",
                 "(rho >= %s; tertile agreement >= %s%%), indicating that the exposure ordering was robust to the choice of (tau, delta)."),
          rho_py, rho_cy, rho_min, agr_min),
  "",
  "## Limitations (must be written)",
  "",
  "- (tau, delta) adopt the published lung-cancer values and were not re-estimated in this cohort; because TCGA has no healthy controls and time since cessation is frequently missing among former smokers, (tau, delta) are not separately identifiable in this cohort.",
  "- CSI relies on self-reported smoking amount, supported by QC of agreement between the GDC official pack-years and reconstructed pack-years.",
  "- Records where the quit year is later than the diagnosis year were set to NA and flagged as `quit_year_irregular`."
)
writeLines(snip, file.path(out_dir, "step1b_paper_snippets.md"))

cat("\nDone. Outputs in ", out_dir, "\n", sep = "")
