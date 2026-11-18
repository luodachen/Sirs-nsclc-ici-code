############################################################
## 05_Step1c_missingness_diagnosis_and_tsc.R
## Project: smoking exposure quantification x NSCLC immune microenvironment x
##          ICI response prediction
##
## Role: a sample-size and convention check before finalizing Step1, answering
##   two questions:
##   (1) in which field and which subgroup do ever-smoker samples drop out?
##       (attrition funnel)
##   (2) when former smokers' time since cessation (tsc) is missing, does
##       "treat missing as 0" (the upper-bound convention) systematically
##       overestimate CSI? The three conventions are compared and the choice is
##       written to step1c_tsc_decision.csv as a reviewable rule.
##
## Data-version guard
##   The master table must use the CPD definition of script 02
##   (CPD = pack-years x 20 / smoking years; ever-smoker median ~20). A table
##   with median CPD < 5 indicates a stale table built from the deprecated Xena
##   cigarettes_per_day column, and the script stops and prompts to re-run 02.
##   Note: cig_years == pack-years x 20 holds by construction and therefore
##   cannot be used to identify the version; only the CPD magnitude and the
##   cpd_src column can.
##
## Dependencies: tidyverse
## Input:  step1_output/smoking_index_clinical.csv   (output of script 02)
## Output: step1_output/step1c_*.csv, step1_output/figs/step1c_*.pdf
## Author: see manuscript author list   Date: 2026-10-09
############################################################

suppressPackageStartupMessages({ library(tidyverse) })

## ---------- 0. Paths and parameters -----------------------------------------
proj_dir <- "D:/Projects/smoking-index-ici"          # consistent with 02/04/07
out_dir  <- file.path(proj_dir, "step1_output")
fig_dir  <- file.path(out_dir, "figs")
invisible(lapply(c(out_dir, fig_dir), function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))

TAU_MAIN   <- 26           # consistent with CSI_MAIN in script 02 (published strategy-A parameter)
DELTA_MAIN <- 0.7

clin <- readr::read_csv(file.path(out_dir, "smoking_index_clinical.csv"), show_col_types = FALSE)
cat("Patients read:", nrow(clin), "\n")
need_cols <- c("smoking_grp", "cohort", "cpd", "duration", "years_since_quit",
               "csi_lc", "py_gdc", "cig_years")
miss_cols <- setdiff(need_cols, names(clin))
if (length(miss_cols)) {
  stop("Master table is missing columns: ", paste(miss_cols, collapse = ", "),
       "\n  -> run script 02 first (it outputs the cpd / cpd_src / csi_lc columns)")
}

## ---------- 0.1 Data-version guard -------------------------------------------
## The master table must use the CPD definition of script 02 (CPD = pack-years x 20 / years); ever-smoker median CPD should be ~20.
## If the median < 5, a stale table built from the deprecated Xena cigarettes_per_day column was read; abort immediately.
ev <- clin %>% filter(smoking_grp %in% c("current", "former"))
cpd_med <- suppressWarnings(median(ev$cpd, na.rm = TRUE))
has_src <- "cpd_src" %in% names(clin)

cat("\n================ Data-version guard ================\n")
cat("Median cpd of ever-smokers =", round(cpd_med, 3), "(expected ~20)\n")
if (has_src) { cat("cpd source composition:\n"); print(table(ev$cpd_src, useNA = "ifany")) }
if (is.finite(cpd_med) && cpd_med < 5) {
  stop(sprintf(paste0("The Step1 master table appears to use a stale CPD definition (median CPD = %.2f; expected ~20).\n",
                      "  Cause: the table was built from the deprecated Xena cigarettes_per_day column, which is spurious (= pack-years/18.25).\n",
                      "  Fix: re-run script 02 (it overwrites the old CSV), then run this script.\n",
                      "  Note: see docs/06_data_source_audit_smoking_fields.md"), cpd_med), call. = FALSE)
} else {
  cat("Guard passed: master table uses the CPD definition of script 02 (CPD = pack-years x 20 / years).\n")
}
## Note: do not use `cig_years == pack-years x 20` to identify the version -- it is a **constructive identity**
## (CPD is back-computed from pack-years and years); the version can only be identified from the CPD magnitude and cpd_src.
rm(need_cols, miss_cols)

## ---------- 1. Missingness funnel: stepwise localization of sample loss ------
funnel <- tibble(
  step = c("1_all ever-smoker",
           "2_and cpd non-missing",
           "3_and duration non-missing",
           "4_and duration > 0",
           "5_and csi_lc finite (= eligible for analysis)"),
  n = c(nrow(ev),
        sum(!is.na(ev$cpd)),
        sum(!is.na(ev$cpd) & !is.na(ev$duration)),
        sum(!is.na(ev$cpd) & !is.na(ev$duration) & ev$duration > 0),
        sum(is.finite(ev$csi_lc)))
)
funnel$lost_from_prev <- c(NA_integer_, -diff(funnel$n))
cat("\n== Sample missingness funnel (ever-smoker) ==\n"); print(as.data.frame(funnel))
write_csv(funnel, file.path(out_dir, "step1c_missing_funnel.csv"))

lost_detail <- ev %>%
  mutate(
    has_cpd  = !is.na(cpd),
    has_dur  = !is.na(duration) & duration > 0,
    usable   = is.finite(csi_lc),
    miss_reason = case_when(
      usable                 ~ "usable",
      !has_cpd               ~ "missing CPD",
      !has_dur               ~ "missing smoking duration",
      TRUE                   ~ "other (non-finite value)")
  ) %>%
  count(cohort, smoking_grp, miss_reason) %>%
  pivot_wider(names_from = miss_reason, values_from = n, values_fill = 0)
cat("\n== Missing reason x cohort x smoking status ==\n"); print(as.data.frame(lost_detail))
write_csv(lost_detail, file.path(out_dir, "step1c_missing_reason.csv"))

## ---------- 2. tsc (time since cessation) missingness diagnosis --------------
tsc_diag <- ev %>%
  mutate(is_former = smoking_grp == "former",
         tsc_missing = is_former & is.na(years_since_quit)) %>%
  group_by(cohort, smoking_grp) %>%
  summarise(n = n(),
            tsc_missing_n = sum(tsc_missing),
            tsc_missing_pct = round(mean(tsc_missing) * 100, 1),
            tsc_median_known = round(median(years_since_quit, na.rm = TRUE), 1),
            .groups = "drop")
cat("\n== tsc (time since cessation) missingness diagnosis -- determines whether CSI is overestimated ==\n"); print(as.data.frame(tsc_diag))
write_csv(tsc_diag, file.path(out_dir, "step1c_tsc_missing_diag.csv"))

worst <- suppressWarnings(max(tsc_diag$tsc_missing_pct[tsc_diag$smoking_grp == "former"], na.rm = TRUE))
cat(sprintf("\nHighest tsc missing rate among former smokers is %.1f%%\n", worst))

## ---------- 3. CSI under three tsc conventions ------------------------------
csi_lc <- function(int, dur, tsc, tau = TAU_MAIN, delta = DELTA_MAIN, na_zero = TRUE) {
  if (na_zero) tsc <- ifelse(is.na(tsc), 0, tsc)   # current smokers tsc = 0 (correct)
  tsc_star <- pmax(tsc - delta, 0)
  dur_star <- pmax(dur + tsc - delta, 0) - tsc_star
  (1 - 0.5^(dur_star / tau)) * 0.5^(tsc_star / tau) * log(int + 1)
}

med_tsc <- ev %>%
  filter(smoking_grp == "former", !is.na(years_since_quit)) %>%
  group_by(cohort) %>%
  summarise(med_tsc = median(years_since_quit), .groups = "drop")
cat("\nCompromise imputation value (median of known tsc among former smokers, by cohort):\n"); print(as.data.frame(med_tsc))

ev2 <- ev %>%
  left_join(med_tsc, by = "cohort") %>%
  mutate(
    tsc_status = case_when(
      smoking_grp == "current"      ~ "current (tsc = 0, true value)",
      !is.na(years_since_quit)      ~ "former, tsc known",
      TRUE                          ~ "former, tsc missing"),
    tsc_upper  = ifelse(is.na(years_since_quit), 0, years_since_quit),          # (1) upper bound (current practice)
    tsc_strict = ifelse(smoking_grp == "current", 0,
                        ifelse(is.na(years_since_quit), NA_real_, years_since_quit)),  # (2) strict
    tsc_mid    = ifelse(is.na(years_since_quit), coalesce(med_tsc, 0), years_since_quit),   # (3) mid/compromise
    csi_upper  = csi_lc(cpd, duration, tsc_upper),
    csi_strict = csi_lc(cpd, duration, tsc_strict, na_zero = FALSE),
    csi_mid    = csi_lc(cpd, duration, tsc_mid),
    ## Self-check: the master table csi_lc should equal the local recomputation (upper convention); the max absolute deviation should be 0
    csi_selfcheck_diff = abs(csi_lc - csi_lc(cpd, duration, tsc_upper))
  )
cat(sprintf("\nSelf-check: max absolute deviation between the master table csi_lc and the local recomputation = %.3e (should be 0)\n",
            max(ev2$csi_selfcheck_diff, na.rm = TRUE)))

## ---------- 4. Comparison of the three conventions --------------------------
scen <- c(upper = "Upper: missing tsc as 0 (current)",
          strict = "Strict: missing tsc dropped",
          mid = "Mid: missing tsc imputed with the cohort median")
summ <- map_dfr(names(scen), function(s) {
  v <- ev2[[paste0("csi_", s)]]
  tibble(scenario = scen[[s]],
         n_usable = sum(is.finite(v)),
         pct_usable = round(mean(is.finite(v)) * 100, 1),
         median = round(suppressWarnings(median(v, na.rm = TRUE)), 4),
         q1 = round(suppressWarnings(quantile(v, .25, na.rm = TRUE)), 4),
         q3 = round(suppressWarnings(quantile(v, .75, na.rm = TRUE)), 4),
         mean = round(suppressWarnings(mean(v, na.rm = TRUE)), 4))
})
cat("\n== CSI distribution under the three tsc conventions ==\n"); print(as.data.frame(summ))
write_csv(summ, file.path(out_dir, "step1c_tsc_scenarios.csv"))

tert3 <- function(x) {
  out <- rep(NA_character_, length(x)); ok <- is.finite(x)
  out[ok] <- as.character(cut(dplyr::ntile(x[ok], 3), breaks = 0:3,
                              labels = c("T1_Low", "T2_Mid", "T3_High")))
  out
}
## Tertile columns are named uniformly *_t3 (consistent with the 02 master table and the any_of() read convention in 07)
ev2 <- ev2 %>% mutate(across(c(csi_upper, csi_strict, csi_mid),
                             tert3, .names = "{.col}_t3"))

pair_cmp <- function(a, b, na_lab) {
  ia <- is.finite(a); ib <- is.finite(b); ok <- ia & ib
  tibble(n_pair = sum(ok),
         spearman_rho = round(suppressWarnings(cor(a[ok], b[ok], method = "spearman")), 4),
         kendall_tau  = round(suppressWarnings(cor(a[ok], b[ok], method = "kendall")), 4),
         tertile_agree = round(mean(as.character(tert3(a[ok])) == as.character(tert3(b[ok]))) * 100, 1),
         median_a = round(suppressWarnings(median(a[ok])), 4),
         median_b = round(suppressWarnings(median(b[ok])), 4),
         ratio_median = round(suppressWarnings(median(a[ok])) / suppressWarnings(median(b[ok])), 4))
}
cmp <- bind_rows(
  pair_cmp(ev2$csi_strict, ev2$csi_upper, "strict vs upper") %>% mutate(comparison = "Strict vs Upper (convention)"),
  pair_cmp(ev2$csi_mid,    ev2$csi_upper, "mid vs upper")    %>% mutate(comparison = "Mid vs Upper (convention)")
) %>% select(comparison, everything())
cat("\n== Convention agreement (see the comparison column for reference) ==\n"); print(as.data.frame(cmp))
write_csv(cmp, file.path(out_dir, "step1c_tsc_scenario_agreement.csv"))

## ---------- 5. Decision table (expressing "which convention to choose" as auditable rules) ----
rule_tsc <- dplyr::case_when(
  worst > 20 ~ "use the mid convention for the primary analysis (missing rate > 20%; the upper convention clearly overestimates)",
  cmp$tertile_agree[cmp$comparison == "Mid vs Upper (convention)"] < 95 ~
    "use the mid convention for the primary analysis (missing rate is low, but the tertile agreement between conventions is < 95%)",
  TRUE ~ "use the upper convention for the primary analysis (csi_lc), with the mid convention for sensitivity"
)
decision <- tibble(
  Item = c("tsc convention", "CPD convention", "CPD definition (for Methods)", "Incomparable metrics", "Self-check"),
  Decision = c(rule_tsc,
               "CPD = pack-years x 20 / years (single definition)",
               "CPD = pack-years x 20 / smoking duration (because the official GDC field has 0% coverage)",
               "cig_years vs pack-years: CPD cannot be observed independently in TCGA, so cig_years = pack-years x 20 (constructive) and the two cannot be compared head-to-head",
               sprintf("max absolute deviation between the master table csi_lc and the local recomputation = %.3e (should be 0)",
                       max(ev2$csi_selfcheck_diff, na.rm = TRUE))),
  Basis = c(sprintf("former highest tsc missing rate %.1f%%; tertile agreement mid vs upper %.1f%%",
                   worst, cmp$tertile_agree[cmp$comparison == "Mid vs Upper (convention)"]),
           "the Xena cigarettes_per_day column is not used (it is a rescaling of pack-years)",
           "that Xena column == pack-years/18.25 (782/782); if it were cigarettes/day we would expect pack-years ~ (c/20) x years, but the observed median deviates by about 9.9-fold",
           "audit document 06, Section 2; downstream uses pack-years / logcig-years / CSI, three axes with different information structure",
           "Section 3 self-check column csi_selfcheck_diff")
)
cat("\n================ Decision table ================\n"); print(as.data.frame(decision))
write_csv(decision, file.path(out_dir, "step1c_tsc_decision.csv"))

## ---------- 6. Per-sample export (for convention sensitivity in 07 Step2) -----
by_sample <- ev2 %>%
  transmute(sample15, patient, cohort, smoking_grp,
            cpd, duration, years_since_quit, quit_year_irregular,
            csi_upper, csi_upper_t3,
            csi_strict, csi_strict_t3,
            csi_mid, csi_mid_t3)
write_csv(by_sample, file.path(out_dir, "step1c_csi_by_sample.csv"))
cat("\nPer-sample convention columns exported -> step1c_csi_by_sample.csv (script 07 merges them automatically and runs endpoint sensitivity)\n")

## ---------- 7. Figures ------------------------------------------------------
set_plot_device <- function() {
  ok <- tryCatch({
    grDevices::cairo_pdf(tempfile(fileext = ".pdf")); TRUE
  }, error = function(e) FALSE)
  if (ok) {
    options(ggplot2.device = grDevices::cairo_pdf)   # ggsave default device
    cat("[plot] using cairo_pdf (renders Chinese)\n")
  } else {
    cat("[plot] no cairo on this machine, using pdf(); all in-plot labels are already English, so no more mbcsToSbcs warnings\n")
  }
}
set_plot_device()

th <- theme_bw(base_size = 11)
p1 <- ggplot(funnel, aes(factor(step, levels = step), n)) +
  geom_col(fill = "steelblue", width = .6) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2) +
  labs(title = "Ever-smoker sample funnel", x = NULL, y = "n") +
  th + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "step1c_missing_funnel.pdf"), p1, width = 7, height = 4)

p2 <- ev2 %>% count(cohort, smoking_grp, tsc_status) %>%
  ggplot(aes(smoking_grp, n, fill = tsc_status)) +
  geom_col(position = "stack", width = .6) +
  facet_wrap(~cohort) +
  labs(title = "tsc (time since cessation) availability", x = NULL, y = "n", fill = NULL) + th
ggsave(file.path(fig_dir, "step1c_tsc_status.pdf"), p2, width = 8, height = 4)

p3 <- ev2 %>%
  filter(smoking_grp == "former") %>%
  select(cohort, csi_upper, csi_strict, csi_mid) %>%
  pivot_longer(-cohort) %>%
  filter(is.finite(value)) %>%
  ggplot(aes(name, value, fill = name)) +
  geom_boxplot(outlier.size = .3) +
  facet_wrap(~cohort, scales = "free_y") +
  labs(title = "CSI under three tsc conventions (former smokers)",
       x = NULL, y = sprintf("CSI (tau=%s, delta=%s)", TAU_MAIN, DELTA_MAIN)) +
  th + theme(legend.position = "none")
ggsave(file.path(fig_dir, "step1c_tsc_conventions.pdf"), p3, width = 8, height = 4)

## ---------- 8. Conclusions --------------------------------------------------
cat("\n================ Conclusions ================\n")
cat("(1) Sample loss: all ever-smokers", funnel$n[1], "-> analyzable", tail(funnel$n, 1),
    sprintf("(%.1f%% retained); main loss field: %s\n",
            tail(funnel$n, 1) / funnel$n[1] * 100,
            funnel$step[which.max(funnel$lost_from_prev)]))
cat("(2) tsc missing: highest among former = ", worst, "%\n", sep = "")
cat("(3) Convention decision: ", rule_tsc, "\n", sep = "")
cat("(4) CPD convention: CPD = pack-years x 20 / years (single definition).\n")
cat("\nOutputs: step1c_*.csv / figs/step1c_*.pdf\n")
