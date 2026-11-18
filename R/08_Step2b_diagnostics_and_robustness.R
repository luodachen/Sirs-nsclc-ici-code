############################################################
## 08_Step2b_diagnostics_and_robustness.R -- Step2 diagnostics and robustness
## Project: smoking exposure quantification x NSCLC immune microenvironment x ICI response prediction
##
## Role: the diagnostics layer between Step2 (07) and Step3 modeling. Step2 has run, but before modeling
##       the 4 questions "a reviewer will certainly ask" must be answered, otherwise modeling rests on untested associations.
##
## Background (contradictions exposed by actual Step2 results):
##   OK SBS4 weight ever vs never: 0.278 vs 0.035, Wilcoxon p = 8.0e-06 (very strong positive control)
##   X but duration x SBS4 = -0.129 (p=0.0025), log_cig_years x SBS4 = -0.108 (p=0.014)
##     -- the longer the exposure, the LOWER the smoking-signature "proportion", **counter-intuitive**
##   X and ever vs never is not significant on TIS/CYT/CD8A/CD274 (p>0.3)
##
## Core hypothesis (to be tested by this script):
##   SBS4 weight is a **proportion** (Signature.4 share of this sample's mutations) and is subject to "denominator dilution":
##   smokers have higher TMB -> larger total mutation count -> if other signatures (SBS1 clock / SBS2&13 APOBEC) grow faster,
##   then the SBS4 **proportion** is diluted. So the "negative correlation" may be an artefact of a proportional metric, not biological truth.
##   With **absolute exposure** (SBS4 weight x TMB, ~SBS4 mutations per Mb) the negative correlation should flip positive.
##
## This script answers 4 questions (across 5 analysis blocks):
##   Q1  the SBS4 "proportion vs absolute" question -> Section 3 (most critical)
##   Q2  whether histology (LUAD/LUSC) is a confounder -> Section 4 stratification + Section 6 partial correlation
##   Q3  whether immune-phenotype associations are robust (including the clumped-data pitfall) -> Section 5 multi-method triangulation
##   Q4  whether associations are confounded by clinical variables -> Section 6 partial correlation (age/sex/stage/histology)
##
## Input: step2_output/step2_immune_tme.csv (output of 07)
## Output: step2_output/step2b_*.csv, figs/step2b_*.pdf
##       + 08_Step2b_diagnostic_summary.md (quick summary of conclusions, for direct use in manuscript writing)
##         Filenames use ASCII: under a non-UTF-8 Windows locale R cannot write Chinese filenames (Permission denied)
##
## Dependencies: tidyverse, data.table (ppcor optional; if missing, a self-implemented partial correlation is used)
## Usage: setwd("D:/Projects/smoking-index-ici"); source("08_Step2b_diagnostics_and_robustness.R")
## Author: see manuscript author list   Date: 2026-10-09
############################################################

## 0. Configuration ----------------------------------------------------------
proj_dir <- "D:/Projects/smoking-index-ici"
out_dir  <- file.path(proj_dir, "step2_output"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir  <- file.path(proj_dir, "figs");         dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({ library(tidyverse); library(data.table) })

## Plot device: consistent with 05/07, use cairo_pdf uniformly (available on this machine), avoiding CJK font warnings
if (isTRUE(capabilities("cairo"))) options(ggplot2.device = grDevices::cairo_pdf)

has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAS_PPCOR <- has_pkg("ppcor")
cat("Dependencies: ppcor =", HAS_PPCOR, " (a self-implemented partial correlation is used if missing)\n\n")

## 1. Read the Step2 master table + variable preparation ----------------------
step2 <- fread(file.path(out_dir, "step2_immune_tme.csv"),
               na.strings = c("", "NA")) %>% as_tibble()
cat("Step2 master table:", nrow(step2), "samples x", ncol(step2), "columns\n")

## Pitfall 1: age_y is actually in "days" (range 12179-32872, median 67.8 years)
stopifnot(median(step2$age_y, na.rm = TRUE) > 1000)   # without this assertion, treating "days" as "years" turns age into noise
step2 <- step2 %>% mutate(age = age_y / 365.25)

## Pitfall 2: stage is a mix of "IA/IB/IIA/IIB/IIIA/IIIB/IV/I/II/III" -> must collapse to I/II/III/IV
step2 <- step2 %>%
  mutate(stage_grp = case_when(
    grepl("^Stage\\s*IV",  stage) ~ "IV",
    grepl("^Stage\\s*III", stage) ~ "III",
    grepl("^Stage\\s*II",  stage) ~ "II",
    grepl("^Stage\\s*I",   stage) ~ "I",
    TRUE ~ NA_character_))

## Pitfall 3: q_T.cells.CD4 has 787/1017 == 0 (quanTIseq floors the value to 0 in low-CD4 samples);
##   this is "clumped data" (zero-inflated) -> Spearman is distorted; must use multi-method triangulation (Section 5)
qc_cd4 <- tibble(
  metric = c("number of samples with q_T.cells.CD4 == 0", "share of non-missing (%)"),
  value  = c(sum(step2$q_T.cells.CD4 == 0, na.rm = TRUE),
             round(100 * mean(step2$q_T.cells.CD4 == 0, na.rm = TRUE), 1)))
cat("\n== Clumped-data check ==\n"); print(as.data.frame(qc_cd4))

## ever-smoker analysis set
ev <- step2 %>% filter(smoking_grp %in% c("current", "former"))

## Exposure axes and endpoint definitions
EXPOSURES <- c("py_gdc", "log_cig_years", "csi_lc", "cpd", "duration")
## Endpoints: CD4/CD8/macrophage from three deconvolution methods (for triangulation) + immune features + TMB/SBS4
ENDPOINTS <- base::intersect(
  c("tmb_per_mb", "sbs4_weight", "sbs4_abs",
    "tis_score", "cyt_score", "cd274_log2tpm", "cd8a_log2tpm",
    "q_T.cells.CD4", "q_T.cells.CD8", "q_Tregs", "q_Macrophages.M2",
    "e_CD4_Tcells", "e_CD8_Tcells", "e_Macrophages",
    "mcp_T cells", "mcp_CD8 T cells", "mcp_Cytotoxic lymphocytes"),
  names(step2))

## 2. Construct SBS4 absolute exposure (core variable of this script) ----------
##   sbs4_weight  = **proportion** of SBS4 among all mutations of this sample (0-1)
##   tmb_per_mb   = nonsynonymous mutation density (mut/Mb)
##   sbs4_abs     = SBS4 weight x TMB  ~ number of SBS4 mutations per Mb (absolute exposure intensity)
##
##   Note: strictly, the SBS4 weight is fitted on **all SNVs** while TMB counts only nonsynonymous SNV+indel,
##       so the denominators differ -> sbs4_abs is an "approximate absolute quantity". Its purpose here is **direction self-verification**
##       (whether the negative correlation is caused by proportional normalization); absolute precision is not sought; state this in Methods.
step2 <- step2 %>% mutate(sbs4_abs = sbs4_weight * tmb_per_mb)
cat(sprintf("\nConstructed sbs4_abs: SBS4 weight x TMB (approximate absolute exposure)\n  non-missing %d cases, median %.3f, IQR [%.3f, %.3f]\n",
            sum(is.finite(step2$sbs4_abs)),
            median(step2$sbs4_abs, na.rm = TRUE),
            quantile(step2$sbs4_abs, .25, na.rm = TRUE),
            quantile(step2$sbs4_abs, .75, na.rm = TRUE)))
ev <- step2 %>% filter(smoking_grp %in% c("current", "former"))

## 3. Q1: the SBS4 "proportion vs absolute" question (the most critical block) ----
##   If duration/log_cig_years is negative for sbs4_weight but turns positive for sbs4_abs ->
##   this proves the negative correlation is an **arithmetic effect of proportional normalization**, not "longer smoking means a smaller smoking signature".
cat("\n================ Q1: SBS4 proportion vs absolute ================\n")
## Note: `pivot_wider(values_from = c(rho, p_adj), names_glue=...)` with "two metric rows
##   for the same exposure" degenerates into two half-filled rows (each row fills only one metric) -- encountered in practice.
##   Here we instead **explicitly split by metric and join**, making the structure clear.
sbs_cmp <- tidyr::expand_grid(exposure = EXPOSURES, metric = c("sbs4_weight", "sbs4_abs")) %>%
  purrr::pmap_dfr(function(exposure, metric) {
    d <- ev %>% filter(is.finite(.data[[exposure]]), is.finite(.data[[metric]]))
    if (nrow(d) < 30) return(tibble())
    ct <- suppressWarnings(cor.test(d[[exposure]], d[[metric]], method = "spearman"))
    tibble(exposure = exposure, metric = metric, n = nrow(d),
           rho = unname(ct$estimate), p_value = ct$p.value)
  }) %>%
  group_by(metric) %>% mutate(p_adj = p.adjust(p_value, method = "BH")) %>% ungroup() %>%
  mutate(sig = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                         p_adj < 0.05 ~ "*", TRUE ~ ""))
## Explicit split then join (avoids the names_glue pitfall)
sbs_cmp <- sbs_cmp %>%
  filter(metric == "sbs4_weight") %>%
  select(exposure, n_weight = n, rho_weight = rho, p_adj_weight = p_adj, sig_weight = sig) %>%
  full_join(sbs_cmp %>% filter(metric == "sbs4_abs") %>%
              select(exposure, n_abs = n, rho_abs = rho, p_adj_abs = p_adj, sig_abs = sig),
            by = "exposure") %>%
  mutate(direction_flip = sign(rho_weight) != sign(rho_abs)) %>%
  select(exposure, n_weight, rho_weight, p_adj_weight, sig_weight,
         n_abs, rho_abs, p_adj_abs, sig_abs, direction_flip)
fwrite(sbs_cmp, file.path(out_dir, "step2b_q1_sbs4_proportion_vs_absolute.csv"))
print(as.data.frame(sbs_cmp %>% mutate(across(where(is.numeric), ~round(.x, 4)))))

## Mechanistic support: relationship between cumulative exposure and TMB (if positive, the "larger denominator" holds)
cat("\n-- Mechanistic support: exposure x TMB (does the denominator grow with exposure) --\n")
denom <- purrr::map_dfr(EXPOSURES, function(v) {
  d <- ev %>% filter(is.finite(.data[[v]]), is.finite(tmb_per_mb))
  if (nrow(d) < 30) return(tibble())
  ct <- suppressWarnings(cor.test(d[[v]], d$tmb_per_mb, method = "spearman"))
  tibble(exposure = v, endpoint = "tmb_per_mb", n = nrow(d),
         rho = unname(ct$estimate), p_value = ct$p.value)
})
print(as.data.frame(denom %>% mutate(rho = round(rho, 3), p_value = signif(p_value, 3))))

## 4. Q2: Stratification by histology (LUAD / LUSC) ----------------------------
##   Median ever-SBS4 is only 0.217 in LUSC vs 0.416-0.468 in LUAD -> histology is a strong confounder;
##   must verify whether the "pooled-cohort association" also holds within each histology.
cat("\n================ Q2: Histology stratification ================\n")
by_h <- tidyr::expand_grid(exposure = EXPOSURES,
                           endpoint = c("tmb_per_mb", "sbs4_weight", "sbs4_abs",
                                        "q_T.cells.CD4", "q_Tregs"),
                           cohort = c("TCGA-LUAD", "TCGA-LUSC", "Pooled")) %>%
  purrr::pmap_dfr(function(exposure, endpoint, cohort) {
    d <- if (cohort == "Pooled") ev else dplyr::filter(ev, cohort == .env$cohort)
    d <- d %>% filter(is.finite(.data[[exposure]]), is.finite(.data[[endpoint]]))
    if (nrow(d) < 30) return(tibble())
    ct <- suppressWarnings(cor.test(d[[exposure]], d[[endpoint]], method = "spearman"))
    tibble(exposure = exposure, endpoint = endpoint, cohort = cohort, n = nrow(d),
           rho = unname(ct$estimate), p_value = ct$p.value)
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                         p_adj < 0.05 ~ "*", TRUE ~ ""))
fwrite(by_h, file.path(out_dir, "step2b_q2_histology_stratified.csv"))

cat("-- Direction-consistency check (do the rho values for the same exposure-endpoint have the same sign in both cohorts) --\n")
## Note: the column names produced by pivot_wider are `rho_TCGA-LUAD` (a **hyphen**, not a dot) --
##   encountered in practice: writing rho_TCGA.LUAD gives object not found. Here an explicit join avoids guessing column names.
dir_consist <- by_h %>% filter(cohort == "TCGA-LUAD") %>%
  select(exposure, endpoint, rho_LUAD = rho, sig_LUAD = sig) %>%
  inner_join(by_h %>% filter(cohort == "TCGA-LUSC") %>%
               select(exposure, endpoint, rho_LUSC = rho, sig_LUSC = sig),
             by = c("exposure", "endpoint")) %>%
  mutate(same_direction = sign(rho_LUAD) == sign(rho_LUSC),
         both_sig = sig_LUAD != "" & sig_LUSC != "")
print(as.data.frame(dir_consist %>% mutate(across(where(is.numeric), ~round(.x, 3)))))

## 5. Q3: Immune-phenotype triangulation (breaking the "clumped data" pitfall) --
##   Problem: 77% of q_T.cells.CD4 is 0 -> Spearman only reflects the "0 vs non-0" clump.
##   Remedy: for the same exposure, compute rho against the CD4/CD8/macrophage estimates of
##         (1) quanTIseq, (2) EPIC, (3) MCP-counter; if all three agree in direction the conclusion is robust; if one is unique, methodology is suspect.
cat("\n================ Q3: Multi-method triangulation (CD4/CD8/Macrophage) ================\n")
tri_map <- tibble::tribble(
  ~cell,        ~quantiseq,        ~epic,           ~mcp,
  "CD4 T",      "q_T.cells.CD4",   "e_CD4_Tcells",  "mcp_T cells",
  "CD8 T",      "q_T.cells.CD8",   "e_CD8_Tcells",  "mcp_CD8 T cells",
  "Macrophage", "q_Macrophages.M2","e_Macrophages", "mcp_Monocytic lineage")
tri <- tri_map %>%
  tidyr::pivot_longer(-cell, names_to = "method", values_to = "col") %>%
  filter(col %in% names(step2)) %>%
  purrr::pmap_dfr(function(cell, method, col) {
    purrr::map_dfr(EXPOSURES, function(v) {
      d <- ev %>% filter(is.finite(.data[[v]]), is.finite(.data[[col]]))
      if (nrow(d) < 30) return(tibble())
      ct <- suppressWarnings(cor.test(d[[v]], d[[col]], method = "spearman"))
      tibble(cell = cell, method = method, exposure = v, n = nrow(d),
             rho = unname(ct$estimate), p_value = ct$p.value)
    })
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         sig = case_when(p_adj < 0.05 ~ "*", TRUE ~ ""))
fwrite(tri, file.path(out_dir, "step2b_q3_multimethod_triangulation.csv"))

cat("-- Triangulation consistency summary (by cell x exposure: do the three methods agree in rho direction) --\n")
tri_sum <- tri %>%
  group_by(cell, exposure) %>%
  summarise(n_method = n(),
            rho_min = round(min(rho), 3), rho_max = round(max(rho), 3),
            same_dir = length(unique(sign(rho))) == 1,
            n_sig = sum(sig != ""), .groups = "drop") %>%
  mutate(robust = same_dir & n_sig >= 2)
print(as.data.frame(tri_sum))

## 6. Q4: Partial correlation (controlling for age / sex / stage / histology) ---
##   Self-implemented partial correlation (order = number of control variables), avoiding a ppcor dependency:
##     r_partial(x, y | Z) = correlation of the residuals from regressing x and y on Z
partial_spearman <- function(x, y, Z) {
  ok <- is.finite(x) & is.finite(y) & stats::complete.cases(Z)
  x <- x[ok]; y <- y[ok]; Z <- as.data.frame(Z)[ok, , drop = FALSE]
  ## Use ranks for partial correlation (= partial Spearman): first convert x, y, Z to ranks, then correlate the residuals
  rx <- rank(x); ry <- rank(y)
  Z  <- as.data.frame(lapply(Z, function(z) rank(z)))
  if (ncol(Z) == 0) {
    r  <- stats::cor(rx, ry)
    n  <- length(rx); tt <- r * sqrt((n - 2) / (1 - r^2))
    return(tibble(rho = r, n = n, p_value = 2 * stats::pt(-abs(tt), df = n - 2)))
  }
  fit <- stats::lm(cbind(rx, ry) ~ as.matrix(Z))
  res <- stats::residuals(fit)
  r   <- stats::cor(res[, 1], res[, 2])
  n   <- length(rx); k <- ncol(Z)
  df  <- n - k - 2
  tt  <- r * sqrt(df / (1 - r^2))
  tibble(rho = r, n = n, p_value = 2 * stats::pt(-abs(tt), df = df))
}
## Fallback when missing: if ppcor is present, cross-check once
if (HAS_PPCOR) cat("(ppcor available; a cross-check is done at the end of Section 6)\n")

cat("\n================ Q4: Partial correlation (controlling for confounders) ================\n")
conf_cols <- c("age", "gender", "stage_grp", "cohort")
pp_tab <- tidyr::expand_grid(exposure = EXPOSURES,
                             endpoint = c("tmb_per_mb", "sbs4_weight", "sbs4_abs",
                                          "q_T.cells.CD4", "q_Tregs",
                                          "q_Macrophages.M2", "tis_score", "cyt_score")) %>%
  purrr::pmap_dfr(function(exposure, endpoint) {
    d <- ev %>% filter(is.finite(.data[[exposure]]), is.finite(.data[[endpoint]]),
                       !is.na(gender), !is.na(stage_grp), is.finite(age))
    if (nrow(d) < 30) return(tibble())
    raw <- suppressWarnings(cor.test(d[[exposure]], d[[endpoint]], method = "spearman"))
    m   <- d %>% mutate(gender_m = as.integer(factor(gender)) - 1L,
                        stage_n  = as.integer(factor(stage_grp, levels = c("I","II","III","IV"))),
                        cohort_n = as.integer(factor(cohort)) - 1L)
    Z <- m %>% select(gender_m, stage_n, cohort_n, age)
    pa <- partial_spearman(m[[exposure]], m[[endpoint]], Z)
    tibble(exposure = exposure, endpoint = endpoint,
           n_raw = nrow(d),
           rho_raw = unname(raw$estimate), p_raw = raw$p.value,
           rho_partial = pa$rho, p_partial = pa$p_value)
  }) %>%
  mutate(p_partial_adj = p.adjust(p_partial, method = "BH"),
         attenuated = abs(rho_partial) < abs(rho_raw),
         conclusion_changed = (p_raw < 0.05) != (p_partial < 0.05))
fwrite(pp_tab, file.path(out_dir, "step2b_q4_partial_correlation.csv"))
print(as.data.frame(pp_tab %>% mutate(across(where(is.numeric), ~signif(.x, 3)))))

## Cross-check (only when ppcor is available)
if (HAS_PPCOR) {
  cat("\n-- ppcor cross-check (one case: sbs4_abs x duration) --\n")
  d <- ev %>% filter(is.finite(duration), is.finite(sbs4_abs), !is.na(gender),
                     !is.na(stage_grp), is.finite(age)) %>%
    mutate(gender_m = as.integer(factor(gender)) - 1L,
           stage_n  = as.integer(factor(stage_grp, levels = c("I","II","III","IV"))),
           cohort_n = as.integer(factor(cohort)) - 1L)
  ref <- ppcor::pcor.test(d$duration, d$sbs4_abs, d[, c("gender_m","stage_n","cohort_n","age")],
                          method = "spearman")
  mine <- partial_spearman(d$duration, d$sbs4_abs, d[, c("gender_m","stage_n","cohort_n","age")])
  cat(sprintf("  ppcor rho=%.4f p=%.3g | self-implemented rho=%.4f p=%.3g\n",
              ref$estimate, ref$p.value, mine$rho, mine$p_value))
}

## 7. Figures ----------------------------------------------------------------
## 7.1 Core figure: SBS4 proportion vs absolute (two bars for the same exposure, negative -> positive)
if (nrow(sbs_cmp)) {
  cmp_long <- sbs_cmp %>%
    select(exposure, rho_weight, rho_abs) %>%
    tidyr::pivot_longer(-exposure, names_to = "metric", values_to = "rho") %>%
    mutate(metric = recode(metric,
                           rho_weight = "SBS4 weight (proportion)",
                           rho_abs    = "SBS4 absolute (weight x TMB)"))
  p_cmp <- ggplot(cmp_long, aes(exposure, rho, fill = metric)) +
    geom_col(position = position_dodge(.8), width = .7) +
    geom_hline(yintercept = 0, linewidth = .3) +
    labs(title = "SBS4: proportion vs absolute exposure (ever-smokers)",
         subtitle = "Flip of sign indicates the negative correlation is a normalization artefact",
         x = NULL, y = "Spearman rho", fill = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave(file.path(fig_dir, "step2b_SBS4_proportion_vs_absolute.pdf"), p_cmp, width = 7.5, height = 4.5)
}

## 7.2 Histology-stratified heatmap
p_h <- by_h %>% filter(endpoint %in% c("tmb_per_mb","sbs4_weight","sbs4_abs","q_T.cells.CD4")) %>%
  ggplot(aes(exposure, endpoint, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f%s", rho, sig)), size = 2.8) +
  facet_wrap(~cohort) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
  labs(title = "Dose-response by histology (ever-smokers)", x = NULL, y = NULL, fill = "rho") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(fig_dir, "step2b_histology_stratified.pdf"), p_h, width = 10, height = 4)

## 7.3 Triangulation dot plot
p_tri <- tri %>% filter(exposure %in% c("py_gdc", "csi_lc")) %>%
  ggplot(aes(exposure, rho, color = method, shape = cell)) +
  geom_point(position = position_dodge(.6), size = 2) +
  geom_hline(yintercept = 0, linewidth = .3) +
  labs(title = "Multi-method triangulation (quantiseq / EPIC / MCPcounter)",
       x = NULL, y = "Spearman rho", color = "method", shape = "cell") +
  theme_bw()
ggsave(file.path(fig_dir, "step2b_multimethod_triangulation.pdf"), p_tri, width = 8, height = 4.5)

## 8. Summary of conclusions (auto-generated Markdown, for direct use in manuscript writing) ----
line <- function(...) paste0(...)
md <- c(
  "# 08_Step2b Diagnostic conclusions summary",
  "",
  sprintf("Generated: %s | Data: step2_immune_tme.csv (%d samples)", format(Sys.time(), "%Y-%m-%d %H:%M"), nrow(step2)),
  "",
  "## Q1 SBS4: proportion vs absolute (core)",
  "",
  "| Exposure | SBS4 proportion rho | SBS4 absolute rho | Sign flip |",
  "|---|---|---|---|",
  apply(sbs_cmp, 1, function(r) line("| ", r[["exposure"]], " | ",
      round(as.numeric(r[["rho_weight"]]), 3), " | ",
      round(as.numeric(r[["rho_abs"]]), 3), " | ",
      ifelse(isTRUE(as.logical(r[["direction_flip"]])), "**Yes**", "No"), " |")),
  "",
  "> Interpretation: if a sign flip occurs (proportion negative -> absolute positive), the negative correlation observed in Step2 **is a proportional-normalization artefact**,",
  "> not biological truth; the manuscript should adopt the absolute exposure sbs4_abs as the primary SBS4 metric.",
  "",
  "## Q2 Histology stratification",
  "",
  "Direction consistency:",
  apply(dir_consist, 1, function(r) line("- ", r[["exposure"]], " x ", r[["endpoint"]],
      ": LUAD rho=", round(as.numeric(r[["rho_LUAD"]]),3),
      ", LUSC rho=", round(as.numeric(r[["rho_LUSC"]]),3),
      " -> ", ifelse(isTRUE(as.logical(r[["same_direction"]])), "same direction", "**opposite direction**"))),
  "",
  "## Q3 Multi-method triangulation",
  "",
  apply(tri_sum, 1, function(r) line("- ", r[["cell"]], " x ", r[["exposure"]],
      ": three-method rho range [", r[["rho_min"]], ", ", r[["rho_max"]], "] ",
      ifelse(isTRUE(as.logical(r[["robust"]])), "**robust**", "not robust"))),
  "",
  "## Q4 Partial correlation (controlling for age/sex/stage/histology)",
  "",
  "Associations whose conclusion changed:",
  if (any(pp_tab$conclusion_changed)) {
    apply(pp_tab %>% filter(conclusion_changed), 1, function(r) line("- ", r[["exposure"]], " x ",
      r[["endpoint"]], ": raw p=", signif(as.numeric(r[["p_raw"]]),3),
      " -> partial p=", signif(as.numeric(r[["p_partial"]]),3)))
  } else "- None (all associations keep their conclusion after controlling for confounders)",
  "",
  "## Data-pitfall log",
  "",
  "- `age_y` is actually in days -> converted via `age = age_y / 365.25` (median 67.8 years)",
  "- `stage` mixes IA/IB/IIA/IIB/IIIA/IIIB/IV/I/II/III -> collapsed to `stage_grp` = I/II/III/IV",
  sprintf("- `q_T.cells.CD4` is 0 for %d/%d (clumped data) -> cross-validated with three deconvolution methods",
          sum(step2$q_T.cells.CD4 == 0, na.rm = TRUE), sum(is.finite(step2$q_T.cells.CD4))),
  "- `sbs4_abs` = SBS4 weight x TMB is an **approximate** absolute quantity (SBS4 is based on all SNVs, TMB on nonsynonymous only)",
  "")
## Note: filenames use ASCII. Under Windows + non-UTF-8 locale, R's file()/writeLines()
##   cannot reliably handle **Chinese filenames** (in practice reports `cannot open the connection` / Permission denied,
##   the same root cause as readLines failing on Chinese filenames). Chinese in the *content* is fine; only **filenames** must be ASCII.
writeLines(md, file.path(out_dir, "08_Step2b_diagnostic_summary.md"))

cat("\n================ Outputs ================\n")
cat("step2b_q1_sbs4_proportion_vs_absolute.csv\n",
    "step2b_q2_histology_stratified.csv\n",
    "step2b_q3_multimethod_triangulation.csv\n",
    "step2b_q4_partial_correlation.csv\n",
    "08_Step2b_diagnostic_summary.md\n",
    "figs/step2b_*.pdf\n", sep = "")
cat("\nDone. Next: read 08_Step2b_diagnostic_summary.md and use it to decide which exposure axis to use for Step3 modeling.\n")
############################################################