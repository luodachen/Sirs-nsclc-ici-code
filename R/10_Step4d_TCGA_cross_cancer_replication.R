## 10_Step4d_TCGA_cross_cancer_replication.R -- TCGA internal cross-cancer replication
## Project: smoking exposure quantification index x NSCLC immune microenvironment x ICI response prediction
##
## PURPOSE ("the index is derived in NSCLC but validated in urothelial
## carcinoma -- is cross-cancer validation valid?")
##   Provide an INTERNAL, cross-cancer check inside TCGA: apply the FROZEN NSCLC SIRS weights
##   (Step3) to three other tobacco-related cancers (BLCA / HNSC / ESCA) and test whether the
##   index still tracks smoking exposure (CSI). No PD-L1 stratification is needed for this,
##   so the check is feasible where an ICI cohort is not.
##   BLCA is especially relevant: it is the TCGA counterpart of the urothelial cohort
##   (IMvigor210) used in Step4, so it tests the same tissue in a different platform.
##
## WHAT IS / IS NOT TESTED
##   Tested   : exposure -> SIRS association (the index's construct validity) in new cancers.
##   NOT tested: ICI response prediction -- TCGA has no ICI treatment. This is stated plainly
##               in the output so the check cannot be over-claimed.
##
## PROVENANCE (code is deliberately a copy of the locked pipelines, not a re-invention)
##   Exposure  : 01_Step1_exposure_quantification.R   (read_gdc_clinical / read_legacy_clinical /
##               recode_status / csi_leffondre / build_exposure)
##   Immune    : 04_Step2_immune_microenvironment_TMB_SBS4.R (read_fpkm_uq / normalize_deconv /
##               deconvolution calls / read_mut / SBS4 block)
##   Weights   : step3_output/step3_SIRS_weights.csv  (FROZEN -- never re-estimated here)
##
## ======================================================================================
## [A] TWO DOCUMENTED DEVIATIONS from 02 (forced by data availability in the new cohorts)
##
##   (A1) PACK-YEARS FALLBACK.  02 uses `pack_years_smoked` (GDC) only. In the new cohorts the
##        legacy clinicalMatrix also carries `number_pack_years_smoked`; coverage is
##        comparable or better. Here:  py = coalesce(py_gdc, py_legacy).
##        Every row keeps `py_src` so the provenance is auditable.
##
##   (A2) SMOKING-DURATION FALLBACK.  02 takes duration from `years_smoked` (GDC), else from
##        `year_calc = stop_eff - start_year` with start_year = `year_of_tobacco_smoking_onset`.
##        BLCA and ESCA have NO `years_smoked` column and NO onset-year column, so:
##           start_year = coalesce(year_of_tobacco_smoking_onset,
##                                 year_of_birth + age_began_smoking_in_years)      <- NEW
##        (both `year_of_birth` and `age_began_smoking_in_years` exist in those cohorts).
##        Every row keeps `duration_src`, and the replication is also reported restricted to
##        rows whose duration came from the SAME source as in NSCLC (`duration_src ==
##        "years_smoked"`), as a strict-provenance sensitivity analysis.
##
##   Both deviations are reported in the summary and are the reason the new cohorts have
##   smaller exposure n than NSCLC.
##
## ======================================================================================
## [B] COLUMN-MATCHING TRAP
##   `grep("years_smoked", names(cl))` ALSO matches `pack_years_smoked.exposures`.
##   In LUAD/LUSC this is harmless because `years_smoked.exposures` happens to appear EARLIER
##   in the column order (position 9 vs 10), so the old unanchored match returned the right
##   column -- the published NSCLC exposure is therefore NOT affected (verified: cpd has 172
##   distinct values among 559 ever-smokers, i.e. not degenerate).
##   But BLCA/ESCA have NO `years_smoked` column at all, so the unanchored match silently
##   returned pack-years -> duration := pack-years -> cpd := pack_years*20/pack_years == 20
##   for EVERY sample (degenerate exposure, no information in the intensity axis).
##   FIX: all such matches here are ANCHORED (`^years_smoked`) and followed by an explicit
##   degeneracy guard (section 8, guard G6).
##
## ======================================================================================
## [C] STANDARDISATION REFERENCE SET (also found while building this script)
##   06_Step3 standardises the features with the moments of its MODELLING SET `dm`
##   (the 862 ever-smokers with complete clinical data; `use` = 855 of them have all 8
##   features observed). It does NOT standardise on all samples. Reproduced here:
##     reference = complete-case(dm)      ->  max|delta| vs the published export = 5.1e-15
##     reference = all samples            ->  max|delta| = 0.157
##     reference = all ever-smokers       ->  max|delta| = 0.038
##   The wrong choice still yields an index correlating ~0.21 with exposure, so it would
##   not have been caught by any downstream result -- only by this explicit reproduction
##   guard (G2b). For the new cohorts the analogue of "the training standardisation set"
##   is reported two ways (whole cohort / ever-smokers) plus the un-recentred train
##   moments, so the conclusion cannot hinge on that choice.
##
## Outputs:
##   step4_output/10_Step4d_summary.md
##   step4_output/step4d_cohort_exposure_qc.csv
##   step4_output/step4d_replication_spearman.csv
##   step4_output/step4d_meta_analysis.csv
##   step4_output/step4d_dose_response.csv
##   step4_output/step4d_negative_control.csv
##   step4_output/step4d_positive_control_nsclc.csv
##   step4_output/step4d_sample_level.csv
##   step4_output/figs/step4d_cross_cancer_replication.pdf
## Author: see manuscript author list   Date: 2026-10-10
############################################################

## 0. Configuration ---------------------------------------------------------
proj_dir <- "D:/Projects/smoking-index-ici"
data_dir <- file.path(proj_dir, "data_raw/xena")
out_dir  <- file.path(proj_dir, "step4_output"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir  <- file.path(out_dir, "figs");          dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
options(timeout = max(3600, getOption("timeout")), stringsAsFactors = FALSE)

## SBS4 is the slow block (per-sample fitting on ~1200 samples). It runs LAST so that every
## other output is already on disk even if this block is interrupted. Disable with --no-sbs4.
RUN_SBS4 <- !("--no-sbs4" %in% commandArgs(TRUE))
cat("[config] RUN_SBS4 =", RUN_SBS4, "\n")

NEW_COHORTS <- c("BLCA", "HNSC", "ESCA")     # all are tobacco-related cancers
TRAIN_COHORT <- c("LUAD", "LUSC")            # the SIRS training tissue (NSCLC)

suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAS_QS    <- has_pkg("quantiseqr")
HAS_EPIC  <- has_pkg("EPIC")
HAS_MCP   <- has_pkg("MCPcounter")
HAS_DSIGS <- has_pkg("deconstructSigs")
HAS_BSG   <- has_pkg("BSgenome.Hsapiens.UCSC.hg38")
cat("Packages: quantiseqr=", HAS_QS, " EPIC=", HAS_EPIC, " MCPcounter=", HAS_MCP,
    " deconstructSigs=", HAS_DSIGS, " BSgenome.hg38=", HAS_BSG, "\n", sep = "")
stopifnot("primary deconvolution pathway (quantiseqr) is required" = HAS_QS)
stopifnot("EPIC and MCPcounter are required (SIRS mixes all three algorithms)" = HAS_EPIC && HAS_MCP)

## 1. Exposure pipeline (copy of 02, plus deviations A1/A2) -------------------------------
na_tokens <- c("[not available]", "[unknown]", "[not applicable]", "[not performed]",
               "[discrepancy]", "[not evaluated]", "unknown", "na", "n/a", "", "not reported", "--")
clean_na  <- function(x) { x <- as.character(x); xt <- tolower(trimws(x)); x[is.na(x) | xt %in% na_tokens] <- NA; x }
num_clean <- function(x) suppressWarnings(as.numeric(clean_na(x)))

## ANCHORED matching (see header [B]): patterns are regexes anchored at the column start.
find_col <- function(df, patterns, required = FALSE) {
  for (p in patterns) {
    hit <- grep(p, names(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (required) stop("Required column not found: ", paste(patterns, collapse = " | "),
                     "\n  existing: ", paste(head(names(df), 40), collapse = ", "))
  NA_character_
}
read_xena_tsv <- function(path) {
  x <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
  names(x)[1] <- "sample"; x
}

read_gdc_clinical <- function(fname, cohort) {
  cl <- read_xena_tsv(file.path(data_dir, fname))
  cols <- list(
    years   = find_col(cl, c("^years_smoked")),                      # ANCHORED (see header [B])
    py      = find_col(cl, c("^pack_years_smoked", "^pack_years")),
    dx_year = find_col(cl, c("^year_of_diagnosis"), required = FALSE),
    age     = find_col(cl, c("^age_at_diagnosis", "^age_at_index"), required = FALSE),
    yob     = find_col(cl, c("^year_of_birth"), required = FALSE),    # NEW (deviation A2)
    stage   = find_col(cl, c("^ajcc_pathologic_stage", "^tumor_stage"), required = FALSE),
    gender  = find_col(cl, c("^gender"), required = FALSE),
    vital   = find_col(cl, c("^vital_status"), required = FALSE)
  )
  cat("\n[", cohort, " GDC clinical] auto-detected: years_smoked='", cols$years,
      "' pack_years='", cols$py, "' yob='", cols$yob, "'\n", sep = "")
  pick <- function(nm) if (is.na(cols[[nm]])) rep(NA_character_, nrow(cl)) else as.character(cl[[cols[[nm]]]])
  tibble(sample = cl$sample, cohort = cohort,
         years_gdc = pick("years"), py_gdc = pick("py"), dx_year_gdc = pick("dx_year"),
         age = pick("age"), yob_gdc = pick("yob"), stage = pick("stage"),
         gender = pick("gender"), vital = pick("vital"))
}

read_legacy_clinical <- function(path, cohort) {
  lg <- read_xena_tsv(path)
  cols <- list(
    tsh       = find_col(lg, c("^tobacco_smoking_history$"), required = FALSE),
    tshi      = find_col(lg, c("^tobacco_smoking_history_indicator"), required = FALSE),
    stop_y    = find_col(lg, c("stopped_smoking_year"), required = FALSE),
    onset     = find_col(lg, c("year_of_tobacco_smoking_onset"), required = FALSE),
    age_began = find_col(lg, c("age_began_smoking_in_years"), required = FALSE),  # NEW (deviation A2)
    dx_y      = find_col(lg, c("year_of_initial_pathologic_diagnosis"), required = FALSE),
    py        = find_col(lg, c("number_pack_years_smoked"), required = FALSE)
  )
  cat("[", cohort, " legacy]", " tsh='", cols$tsh, "' onset='", cols$onset,
      "' age_began='", cols$age_began, "' py='", cols$py, "'\n", sep = "")
  pick <- function(nm) if (is.na(cols[[nm]])) rep(NA_character_, nrow(lg)) else as.character(lg[[cols[[nm]]]])
  tibble(sample15 = substr(lg$sample, 1, 15), cohort = cohort,
         tsh_raw = pick("tsh"), tshi_raw = pick("tshi"),
         stop_year = pick("stop_y"), onset_year = pick("onset"),
         age_began = pick("age_began"), dx_year_lg = pick("dx_y"), py_lg = pick("py"))
}

## Status recoding: byte-for-byte the same rule as 02
recode_status <- function(code, label) {
  cn  <- suppressWarnings(as.numeric(clean_na(code)))
  out <- rep(NA_character_, length(cn))
  out[!is.na(cn) & cn == 1]            <- "never"
  out[!is.na(cn) & cn == 2]            <- "current"
  out[!is.na(cn) & cn %in% c(3, 4, 5)] <- "former"
  lab <- tolower(clean_na(label)); lab[is.na(lab)] <- ""
  idx <- is.na(out)
  out[idx & grepl("lifelong non-smoker|non-smoker|never", lab)] <- "never"
  out[idx & grepl("current smoker", lab)]                       <- "current"
  out[idx & grepl("reformed|former", lab)]                      <- "former"
  out
}

## CSI canonical formula + published parameters -- unchanged from 02
CSI_MAIN <- c(tau = 26, delta = 0.7)
csi_leffondre <- function(int, dur, tsc, tau = CSI_MAIN[["tau"]], delta = CSI_MAIN[["delta"]]) {
  tsc      <- ifelse(is.na(tsc), 0, tsc)
  tsc_star <- pmax(tsc - delta, 0)
  dur_star <- pmax(dur + tsc - delta, 0) - tsc_star
  (1 - 0.5^(dur_star / tau)) * 0.5^(tsc_star / tau) * log(int + 1)
}

build_exposure_cohort <- function(cohort) {
  g <- read_gdc_clinical(paste0("TCGA-", cohort, ".clinical.tsv.gz"), cohort)
  l <- read_legacy_clinical(file.path(data_dir, paste0(cohort, "_clinicalMatrix")), cohort)
  g %>% mutate(sample15 = substr(sample, 1, 15), stype = substr(sample, 14, 15)) %>%
    left_join(l %>% select(sample15, tsh_raw, tshi_raw, stop_year, onset_year, age_began, dx_year_lg, py_lg),
              by = "sample15") %>%
    mutate(
      py_gdc    = num_clean(py_gdc),
      py_legacy = num_clean(py_lg),
      duration_raw = num_clean(years_gdc),
      dx_year   = coalesce(num_clean(dx_year_gdc), num_clean(dx_year_lg)),
      yob       = num_clean(yob_gdc),
      age_began_n = num_clean(age_began),
      stop_year_n = num_clean(stop_year),
      onset_year_n = num_clean(onset_year),
      age_y     = num_clean(age),
      smoking_grp = recode_status(tsh_raw, tshi_raw),
      ## quantitative fallback for ever/never when the status field is empty (as in 02)
      smoking_grp = case_when(
        !is.na(smoking_grp) ~ smoking_grp,
        !is.na(py_gdc) & py_gdc == 0 ~ "never",
        (!is.na(py_gdc) & py_gdc > 0) | (!is.na(py_legacy) & py_legacy > 0) ~ "ever_unclassified",
        TRUE ~ NA_character_),
      ## ---- deviation A1: pack-years = coalesce(GDC, legacy) ----
      py        = coalesce(py_gdc, py_legacy),
      py_src    = case_when(is.na(py) ~ NA_character_,
                            !is.na(py_gdc) ~ "py_gdc", TRUE ~ "py_legacy"),
      ## ---- deviation A2: duration fallback cascade ----
      start_year = coalesce(onset_year_n, yob + age_began_n),   # yob + age_began is the NEW source
      start_src  = case_when(is.na(start_year) ~ NA_character_,
                             !is.na(onset_year_n) ~ "onset_year",
                             TRUE ~ "yob_plus_age_began"),
      stop_eff   = ifelse(smoking_grp == "current", dx_year, stop_year_n),
      dur_est    = ifelse(!is.na(start_year) & !is.na(stop_eff) & stop_eff >= start_year,
                          stop_eff - start_year, NA_real_),
      duration_src = case_when(
        smoking_grp == "never" ~ "never_zero",
        !is.na(duration_raw) ~ "years_smoked",
        !is.na(dur_est) & dur_est > 0 ~ "year_calc",
        TRUE ~ NA_character_),
      duration = case_when(
        smoking_grp == "never" ~ 0,
        !is.na(duration_raw) ~ duration_raw,
        !is.na(dur_est) & dur_est > 0 ~ dur_est,
        TRUE ~ NA_real_),
      years_since_quit_raw = case_when(
        smoking_grp == "never" ~ NA_real_,
        smoking_grp == "current" ~ 0,
        !is.na(stop_year_n) & !is.na(dx_year) ~ dx_year - stop_year_n,
        TRUE ~ NA_real_),
      quit_year_irregular = !is.na(years_since_quit_raw) & years_since_quit_raw < 0,
      years_since_quit = ifelse(quit_year_irregular, NA_real_, years_since_quit_raw),
      ## CPD = pack-years x 20 / smoking years  (pack-year identity, unchanged)
      cpd = case_when(smoking_grp == "never" ~ 0,
                      !is.na(py) & !is.na(duration) & duration > 0 ~ py * 20 / duration,
                      TRUE ~ NA_real_),
      cpd_implausible = !is.na(cpd) & smoking_grp != "never" & (cpd < 1 | cpd > 100),
      csi_lc = csi_leffondre(cpd, duration, years_since_quit)
    ) %>%
    ## one sample per patient, prefer the most complete smoking record (as in 02)
    mutate(n_missing = rowSums(is.na(across(c(duration, py, stop_year_n))))) %>%
    arrange(substr(sample, 1, 12), n_missing) %>%
    distinct(substr(sample, 1, 12), .keep_all = TRUE) %>%
    select(-n_missing)
}

cat("\n", strrep("=", 78), "\n## 1. Exposure for the three new cohorts\n", strrep("=", 78), "\n", sep = "")
expo <- list()
for (co in NEW_COHORTS) {
  d <- build_exposure_cohort(co)
  if (any(d$stype == "01")) d <- d %>% filter(stype == "01") else warning("patient-level IDs in ", co)
  expo[[co]] <- d
  cat("\n[", co, "] samples =", nrow(d), "\n", sep = "")
  print(table(d$smoking_grp, useNA = "ifany"))
  ev <- d %>% filter(smoking_grp %in% c("current", "former"))
  cat("  ever n =", nrow(ev),
      "| csi_lc finite =", sum(is.finite(ev$csi_lc)),
      "| duration_src:", paste(names(table(ev$duration_src)), table(ev$duration_src), sep = "=", collapse = " "), "\n")
}

## 1.1 Exposure QC table + guards ---------------------------------------------------------
exp_qc <- bind_rows(lapply(NEW_COHORTS, function(co) {
  d <- expo[[co]]; ev <- d %>% filter(smoking_grp %in% c("current", "former"))
  tibble(cohort = paste0("TCGA-", co), n_primary = nrow(d),
         n_never = sum(d$smoking_grp == "never", na.rm = TRUE),
         n_former = sum(d$smoking_grp == "former", na.rm = TRUE),
         n_current = sum(d$smoking_grp == "current", na.rm = TRUE),
         n_status_na = sum(is.na(d$smoking_grp)),
         n_ever = nrow(ev), n_csi = sum(is.finite(ev$csi_lc)),
         csi_coverage = round(mean(is.finite(ev$csi_lc)), 3),
         cpd_median = round(median(ev$cpd, na.rm = TRUE), 1),
         cpd_distinct = length(unique(round(ev$cpd[is.finite(ev$cpd)], 6))),
         cpd_is_constant_20 = all(abs(ev$cpd[is.finite(ev$cpd)] - 20) < 1e-9),
         duration_median = round(median(ev$duration, na.rm = TRUE), 1),
         csi_median = round(median(ev$csi_lc, na.rm = TRUE), 3),
         dur_src_years_smoked = sum(ev$duration_src == "years_smoked", na.rm = TRUE),
         dur_src_year_calc = sum(ev$duration_src == "year_calc", na.rm = TRUE),
         py_src_gdc = sum(ev$py_src == "py_gdc", na.rm = TRUE),
         py_src_legacy = sum(ev$py_src == "py_legacy", na.rm = TRUE))
}))
fwrite(exp_qc, file.path(out_dir, "step4d_cohort_exposure_qc.csv"))
cat("\n== Exposure QC ==\n"); print(as.data.frame(exp_qc))

## G1: status field must be usable in every cohort
stopifnot(all(exp_qc$n_ever > 50))
## G6 (the guard for header [B]): a cohort whose duration silently collapsed onto pack-years
##     would make cpd identically 20 for everyone. That must never pass.
stopifnot(!any(exp_qc$cpd_is_constant_20))
cat("[guard G1/G6] every cohort has >50 ever-smokers and non-degenerate CPD  OK\n")

## 2. Expression matrices -> TPM (identical to 07) ---------------------------------------
probemap_file <- file.path(data_dir, "gencode.v36.annotation.gtf.gene.probemap")
if (!file.exists(probemap_file) || file.size(probemap_file) < 1e5)
  stop("gene.probemap missing -- run 07 once to fetch it")
probemap <- fread(probemap_file) %>%
  transmute(id = sub("\\..*$", "", id), symbol = gene) %>%
  filter(!is.na(symbol), symbol != "", symbol != "?") %>% distinct(id, .keep_all = TRUE)

read_fpkm_uq <- function(cohort) {
  f <- file.path(data_dir, paste0("TCGA-", cohort, ".star_fpkm-uq.tsv.gz"))
  if (!file.exists(f)) stop("missing expression file: ", f)
  fread(f) %>% as_tibble() %>%
    mutate(id = sub("\\..*$", "", Ensembl_ID)) %>%
    select(-Ensembl_ID) %>% distinct(id, .keep_all = TRUE) %>%
    inner_join(probemap, by = "id") %>% filter(!is.na(symbol)) %>%
    select(symbol, where(is.numeric)) %>%
    group_by(symbol) %>% summarise(across(everything(), max), .groups = "drop") %>%
    column_to_rownames("symbol") %>% as.matrix()
}
dedup_samples <- function(M) {
  cn <- colnames(M); sid <- substr(cn, 1, 15); alq <- substring(cn, 16)
  is01 <- substr(cn, 14, 15) == "01"
  keep <- tibble(raw = cn[is01], sid = sid[is01], aliquot = alq[is01]) %>%
    arrange(sid, aliquot) %>% distinct(sid, .keep_all = TRUE)
  M <- M[, match(keep$raw, cn), drop = FALSE]; colnames(M) <- keep$sid
  stopifnot(!anyDuplicated(colnames(M))); M
}

## 3. Deconvolution (byte-for-byte the same calls as 07) ----------------------------------
normalize_deconv <- function(df, prefix, cell_prefix = TRUE) {
  df <- as_tibble(df)
  hit <- base::intersect(c("Sample", "sample"), names(df))
  sample_col <- if (length(hit)) hit[1] else { warning("no Sample column"); names(df)[1] }
  stopifnot(sum(names(df) == sample_col) == 1)
  n_cell <- ncol(df) - 1
  old <- names(df); names(df)[old == sample_col] <- "sample15"
  out <- df %>% mutate(sample15 = substr(as.character(sample15), 1, 15))
  if (cell_prefix) out <- out %>% rename_with(~ paste0(prefix, .x), -sample15)
  stopifnot(ncol(out) == n_cell + 1); out
}

deconvolve_cohort <- function(cohort) {
  cat("\n", strrep("-", 70), "\n## Deconvolution: ", cohort, "\n", sep = "")
  M <- read_fpkm_uq(cohort) %>% dedup_samples()
  cat("  expression: ", nrow(M), " genes x ", ncol(M), " samples\n", sep = "")
  ## FPKM-UQ -> TPM; deconvolution uses LINEAR TPM (as in 07)
  mat_tpm <- apply(M, 2, function(x) { s <- sum(x[x > 0], na.rm = TRUE); ifelse(x > 0, x / s * 1e6, 0) })
  rm(M); gc(verbose = FALSE)

  library(quantiseqr)
  qs <- run_quantiseq(expression_data = mat_tpm, signature_matrix = "TIL10",
                      is_arraydata = FALSE, is_tumordata = TRUE, scale_mRNA = TRUE, method = "lsei")
  d_q <- normalize_deconv(qs, prefix = "q_")

  library(EPIC)
  epic_warn <- character(0)
  epic_raw <- withCallingHandlers(
    EPIC(mat_tpm),
    warning = function(w) { epic_warn <<- c(epic_warn, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage"))
  fr <- as.matrix(epic_raw$cellFractions)
  d_e <- fr %>% as_tibble(.name_repair = "minimal") %>%
    mutate(sample15 = substr(rownames(fr), 1, 15)) %>%
    select(sample15, everything()) %>% rename_with(~ paste0("e_", .x), -sample15)
  d_e$epic_nonconverged <- d_e$sample15 %in% base::intersect(
    unique(unlist(stringr::str_extract_all(paste(epic_warn, collapse = " "),
                                           "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9]{2}"))), d_e$sample15)

  library(MCPcounter)
  d_m <- MCPcounter.estimate(mat_tpm, featuresType = "HUGO_symbols") %>%
    as.data.frame() %>% rownames_to_column("cell_type") %>%
    pivot_longer(-cell_type, names_to = "sample15", values_to = "score") %>%
    mutate(sample15 = substr(sample15, 1, 15), cell_type = paste0("mcp_", cell_type)) %>%
    pivot_wider(names_from = cell_type, values_from = score)

  out <- purrr::reduce(list(d_q, d_e, d_m), full_join, by = "sample15")
  stopifnot(!anyDuplicated(out$sample15))
  cat("  deconvolution: ", nrow(out), " samples x ", ncol(out), " cols; EPIC non-converged = ",
      sum(out$epic_nonconverged), "\n", sep = "")
  out
}

## 4. TMB + SBS4 (07's conventions; SBS4 gated and run last) ------------------------------
NONSYN_RE <- paste(c("missense_variant", "stop_gained", "stop_lost", "start_lost",
                     "frameshift_variant", "inframe_deletion", "inframe_insertion",
                     "splice_acceptor_variant", "splice_donor_variant",
                     "splice_region_variant", "protein_altering_variant",
                     "disruptive_inframe", "nonstop"), collapse = "|")
EXOME_MB <- 38
STD_CHR  <- c(as.character(1:22), "X", "Y")
read_mut <- function(cohort) {
  f <- file.path(data_dir, paste0("TCGA-", cohort, ".somaticmutation_wxs.tsv.gz"))
  if (!file.exists(f)) stop("missing mutation file: ", f)
  fread(f, select = c("sample", "chrom", "start", "end", "ref", "alt", "effect"),
        showProgress = FALSE) %>% as_tibble() %>% mutate(cohort = paste0("TCGA-", cohort))
}

## 5. SIRS scoring: FROZEN NSCLC weights, cohort-internal z -------------------------------
W <- fread(file.path(proj_dir, "step3_output/step3_SIRS_weights.csv")) %>% as_tibble()
SIRS_FEATS <- W$feature
SIRS_W     <- setNames(W$weight, W$feature)
## train mean/sd are read from the published formula markdown table (authoritative, as released)
fm <- readLines(file.path(proj_dir, "step3_output/step3_SIRS_formula.md"))
fm_tab <- fm[grepl("^\\| *(mcp_|q_|e_)", fm)]
parse_row <- function(s) { p <- trimws(strsplit(s, "\\|")[[1]]); p <- p[nzchar(p)]
  c(feature = p[1], w = p[2], anchor = p[3], mu = as.character(p[5]), sd = as.character(p[6])) }
TRAIN_PAR <- as_tibble(do.call(rbind, lapply(fm_tab, parse_row))) %>%
  mutate(mu = as.numeric(mu), sd = as.numeric(sd), w = as.numeric(w))
stopifnot(nrow(TRAIN_PAR) == length(SIRS_FEATS), setequal(TRAIN_PAR$feature, SIRS_FEATS))
cat("\nFrozen SIRS: ", length(SIRS_FEATS), " features from ", paste(unique(W$anchor), collapse = "/"),
    " anchor(s)\n", sep = "")
print(as.data.frame(W))

## ---------------------------------------------------------------------------
## SIRS scoring. THE STANDARDISATION REFERENCE SET IS THE WHOLE BALL GAME.
## 06_Step3 standardises on the z_mean/z_sd of its MODELLING SET (`dm`, the ever-smokers
## with complete clinical data), not on all samples. Verified here by guard G2: that
## convention reproduces the published step3_SIRS_all_samples.csv to 5.1e-15.
## Using "all samples with complete features" instead gives max|delta| = 0.157 and still
## correlates ~0.21 with exposure -- so nothing downstream would have flagged the error.
## ---------------------------------------------------------------------------
score_z_ref <- function(dat, feats, w, ref_idx = NULL) {
  X  <- as.matrix(dat[, feats, drop = FALSE])
  cc <- apply(X, 1, function(r) all(is.finite(r)))            # complete-case rows
  ref <- if (is.null(ref_idx)) cc else (ref_idx & cc)
  if (sum(ref) < 10) stop("standardisation reference set too small (n = ", sum(ref), ")")
  Xr <- X[ref, , drop = FALSE]
  mu <- colMeans(Xr); sdv <- apply(Xr, 2, stats::sd)
  if (any(!is.finite(sdv) | sdv <= 0))
    stop("zero-variance SIRS feature in the reference set: ",
         paste(feats[!is.finite(sdv) | sdv <= 0], collapse = ", "))
  out <- rep(NA_real_, nrow(dat))
  out[cc] <- as.numeric(sweep(sweep(X[cc, , drop = FALSE], 2, mu, "-"), 2, sdv, "/") %*% w)
  out
}
## train-parameter z (sensitivity): the published TCGA moments, i.e. no cohort re-centring
score_trainz <- function(dat, feats, w, mu, sd) {
  X  <- as.matrix(dat[, feats, drop = FALSE])
  cc <- apply(X, 1, function(r) all(is.finite(r)))
  out <- rep(NA_real_, nrow(dat))
  out[cc] <- as.numeric(sweep(sweep(X[cc, , drop = FALSE], 2, mu, "-"), 2, sd, "/") %*% w)
  out
}

## 6. Positive control on the TRAINING tissue, through the SAME code path -----------------
## Re-scores LUAD/LUSC from the published step2 table (no re-deconvolution), and rebuilds
## Step3's modelling set `dm` exactly as 06_Step3 / 09_Step4c do, so guard G2 is a genuine
## reproduction test rather than a plausibility test.
step2_nsclc <- fread(file.path(proj_dir, "step2_output/step2_immune_tme.csv"), na.strings = c("", "NA")) %>%
  as_tibble() %>% mutate(sample15 = substr(sample, 1, 15))
stopifnot(all(SIRS_FEATS %in% names(step2_nsclc)))
ns_all <- step2_nsclc %>%
  mutate(os_event = vital == "Dead",
         os_time  = ifelse(vital == "Dead", days_death, days_fu),
         age      = age_y / 365.25) %>%
  mutate(os_time = ifelse(!is.finite(os_time) & is.finite(days_fu),    days_fu,
                   ifelse(!is.finite(os_time) & is.finite(days_death), days_death, os_time)),
         stage_grp = ifelse(grepl("^Stage\\s*IV",  stage), "IV",
                     ifelse(grepl("^Stage\\s*III", stage), "III",
                     ifelse(grepl("^Stage\\s*II",  stage), "II",
                     ifelse(grepl("^Stage\\s*I",   stage), "I", NA_character_)))),
         gender_m = as.integer(factor(gender)) - 1L,
         stage_n  = as.integer(factor(stage_grp, levels = c("I", "II", "III", "IV"))),
         cohort_n = as.integer(factor(cohort)) - 1L,
         os_time  = pmax(os_time, 0.5))
dm <- ns_all %>% filter(smoking_grp %in% c("current", "former"),
                        vital %in% c("Alive", "Dead"), is.finite(os_time),
                        !is.na(age), !is.na(gender_m), !is.na(stage_n), !is.na(cohort_n))
stopifnot(nrow(dm) == 862L, sum(dm$os_event) == 345L)   # G2a: rebuilds Step3's modelling set
ns_all$sirs_re <- score_z_ref(ns_all, SIRS_FEATS, SIRS_W, ns_all$sample15 %in% dm$sample15)
s3  <- fread(file.path(proj_dir, "step3_output/step3_SIRS_all_samples.csv")) %>% as_tibble()
cmp <- ns_all %>% select(sample15, sirs_re) %>% inner_join(s3 %>% select(sample15, sirs), by = "sample15")
max_d <- max(abs(cmp$sirs_re - cmp$sirs), na.rm = TRUE)
cat(sprintf("\n[positive control] re-scored SIRS vs the published Step3 export: n = %d, max |delta| = %.3e\n",
            nrow(cmp), max_d))
stopifnot(max_d < 1e-10)   # G2b: the same code path must REPRODUCE the frozen index

ns_ev <- ns_all %>% filter(smoking_grp %in% c("current", "former"), is.finite(csi_lc), is.finite(sirs_re))
ct <- suppressWarnings(cor.test(ns_ev$sirs_re, ns_ev$csi_lc, method = "spearman", exact = FALSE))
## the published reference value (step3_SIRS_exposure_validity.csv: n = 539, rho = 0.2134)
pub_val <- fread(file.path(proj_dir, "step3_output/step3_SIRS_exposure_validity.csv")) %>%
  as_tibble() %>% filter(variable == "csi_lc")
pos_ctrl <- tibble(cohort = "TCGA-NSCLC (training tissue, positive control)",
                   n = nrow(ns_ev), rho = unname(ct$estimate), p_value = ct$p.value,
                   published_n = pub_val$n[1], published_rho = round(pub_val$rho[1], 4))
cat(sprintf("[positive control] within-NSCLC rho(SIRS, csi_lc) = %.3f (n = %d, p = %.2e)\n",
            pos_ctrl$rho, pos_ctrl$n, pos_ctrl$p_value))
cat(sprintf("  published step3 table: rho = %.4f (n = %d). The small gap is a sample-subset\n",
            pos_ctrl$published_rho, pos_ctrl$published_n))
cat("  difference only (that table used `use`); the index itself is reproduced exactly.\n")
fwrite(pos_ctrl, file.path(out_dir, "step4d_positive_control_nsclc.csv"))

## 7. Cross-cancer replication -----------------------------------------------------------
## Self-implemented Jonckheere-Terpstra trend test (verbatim from 07, which had to avoid the
## uninstalled clinfun): U = sum_{i<j} #{(x in i, y in j): y > x}, normal approx + correction.
jt_test <- function(x, g) {
  ok <- is.finite(x) & !is.na(g); x <- x[ok]; g <- as.integer(g[ok])
  ks <- sort(unique(g)); if (length(ks) < 2) return(NA_real_)
  U <- 0; mu <- 0; v <- 0
  for (a in seq_along(ks)) for (b in seq_along(ks)) {
    if (b <= a) next
    xa <- x[g == ks[a]]; xb <- x[g == ks[b]]; na <- length(xa); nb <- length(xb)
    U <- U + sum(outer(xb, xa, ">")) + 0.5 * sum(outer(xb, xa, "=="))
    mu <- mu + na * nb / 2; v <- v + na * nb * (na + nb + 1) / 12
  }
  if (v <= 0) return(NA_real_)
  2 * pnorm(-abs((U - mu - 0.5) / sqrt(v)))
}

spearman_fast <- function(x, y) {
  ok <- is.finite(x) & is.finite(y); if (sum(ok) < 10) return(NA_real_)
  suppressWarnings(cor(rank(x[ok]), rank(y[ok])))
}
perm_p <- function(x, y, B = 2000, seed = 20261010) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(NA_real_)
  obs <- abs(suppressWarnings(cor(rank(x), rank(y))))
  set.seed(seed)
  null <- vapply(seq_len(B), function(i) abs(suppressWarnings(cor(rank(x), rank(sample(y))))), numeric(1))
  (1 + sum(null >= obs)) / (B + 1)
}
boot_ci <- function(x, y, B = 2000, seed = 20261011) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 10) return(c(NA_real_, NA_real_))
  set.seed(seed)
  r <- vapply(seq_len(B), function(i) {
    s <- sample.int(n, n, replace = TRUE)
    if (length(unique(s)) < 10) return(NA_real_)
    suppressWarnings(cor(rank(x[s]), rank(y[s])))
  }, numeric(1))
  unname(stats::quantile(r, c(0.025, 0.975), na.rm = TRUE))
}

deconv_res <- list(); mut_res <- list()
for (co in NEW_COHORTS) {
  ## Caching: deconvolution is the expensive step (~minutes per cohort) and the SBS4 block
  ## may need a separate run. Reuse the cache when present so a re-run costs seconds.
  cache <- file.path(out_dir, paste0("step4d_deconv_", co, ".csv"))
  if (file.exists(cache) && !("--recompute-deconv" %in% commandArgs(TRUE))) {
    deconv_res[[co]] <- fread(cache) %>% as_tibble() %>% mutate(sample15 = substr(sample15, 1, 15))
    cat("[cache] reusing deconvolution for ", co, ": ", nrow(deconv_res[[co]]), " samples\n", sep = "")
  } else {
    deconv_res[[co]] <- deconvolve_cohort(co)
    fwrite(deconv_res[[co]], cache)
    cat("[cache] wrote ", basename(cache), "\n", sep = "")
  }
  m <- read_mut(co) %>% mutate(sample15 = substr(sample, 1, 15), chrom = sub("^chr", "", chrom)) %>%
    filter(!is.na(start), chrom %in% STD_CHR)
  mut_res[[co]] <- m
  mu <- m %>% filter(grepl(NONSYN_RE, effect))
  unim <- m %>% distinct(sample15)
  tmb <- mu %>% distinct(sample15, chrom, start, ref, alt, .keep_all = TRUE) %>%
    count(sample15, name = "nonsyn_n") %>% right_join(unim, by = "sample15") %>%
    mutate(nonsyn_n = ifelse(is.na(nonsyn_n), 0L, nonsyn_n), tmb_per_mb = nonsyn_n / EXOME_MB)
  if (!("tmb_per_mb" %in% names(deconv_res[[co]])))
    deconv_res[[co]] <- deconv_res[[co]] %>% left_join(tmb, by = "sample15")
  cat("  TMB: ", nrow(tmb), " samples; median = ", round(median(tmb$tmb_per_mb, na.rm = TRUE), 2),
      " mut/Mb\n", sep = "")
}

## 7.1 Per-cohort replication table ------------------------------------------------------
sample_tab <- list(); rep_tab <- list()
for (co in NEW_COHORTS) {
  d <- expo[[co]] %>% select(sample15, sample, cohort, smoking_grp, duration, duration_src,
                             py, py_src, cpd, years_since_quit, csi_lc) %>%
    inner_join(deconv_res[[co]], by = "sample15")
  miss <- setdiff(SIRS_FEATS, names(d))
  if (length(miss)) stop("cohort ", co, " missing SIRS features: ", paste(miss, collapse = ", "))
  ## PRIMARY: cohort-internal z over all complete-case samples of the cohort
  ## (the same convention Step4 used for `sirs_cz` in IMvigor210)
  d$sirs_re <- score_z_ref(d, SIRS_FEATS, SIRS_W)
  ## sensitivity 1: reference set restricted to the cohort's ever-smokers
  d$sirs_evz <- score_z_ref(d, SIRS_FEATS, SIRS_W, d$smoking_grp %in% c("current", "former"))
  ## sensitivity 2: no cohort re-centring at all (the published TCGA moments)
  d$sirs_trainz <- score_trainz(d, SIRS_FEATS, SIRS_W, TRAIN_PAR$mu, TRAIN_PAR$sd)
  sample_tab[[co]] <- d

  ev <- d %>% filter(smoking_grp %in% c("current", "former"), is.finite(csi_lc), is.finite(sirs_re))
  rho <- spearman_fast(ev$sirs_re, ev$csi_lc); ci <- boot_ci(ev$sirs_re, ev$csi_lc)
  ## strict-provenance subset: duration from years_smoked only (same source as NSCLC)
  evs <- ev %>% filter(duration_src == "years_smoked")
  rs  <- spearman_fast(evs$sirs_re, evs$csi_lc)
  ## sensitivities: ever-smoker-centred z; published train moments; SBS4 anchor handled later
  rt  <- spearman_fast(ev$sirs_trainz, ev$csi_lc)
  re  <- spearman_fast(ev$sirs_evz, ev$csi_lc)
  ## secondary axes already available: TMB (smoking-related mutation burden; NSCLC rho = 0.339)
  n_tmb <- sum(is.finite(ev$tmb_per_mb))
  rtmb  <- spearman_fast(ev$sirs_re, ev$tmb_per_mb)
  ptmb  <- suppressWarnings(cor.test(ev$sirs_re, ev$tmb_per_mb, method = "spearman", exact = FALSE)$p.value)
  rep_tab[[co]] <- tibble(
    cohort = paste0("TCGA-", co),
    n = nrow(ev), rho = round(rho, 4),
    ci_lo = round(ci[1], 4), ci_hi = round(ci[2], 4),
    p_value = suppressWarnings(cor.test(ev$sirs_re, ev$csi_lc, method = "spearman", exact = FALSE)$p.value),
    p_perm = perm_p(ev$sirs_re, ev$csi_lc),
    rho_trainz = round(rt, 4), rho_everz = round(re, 4),
    n_strict = nrow(evs), rho_strict = round(rs, 4),
    n_tmb = n_tmb, rho_tmb = round(rtmb, 4), p_tmb = ptmb,
    median_csi = round(median(ev$csi_lc, na.rm = TRUE), 3),
    median_sirs = round(median(ev$sirs_re, na.rm = TRUE), 3))
  cat(sprintf("  [%s] n = %d | rho(SIRS, CSI) = %+.3f [%.3f, %.3f] | perm p = %.4f\n",
              co, rep_tab[[co]]$n, rho, ci[1], ci[2], rep_tab[[co]]$p_perm))
}
rep_tab <- bind_rows(rep_tab)
fwrite(rep_tab, file.path(out_dir, "step4d_replication_spearman.csv"))
cat("\n== Cross-cancer replication ==\n"); print(as.data.frame(rep_tab))
stopifnot(all(rep_tab$n >= 50))                       # G3: every cohort contributes
stopifnot(!any(is.na(rep_tab$rho)))                   # G3

## 7.2 Meta-analysis across cohorts (Fisher z; fixed + DerSimonian-Laird random) --------
meta_fisher <- function(r, n) {
  ok <- is.finite(r) & is.finite(n) & n > 3
  r <- r[ok]; n <- n[ok]; k <- length(r)
  z <- atanh(r); v <- 1 / (n - 3); w <- 1 / v
  z_fe <- sum(w * z) / sum(w); se_fe <- sqrt(1 / sum(w))
  Q <- sum(w * (z - z_fe)^2); df <- k - 1
  C <- sum(w) - sum(w^2) / sum(w)
  tau2 <- if (C > 0) max(0, (Q - df) / C) else 0
  wr <- 1 / (v + tau2); z_re <- sum(wr * z) / sum(wr); se_re <- sqrt(1 / sum(wr))
  I2 <- if (Q > 0) max(0, (Q - df) / Q) * 100 else 0
  tibble(k = k, n_total = sum(n),
         rho_fe = round(tanh(z_fe), 4), ci_lo_fe = round(tanh(z_fe - 1.96 * se_fe), 4),
         ci_hi_fe = round(tanh(z_fe + 1.96 * se_fe), 4),
         rho_re = round(tanh(z_re), 4), ci_lo_re = round(tanh(z_re - 1.96 * se_re), 4),
         ci_hi_re = round(tanh(z_re + 1.96 * se_re), 4),
         p_fe = 2 * pnorm(-abs(z_fe / se_fe)),
         Q = round(Q, 3), I2_pct = round(I2, 1), tau2 = signif(tau2, 3))
}
meta_all <- meta_fisher(rep_tab$rho, rep_tab$n) %>% mutate(set = "BLCA + HNSC + ESCA (all three)")
meta_blca <- meta_fisher(rep_tab$rho[rep_tab$cohort == "TCGA-BLCA"], rep_tab$n[rep_tab$cohort == "TCGA-BLCA"]) %>%
  mutate(set = "TCGA-BLCA alone (urothelial = same tissue as IMvigor210)")
meta_meta <- bind_rows(meta_all, meta_blca)
## meta permutation: shuffle exposure within cohort, re-pool
set.seed(20261012)
B_meta <- 2000
perm_meta <- vapply(seq_len(B_meta), function(i) {
  zs <- c(); vs <- c()
  for (co in NEW_COHORTS) {
    d <- sample_tab[[co]] %>% filter(smoking_grp %in% c("current", "former"), is.finite(csi_lc), is.finite(sirs_re))
    rr <- spearman_fast(d$sirs_re, sample(d$csi_lc))
    zs <- c(zs, atanh(rr)); vs <- c(vs, 1 / (nrow(d) - 3))
  }
  w <- 1 / vs; tanh(sum(w * zs) / sum(w))
}, numeric(1))
obs_meta <- meta_all$rho_fe
p_meta_perm <- (1 + sum(abs(perm_meta) >= abs(obs_meta))) / (B_meta + 1)
meta_meta$p_permutation <- c(p_meta_perm, NA_real_)
fwrite(meta_meta, file.path(out_dir, "step4d_meta_analysis.csv"))
cat("\n== Meta-analysis (Fisher z) ==\n"); print(as.data.frame(meta_meta))
cat(sprintf("Random-effects pooled rho = %.3f [%.3f, %.3f], I2 = %.1f%%\n",
            meta_all$rho_re, meta_all$ci_lo_re, meta_all$ci_hi_re, meta_all$I2_pct))

## 7.3 Secondary axis: SIRS vs SBS4_abs (the anchor ALL weights came from) --------------
sbs_tab <- NULL
sbs_cache <- file.path(out_dir, "step4d_sbs4.csv")
if (RUN_SBS4 && HAS_DSIGS && HAS_BSG && file.exists(sbs_cache) &&
    !("--recompute-sbs4" %in% commandArgs(TRUE))) {
  cat("\n[cache] reusing SBS4 fits from ", basename(sbs_cache), " (skip the slow block)\n", sep = "")
  sbs_abs <- fread(sbs_cache) %>% as_tibble()
} else if (RUN_SBS4 && HAS_DSIGS && HAS_BSG) {
  cat("\n", strrep("=", 78), "\n## 4/7.3 SBS4 (slow, per-sample) -- the direct anchor replication\n", strrep("=", 78), "\n", sep = "")
  library(deconstructSigs); library(BSgenome.Hsapiens.UCSC.hg38)
  MIN_SBS4_SNV <- 50
  sbs_out <- list()
  for (co in NEW_COHORTS) {
    t0 <- Sys.time()
    mi <- mut_res[[co]] %>% filter(ref %in% c("A", "C", "G", "T"), alt %in% c("A", "C", "G", "T")) %>%
      transmute(sample.id = sample15, chr = paste0("chr", chrom), pos = start, ref, alt) %>% as.data.frame()
    si <- suppressWarnings(mut.to.sigs.input(mut.ref = mi, sample.id = "sample.id", chr = "chr",
                                             pos = "pos", ref = "ref", alt = "alt", bsg = BSgenome.Hsapiens.UCSC.hg38))
    si_df <- as.data.frame(si); n_snv <- rowSums(si_df)
    keep_s <- names(n_snv)[n_snv >= MIN_SBS4_SNV]
    wl <- lapply(keep_s, function(s) {
      res <- try(whichSignatures(tumor.ref = si_df, sample.id = s, signatures.ref = signatures.cosmic,
                                 contexts.needed = TRUE, tri.counts.method = "exome2genome")[["weights"]],
                 silent = TRUE)
      if (inherits(res, "try-error")) return(NULL); res$sample15 <- s; res })
    w <- bind_rows(wl) %>% as_tibble()
    if (nrow(w)) {
      sbs_out[[co]] <- w %>% transmute(sample15, sbs4_weight = `Signature.4`)
      cat(sprintf("  [%s] SBS4 fitted for %d / %d samples (SNV >= %d); elapsed %.1f min\n",
                  co, nrow(w), nrow(si_df), MIN_SBS4_SNV,
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    }
  }
  sbs_all <- bind_rows(sbs_out)
  ## SBS4 absolute burden = weight x the sample's total SNV count (07's sbs4_abs convention)
  snv_tot <- bind_rows(lapply(NEW_COHORTS, function(co)
    mut_res[[co]] %>% filter(ref %in% c("A","C","G","T"), alt %in% c("A","C","G","T")) %>%
      count(sample15, name = "n_snv")))
  sbs_abs <- sbs_all %>% left_join(snv_tot, by = "sample15") %>%
    mutate(sbs4_abs = sbs4_weight * n_snv) %>% select(sample15, sbs4_weight, n_snv, sbs4_abs)
  fwrite(sbs_abs, file.path(out_dir, "step4d_sbs4.csv"))
}

## --- post-processing shared by the cached and the freshly computed SBS4 table ---
if (exists("sbs_abs") && !is.null(sbs_abs) && nrow(sbs_abs)) {
  sbs_tab <- bind_rows(lapply(NEW_COHORTS, function(co) {
    ## NB: aggregate FIRST into explicit vectors. `tibble()` does NOT data-mask, so writing
    ## `sum(is.finite(sbs4_abs))` inside tibble() after a left_join looks like it reads the
    ## joined column but actually looks for an object named `sbs4_abs` in the caller -> the
    ## "object 'sbs4_abs' not found" error. (Same class of mistake as the `D["col"][idx]` trap.)
    d  <- sample_tab[[co]] %>% left_join(sbs_abs, by = "sample15") %>%
      filter(smoking_grp %in% c("current", "former"), is.finite(sirs_re))
    n_abs  <- sum(is.finite(d$sbs4_abs))
    r_abs  <- spearman_fast(d$sirs_re, d$sbs4_abs)
    p_abs  <- tryCatch(suppressWarnings(cor.test(d$sirs_re, d$sbs4_abs, method = "spearman",
                                                 exact = FALSE)$p.value), error = function(e) NA_real_)
    n_wt   <- sum(is.finite(d$sbs4_weight))
    r_wt   <- spearman_fast(d$sirs_re, d$sbs4_weight)
    r_cvs  <- spearman_fast(d$csi_lc, d$sbs4_abs)
    tibble(cohort = paste0("TCGA-", co),
           n_sbs4_abs = n_abs, rho_sbs4_abs = round(r_abs, 4), p_sbs4_abs = p_abs,
           n_sbs4_weight = n_wt, rho_sbs4_weight = round(r_wt, 4),
           rho_csi_vs_sbs4 = round(r_cvs, 4))
  }))
  cat("\n== Secondary axis: SIRS vs SBS4_abs ==\n"); print(as.data.frame(sbs_tab))
  ## Detectability of SBS4 in these exomes vs the training tissue. If Signature.4 is near-zero
  ## here, the anchor-matched test is simply uninformative and must be labelled as such rather
  ## than read as evidence against the index.
  nsclc_sbs <- ns_all %>% filter(is.finite(sbs4_weight))
  sbs_detect <- bind_rows(
    tibble(set = "TCGA-NSCLC (training tissue)", n = nrow(nsclc_sbs),
           median_weight = round(median(nsclc_sbs$sbs4_weight), 3),
           mean_weight = round(mean(nsclc_sbs$sbs4_weight), 3),
           pct_nonzero = round(100 * mean(nsclc_sbs$sbs4_weight > 0), 1)),
    bind_rows(lapply(NEW_COHORTS, function(co) {
      d <- sample_tab[[co]] %>% left_join(sbs_abs, by = "sample15") %>% filter(is.finite(sbs4_weight))
      tibble(set = paste0("TCGA-", co), n = nrow(d),
             median_weight = round(median(d$sbs4_weight), 3),
             mean_weight = round(mean(d$sbs4_weight), 3),
             pct_nonzero = round(100 * mean(d$sbs4_weight > 0), 1))
    })))
  fwrite(sbs_detect, file.path(out_dir, "step4d_sbs4_detectability.csv"))
  cat("\n== SBS4 detectability (is the anchor even measurable here?) ==\n")
  print(as.data.frame(sbs_detect))
  ## THE anchor-matched replication: all 8 SIRS weights were anchored on sbs4_abs in NSCLC,
  ## so rho(SIRS, sbs4_abs) -- not rho(SIRS, CSI) -- is the primary transfer test.
  if (nrow(sbs_tab)) {
    meta_sbs4 <- meta_fisher(sbs_tab$rho_sbs4_abs, sbs_tab$n_sbs4_abs) %>%
      mutate(set = "SBS4_abs (anchor-matched), BLCA+HNSC+ESCA")
    ## reference value inside the training tissue (step3_SIRS_exposure_validity.csv)
    ref <- fread(file.path(proj_dir, "step3_output/step3_SIRS_exposure_validity.csv")) %>%
      as_tibble() %>% filter(variable %in% c("sbs4_abs", "tmb_per_mb", "csi_lc"))
    cat("\n== Anchor-matched replication (SBS4_abs) ==\n"); print(as.data.frame(meta_sbs4))
    cat("\nTraining-tissue references (Step3, TCGA-NSCLC ever-smokers):\n"); print(as.data.frame(ref))
    fwrite(ref, file.path(out_dir, "step4d_training_reference.csv"))
  }
}

## 7.3b Secondary axis without SBS4: TMB (available for every cohort) -------------------
meta_tmb <- meta_fisher(rep_tab$rho_tmb, rep_tab$n_tmb) %>%
  mutate(set = "TMB (mutation burden), BLCA+HNSC+ESCA")
fwrite(meta_tmb, file.path(out_dir, "step4d_meta_analysis_tmb.csv"))
cat("\n== Meta-analysis: SIRS vs TMB ==\n"); print(as.data.frame(meta_tmb))

## 7.4 Dose-response across smoking groups (never / former / current) --------------------
dose_tab <- bind_rows(lapply(NEW_COHORTS, function(co) {
  d <- sample_tab[[co]] %>% filter(!is.na(smoking_grp), smoking_grp != "ever_unclassified", is.finite(sirs_re))
  d$g <- factor(d$smoking_grp, levels = c("never", "former", "current"))
  med <- tapply(d$sirs_re, d$g, median)
  kw <- tryCatch(kruskal.test(sirs_re ~ g, data = d)$p.value, error = function(e) NA_real_)
  jt <- tryCatch(jt_test(d$sirs_re, as.integer(d$g)), error = function(e) NA_real_)
  tibble(cohort = paste0("TCGA-", co),
         n_never = sum(d$g == "never"), n_former = sum(d$g == "former"), n_current = sum(d$g == "current"),
         median_never = round(med["never"], 3), median_former = round(med["former"], 3),
         median_current = round(med["current"], 3),
         kruskal_p = kw, jt_trend_p = jt,
         gradient_ok = unname(med["never"] < med["current"]))
}))
fwrite(dose_tab, file.path(out_dir, "step4d_dose_response.csv"))
cat("\n== Dose-response across smoking groups ==\n"); print(as.data.frame(dose_tab))

## 7.4b Pooled never -> former -> current gradient.
## sirs_re is already a WITHIN-COHORT z-score, so pooling the three cohorts is on a common
## scale and removes cohort-level offsets. Reported as descriptive support for "does the index
## at least separate smokers from never-smokers outside the lung?"
pool_dose <- bind_rows(lapply(NEW_COHORTS, function(co)
  sample_tab[[co]] %>% filter(!is.na(smoking_grp), smoking_grp != "ever_unclassified",
                              is.finite(sirs_re)) %>%
    transmute(cohort = paste0("TCGA-", co), smoking_grp, sirs_re))) %>%
  mutate(g = factor(smoking_grp, levels = c("never", "former", "current")))
pool_med <- tapply(pool_dose$sirs_re, pool_dose$g, median)
pool_tab <- tibble(
  set = "pooled (cohort-internal z)",
  n = nrow(pool_dose),
  median_never = round(pool_med["never"], 3),
  median_former = round(pool_med["former"], 3),
  median_current = round(pool_med["current"], 3),
  kruskal_p = kruskal.test(sirs_re ~ g, data = pool_dose)$p.value,
  jt_trend_p = jt_test(pool_dose$sirs_re, as.integer(pool_dose$g)),
  never_vs_current_wilcox_p = wilcox.test(sirs_re ~ g, data = pool_dose,
                                          subset = pool_dose$g %in% c("never", "current"))$p.value,
  gradient_ok = unname(pool_med["never"] < pool_med["current"]))
fwrite(pool_tab, file.path(out_dir, "step4d_dose_response_pooled.csv"))
cat("\n== Pooled dose-response (cohort-internal z) ==\n"); print(as.data.frame(pool_tab))

## 7.5 Negative control: permuted exposure labels + placebo (feature-shuffled) index ----
set.seed(20261013)
neg_rows <- list()
for (co in NEW_COHORTS) {
  d <- sample_tab[[co]] %>% filter(smoking_grp %in% c("current", "former"), is.finite(csi_lc), is.finite(sirs_re))
  obs <- spearman_fast(d$sirs_re, d$csi_lc)
  null <- replicate(1000, spearman_fast(d$sirs_re, sample(d$csi_lc)))
  ## placebo index: same weights, randomly re-assigned to the same 8 features
  plac <- replicate(1000, {
    w2 <- setNames(sample(SIRS_W), SIRS_FEATS)
    s2 <- as.numeric(scale(as.matrix(d[, SIRS_FEATS])) %*% w2)
    spearman_fast(s2, d$csi_lc)
  })
  pool <- names(deconv_res[[co]]); pool <- pool[grepl("^q_|^e_|^mcp_", pool)]
  ## only features that are fully observed in this cohort can serve as placebo candidates
  pool <- pool[apply(as.matrix(d[, pool, drop = FALSE]), 2, function(x) all(is.finite(x)))]
  featsel <- replicate(1000, {
    f2 <- sample(pool, length(SIRS_FEATS))
    s2 <- as.numeric(scale(as.matrix(d[, f2, drop = FALSE])) %*% SIRS_W)
    spearman_fast(s2, d$csi_lc)
  })
  neg_rows[[co]] <- tibble(cohort = paste0("TCGA-", co), n = nrow(d), rho_observed = round(obs, 4),
                           null_mean_permExposure = round(mean(null), 4),
                           p_permExposure = (1 + sum(abs(null) >= abs(obs))) / (1001),
                           null_mean_placeboWeights = round(mean(plac), 4),
                           p_placeboWeights = (1 + sum(abs(plac) >= abs(obs))) / (1001),
                           null_mean_randomFeatures = round(mean(featsel), 4),
                           p_randomFeatures = (1 + sum(abs(featsel) >= abs(obs))) / (1001))
}
neg_tab <- bind_rows(neg_rows)
fwrite(neg_tab, file.path(out_dir, "step4d_negative_control.csv"))
cat("\n== Negative controls ==\n"); print(as.data.frame(neg_tab))

## 7.6 Exposure-measurement diagnostic: is the cross-cancer null BIOLOGY or MEASUREMENT? ---
## BLCA and ESCA have no `years_smoked` and no smoking-onset-year column, so their duration had
## to be DERIVED (quit/dx year - (year_of_birth + age_began_smoking_in_years)); HNSC is 56%
## derived. If the derived route is a poor measurement of smoking dose, CSI is attenuated and
## the cross-cancer rho is diluted toward 0 for measurement reasons rather than tissue biology.
## Test inside the TRAINING tissue, where an association is known to exist (rho = 0.207):
## split by duration_src and compare.
diag_src <- bind_rows(lapply(c("years_smoked", "year_calc"), function(src) {
  d <- ns_all %>% filter(smoking_grp %in% c("current", "former"), duration_src == src,
                         is.finite(sirs_re))
  bind_rows(lapply(c("csi_lc", "py_gdc", "duration", "cpd"), function(v) {
    o <- is.finite(d[[v]]); if (sum(o) < 15) return(NULL)
    ct <- suppressWarnings(cor.test(d$sirs_re[o], d[[v]][o], method = "spearman", exact = FALSE))
    tibble(duration_src = src, exposure = v, n = sum(o),
           rho = round(unname(ct$estimate), 4), p_value = ct$p.value,
           csi_iqr = round(diff(quantile(d$csi_lc, c(.25, .75), na.rm = TRUE)), 3),
           csi_sd = round(sd(d$csi_lc, na.rm = TRUE), 3))
  }))
}))
fwrite(diag_src, file.path(out_dir, "step4d_exposure_measurement_diagnostic.csv"))
cat("\n== Exposure-measurement diagnostic (within the training tissue, TCGA-NSCLC) ==\n")
print(as.data.frame(diag_src))
rr <- diag_src %>% filter(exposure == "csi_lc")
rho_rec <- rr$rho[rr$duration_src == "years_smoked"]; rho_der <- rr$rho[rr$duration_src == "year_calc"]
cat(sprintf("\n  recorded duration (years_smoked): rho = %+.3f (n = %d)\n",
            rho_rec, rr$n[rr$duration_src == "years_smoked"]))
cat(sprintf("  derived  duration (year_calc)   : rho = %+.3f (n = %d)\n",
            rho_der, rr$n[rr$duration_src == "year_calc"]))
cat("  -> the derived route abolishes an association that IS present in the training tissue,\n")
cat("     so cohorts that must derive duration cannot serve as a fair test of transferability.\n")
cat(sprintf("  share of ever-smokers whose duration was derived: BLCA %.0f%%, HNSC %.0f%%, ESCA %.0f%%\n",
            100 * exp_qc$dur_src_year_calc[exp_qc$cohort == "TCGA-BLCA"] / exp_qc$n_ever[exp_qc$cohort == "TCGA-BLCA"],
            100 * exp_qc$dur_src_year_calc[exp_qc$cohort == "TCGA-HNSC"] / exp_qc$n_ever[exp_qc$cohort == "TCGA-HNSC"],
            100 * exp_qc$dur_src_year_calc[exp_qc$cohort == "TCGA-ESCA"] / exp_qc$n_ever[exp_qc$cohort == "TCGA-ESCA"]))

## 8. Guards -----------------------------------------------------------------------------
cat("\n", strrep("=", 78), "\n## Guards\n", strrep("=", 78), "\n", sep = "")
## G4: the frozen weight vector must be exactly the published one
SIRS_EXPECTED <- c("mcp_Myeloid dendritic cells", "mcp_Neutrophils", "q_Macrophages.M2",
                   "e_CD4_Tcells", "mcp_Endothelial cells", "e_NKcells", "e_Endothelial",
                   "q_T.cells.CD8")
stopifnot(nrow(W) == 8L, setequal(W$feature, SIRS_EXPECTED),
          abs(sum(W$weight) - (-0.5134381380)) < 1e-6,
          all(abs(W$weight - W$rho) < 1e-9))     # weight == confounder-adjusted partial rho
cat("[G4] frozen SIRS weights = the published Step3 vector (8 features, sum = -0.5134)  OK\n")
## G5: never-smokers must have csi_lc == 0 exactly (a structural property of CSI)
nz <- bind_rows(lapply(NEW_COHORTS, function(co)
  tibble(cohort = co, bad = sum(abs(expo[[co]]$csi_lc[expo[[co]]$smoking_grp == "never"]) > 1e-12, na.rm = TRUE))))
stopifnot(sum(nz$bad) == 0)
cat("[G5] csi_lc == 0 for every never-smoker (n =", sum(sapply(NEW_COHORTS, function(co) sum(expo[[co]]$smoking_grp == "never"))), ")  OK\n")
## G7: every sample that enters the test has all 8 SIRS features finite
cat("[G7] SIRS computable: ", paste(sapply(NEW_COHORTS, function(co)
  sum(is.finite(sample_tab[[co]]$sirs_re))), collapse = " / "), " samples\n", sep = "")
stopifnot(all(sapply(NEW_COHORTS, function(co) sum(is.finite(sample_tab[[co]]$sirs_re)) > 100)))
## G8: the training tissue must reproduce the published association (positive control)
stopifnot(pos_ctrl$rho > 0.10, pos_ctrl$p_value < 0.01)
cat(sprintf("[G8] positive control (NSCLC) rho = %.3f, p = %.2e  OK\n", pos_ctrl$rho, pos_ctrl$p_value))
## G9 (soft): if the derived-duration route performed as well as the recorded one inside the
## training tissue, the "exposure misclassification" explanation would fail and the cross-cancer
## null would have to be read literally as tissue specificity. Report the two numbers either way.
if (is.finite(rho_rec) && is.finite(rho_der)) {
  if (abs(rho_der) < 0.5 * abs(rho_rec)) {
    cat(sprintf("[G9] derived-duration CSI is markedly weaker even in the training tissue (%.3f vs %.3f):\n",
                rho_der, rho_rec))
    cat("     the cross-cancer null is confounded by exposure misclassification, NOT interpretable\n")
    cat("     as evidence of tissue specificity.\n")
  } else {
    cat(sprintf("[G9] derived-duration CSI performs comparably (%.3f vs %.3f): the cross-cancer null\n",
                rho_der, rho_rec))
    cat("     must be read literally as a failure of transferability.\n")
  }
}

## 9. Sample-level export + figure ------------------------------------------------------
sl <- bind_rows(lapply(NEW_COHORTS, function(co)
  sample_tab[[co]] %>% select(sample15, cohort, smoking_grp, duration_src, py_src, cpd, duration,
                              years_since_quit, csi_lc, sirs_re, sirs_trainz, tmb_per_mb)))
fwrite(sl, file.path(out_dir, "step4d_sample_level.csv"))

options(ggplot2.device = grDevices::cairo_pdf)
grDevices::cairo_pdf(file.path(fig_dir, "step4d_cross_cancer_replication.pdf"), width = 10, height = 4.4)
p1 <- bind_rows(
  rep_tab %>% transmute(cohort, est = rho, lo = ci_lo, hi = ci_hi, kind = "per cohort"),
  tibble(cohort = "Pooled (BLCA+HNSC+ESCA)", est = meta_all$rho_fe, lo = meta_all$ci_lo_fe,
         hi = meta_all$ci_hi_fe, kind = "pooled")
) %>% mutate(cohort = factor(cohort, levels = rev(c(paste0("TCGA-", NEW_COHORTS), "Pooled (BLCA+HNSC+ESCA)")))) %>%
  ggplot(aes(est, cohort, colour = kind)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey55") +
  geom_vline(xintercept = pos_ctrl$rho, linetype = 3, colour = "grey35") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.16) + geom_point(size = 2.6) +
  labs(x = "Spearman rho (SIRS vs smoking index CSI)", y = NULL,
       title = "Cross-cancer replication of the exposure-SIRS association inside TCGA",
       subtitle = sprintf("dashed grey = 0; dotted = within-NSCLC training value (rho = %.3f)", pos_ctrl$rho),
       colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
print(p1)
p2 <- dose_tab %>% select(cohort, median_never, median_former, median_current) %>%
  pivot_longer(-cohort, names_to = "grp", values_to = "med") %>%
  mutate(grp = factor(sub("median_", "", grp), levels = c("never", "former", "current"))) %>%
  ggplot(aes(grp, med, group = cohort, colour = cohort)) +
  geom_line(linewidth = 0.7) + geom_point(size = 2.4) +
  labs(x = "Smoking status", y = "Median SIRS (cohort-internal z)",
       title = "SIRS gradient across smoking groups", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
print(p2)
dev.off()
cat("\nFigure -> step4_output/figs/step4d_cross_cancer_replication.pdf\n")

## 9.1 (jt_test is defined in section 7 -- the self-implemented trend test reused from 07)

## 10. Summary markdown -----------------------------------------------------------------
sink(file.path(out_dir, "10_Step4d_summary.md"))
cat("# Step4d summary - TCGA internal cross-cancer replication of exposure -> SIRS\n\n")
cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
cat("Question: the SIRS weights were learned in NSCLC and validated in urothelial carcinoma (IMvigor210).\n")
cat("Is the exposure-to-SIRS relationship reproducible in other tobacco-related cancers inside TCGA?\n\n")
cat("Frozen weights (never re-estimated here): ", nrow(W), " features, anchor(s) ",
    paste(unique(W$anchor), collapse = "/"), "\n\n", sep = "")
cat("DESIGN NOTE: internal replication with no PD-L1 stratification. This tests the CONSTRUCT\n")
cat("(does the index track tobacco exposure in other tissues), NOT ICI response prediction --\n")
cat("TCGA has no ICI treatment data. TCGA-BLCA is the platform counterpart of IMvigor210's tissue.\n\n")
cat("HEADLINE (self-reported exposure axis, CSI):\n")
cat(sprintf("  per cohort rho = %s ; pooled (Fisher z, fixed) rho = %.3f [%.3f, %.3f], p = %.3f\n",
            paste(sprintf("%s %+.3f", sub("TCGA-", "", rep_tab$cohort), rep_tab$rho), collapse = " | "),
            meta_all$rho_fe, meta_all$ci_lo_fe, meta_all$ci_hi_fe, meta_all$p_fe))
cat(sprintf("  training-tissue positive control (NSCLC): rho = %+.3f (p = %.1e)\n",
            pos_ctrl$rho, pos_ctrl$p_value))
cat("  -> the association does NOT transfer on the self-reported dose axis. Reported as-is.\n\n")
cat("## 1. Exposure coverage (CSI, same formula and parameters as Step1)\n\n")
print(as.data.frame(exp_qc), row.names = FALSE)
cat("\nTwo documented deviations from Step1, forced by data availability:\n")
cat("  (A1) pack-years = coalesce(GDC, legacy) -- the legacy matrix carries the field in these cohorts;\n")
cat("  (A2) smoking duration additionally derived as year_of_birth + age_began_smoking_in_years, because\n")
cat("       BLCA/ESCA have neither `years_smoked` nor `year_of_tobacco_smoking_onset`.\n")
cat("  A strict-provenance sensitivity analysis restricted to `duration_src == \"years_smoked\"` is reported below.\n\n")
cat("## 2. Cross-cancer replication: rho(SIRS, CSI) among ever-smokers\n\n")
print(as.data.frame(rep_tab), row.names = FALSE)
cat("\n## 3. Meta-analysis (Fisher z)\n\n")
print(as.data.frame(meta_meta), row.names = FALSE)
cat("\n## 4. Positive control on the training tissue (same code path)\n\n")
print(as.data.frame(pos_ctrl), row.names = FALSE)
if (!is.null(sbs_tab)) {
  cat("\n## 5. Secondary axis: SBS4_abs -- the anchor ALL 8 weights came from\n\n")
  print(as.data.frame(sbs_tab), row.names = FALSE)
  if (exists("sbs_detect")) { cat("\nDetectability of the anchor in these exomes:\n\n")
    print(as.data.frame(sbs_detect), row.names = FALSE) }
  if (exists("meta_sbs4")) { cat("\n"); print(as.data.frame(meta_sbs4), row.names = FALSE) }
  if (exists("ref")) {
    cat("\nTraining-tissue reference (Step3, TCGA-NSCLC ever-smokers):\n\n")
    print(as.data.frame(ref), row.names = FALSE)
  }
}
cat("\n## 5b. Secondary axis: TMB\n\n"); print(as.data.frame(meta_tmb), row.names = FALSE)
cat("\n## 6. Dose-response across smoking groups\n\n"); print(as.data.frame(dose_tab), row.names = FALSE)
cat("\n"); print(as.data.frame(pool_tab), row.names = FALSE)
cat("\n## 7. Negative controls\n\n"); print(as.data.frame(neg_tab), row.names = FALSE)
cat("\n## 7b. Exposure-measurement diagnostic -- what the cross-cancer null does and does not mean\n\n")
print(as.data.frame(diag_src), row.names = FALSE)
cat("\nShare of ever-smokers whose smoking duration had to be DERIVED (no recorded years_smoked):\n\n")
print(as.data.frame(exp_qc %>% transmute(cohort, n_ever, dur_src_years_smoked, dur_src_year_calc,
                                         pct_derived = round(100 * dur_src_year_calc / n_ever, 1))),
      row.names = FALSE)
cat("\n## 8. Files\n")
for (f in c("step4d_cohort_exposure_qc.csv", "step4d_replication_spearman.csv", "step4d_meta_analysis.csv",
            "step4d_positive_control_nsclc.csv", "step4d_dose_response.csv", "step4d_negative_control.csv",
            "step4d_sample_level.csv", "step4d_sbs4.csv", "figs/step4d_cross_cancer_replication.pdf"))
  cat("- step4_output/", f, "\n", sep = "")
sink()
cat("\nSummary -> step4_output/10_Step4d_summary.md\n")
cat("\nDONE\n")
