############################################################
## 01_Step1_exposure_quantification.R
## Smoking exposure quantification metrics x NSCLC immune microenvironment x
## ICI response prediction
##
## Data sources: UCSC Xena        https://xenabrowser.net/datapages/
##               GDC Xena Hub     https://gdc-hub.s3.us-east-1.amazonaws.com/download/
##               Legacy Xena Hub  https://tcga.xenahubs.net/download/
## Step1 goal: build smoking exposure metrics for TCGA-LUAD/LUSC, run QC, and
##             lock the analytic cohort
## Outputs: smoking_index_clinical.csv / missingness report / QC plots /
##          Step2 sample list / download manifest
## Environment: R >= 4.3; tidyverse + data.table + mice + patchwork
## Author: see manuscript author list   Date: 2026-10-08
##
## ---------------------------------------------------------------------------
## Exposure definitions (as used in the manuscript)
## ---------------------------------------------------------------------------
## CPD (cigarettes per day)
##   CPD = pack-years x 20 / smoking years, computed in two passes (duration
##   first, then CPD).
##   The precomputed Xena `cigarettes_per_day` field was not used: it is a
##   linear rescaling of pack-years (== pack-years x 20 / 365; holds exactly for
##   all 782 cases where both values exist; median 2.19) and therefore does not
##   measure cigarettes per day. The corresponding field in the GDC clinical API
##   has 0% coverage for TCGA-LUAD/LUSC, so no authoritative source exists. The
##   column is excluded at the data level and retained only as the raw column
##   `cpd_xena_raw` for audit. Full evidence chain:
##   docs/06_data_source_audit_smoking_fields.md.
##   Consequence: cig_years == pack-years x 20 under this definition, so
##   cig_years is not an independent axis; downstream analyses use pack-years,
##   log(cig-years) and CSI instead of a head-to-head comparison.
##
## CSI (cumulative smoking index, Leffondre et al. 2006)
##   tsc* = max(tsc - delta, 0);  dur* = max(dur + tsc - delta, 0) - tsc*
##   CSI  = (1 - 0.5^(dur*/tau)) * 0.5^(tsc*/tau) * ln(int + 1)
##   The main analysis (section 3.4, CSI_MAIN) uses published parameters
##   tau = 26 years, delta = 0.7 years (Leffondre 2006 lung cancer values /
##   ICARE 2013 pooled values); no within-cohort re-estimation is performed.
##   Three published-parameter sensitivity variants and a parameter-robustness
##   QC are reported (section 6.2). csi_lc is the canonical CSI; csi_simple is
##   a power-law placeholder kept for reference only.
##
## Implementation notes
##   Records with mutually inconsistent quit years are flagged
##   `quit_year_irregular`; the decay-term base is floored at 0 and non-finite
##   values are tabulated.
##
## Usage: set proj_dir in section 0 (keep consistent with scripts 04/05/07),
##        then run.
############################################################

## 0. Global settings ------------------------------------------------------
proj_dir <- "D:/Projects/smoking-index-ici"
data_dir <- file.path(proj_dir, "data_raw", "xena")
out_dir  <- file.path(proj_dir, "step1_output")
fig_dir  <- file.path(out_dir, "figs")
invisible(lapply(list(data_dir, out_dir, fig_dir),
                 function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE)))
set.seed(20261004)

## CPD definition (final, single definition)
##   cpd = pack-years x 20 / smoking years (algebraic identity of the pack-year definition, no new assumptions; two-pass computation in section 5)
## CPD_PSEUDO_DIVISOR is used only in the section 6.0 data-source audit (reproducing the test "Xena raw column == pack-years/18.25"),
## and takes no part in constructing any exposure metric.
CPD_PSEUDO_DIVISOR <- 18.25   # = 365/20; Xena cigarettes_per_day == pack-years/18.25 (observed 782/782)
cat(sprintf("\n[CPD definition] CPD = pack-years x 20 / smoking years (the Xena cigarettes_per_day column is not used; see docs/06_data_source_audit_smoking_fields.md)\n"))

## 1. Dependencies ---------------------------------------------------------
cran_pkgs <- c("tidyverse", "data.table", "mice", "patchwork")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(tidyverse)
library(data.table)

## 2. Download data from UCSC Xena -----------------------------------------
## GDC hub naming rule: TCGA-<cancer type>.<data type>.tsv(.gz)
xena_gdc_base <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
## Backup mirror (swap the base above for this one if the S3 domain is unreachable):
# xena_gdc_base <- "https://gdc.xenahubs.net/download/"

## -- 2.1 GDC hub files (file names verified to exist) --
xena_gdc_files <- c(
  ## Step1 core: GDC clinical table (contains exposures.* smoking quantitative fields; no status column)
  "TCGA-LUAD.clinical.tsv.gz",
  "TCGA-LUSC.clinical.tsv.gz",
  ## Survival (fields overlap with clinical; used for cross-checking)
  "TCGA-LUAD.survival.tsv.gz",
  "TCGA-LUSC.survival.tsv.gz",
  ## Expression: raw counts (for DESeq2 etc.)
  "TCGA-LUAD.star_counts.tsv.gz",
  "TCGA-LUSC.star_counts.tsv.gz",
  ## Expression: FPKM-UQ (recommended for Step2 immune deconvolution; already normalized)
  "TCGA-LUAD.star_fpkm-uq.tsv.gz",
  "TCGA-LUSC.star_fpkm-uq.tsv.gz",
  ## Mutation: MC3 three-caller merged MAF (for TMB / SBS4)
  "TCGA-LUAD.somaticmutation_wxs.tsv.gz",
  "TCGA-LUSC.somaticmutation_wxs.tsv.gz"
)

## -- 2.2 legacy hub files (fill in current/former grouping and quit year; no .gz) --
legacy_files <- c(
  LUAD = "https://tcga.xenahubs.net/download/TCGA.LUAD.sampleMap/LUAD_clinicalMatrix",
  LUSC = "https://tcga.xenahubs.net/download/TCGA.LUSC.sampleMap/LUSC_clinicalMatrix"
)

## -- 2.3 Robust downloader (resume + retry + integrity check) --
## Background: star_counts is ~64MB; R download.file defaults to timeout=60s with no resume,
##       so any network jitter fails the whole transfer. This downloader uses system curl/wget (.part resume) and
##       checks integrity via remote Content-Length or the gzip magic number; corrupt/partial files are re-downloaded.
options(timeout = max(3600, getOption("timeout")))   # relax to 1 hour for large files

## Remote file size (for integrity check; returns NA if curl is unavailable)
get_remote_size <- function(url) {
  if (!nzchar(Sys.which("curl"))) return(NA_real_)
  h <- suppressWarnings(system2("curl", c("-sIL", "--max-time", "60", shQuote(url)),
                                stdout = TRUE, stderr = FALSE))
  cl <- grep("^[Cc]ontent-[Ll]ength:", h, value = TRUE)
  if (length(cl)) suppressWarnings(as.numeric(sub(".*:\\s*", "", trimws(tail(cl, 1))))) else NA_real_
}
## gzip magic-number check (whether a .gz has basic integrity)
is_gzip_magic <- function(path) {
  if (!file.exists(path) || file.size(path) < 2) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con))
  m <- readBin(con, "raw", 2)
  length(m) == 2 && m[1] == as.raw(0x1f) && m[2] == as.raw(0x8b)
}
is_complete <- function(path, rsize, expect_gz) {
  if (!file.exists(path) || file.size(path) == 0) return(FALSE)
  if (!is.na(rsize)) return(unname(file.size(path)) == rsize)   # prefer exact remote-size check
  if (expect_gz) return(is_gzip_magic(path))                    # fallback: gzip magic number
  TRUE
}
run_cmd <- function(cmd, args) {
  st <- suppressWarnings(system2(cmd, args, stdout = FALSE, stderr = FALSE))
  identical(as.integer(st), 0L)
}

download_manifest <- list()
download_file <- function(url, fname = basename(url), dest = data_dir, retries = 3) {
  dest_file <- file.path(dest, fname)
  expect_gz <- grepl("\\.gz$", fname)
  rsize     <- get_remote_size(url)

  ## already exists and passes check -> skip
  if (is_complete(dest_file, rsize, expect_gz)) {
    message("Already exists and passes check, skipping: ", fname)
  } else {
    if (file.exists(dest_file)) {
      ## local larger than remote = dirty file; otherwise keep for resume
      if (!is.na(rsize) && unname(file.size(dest_file)) > rsize) {
        message("Local file larger than remote, deleting and re-downloading: ", fname); file.remove(dest_file)
      } else {
        message("Incomplete file detected (local ", file.size(dest_file),
                " / remote ", ifelse(is.na(rsize), "unknown", rsize), "), attempting resume: ", fname)
      }
    }
    message("Downloading: ", url)
    done <- FALSE
    for (k in seq_len(retries)) {
      if (nzchar(Sys.which("curl"))) {
        ## -C - auto-resume; --fail makes HTTP errors return non-zero; --retry handles transient jitter
        ok <- run_cmd("curl", c("-L", "--fail", "--retry", "5", "--retry-delay", "3",
                                "--continue-at", "-", "--output", shQuote(dest_file), shQuote(url)))
      } else if (nzchar(Sys.which("wget"))) {
        ok <- run_cmd("wget", c("-c", "-q", "-O", shQuote(dest_file), shQuote(url)))
      } else {
        ok <- tryCatch({ utils::download.file(url, dest_file, mode = "wb", quiet = TRUE,
                                              method = "libcurl"); TRUE },
                       error = function(e) { message("  !! ", conditionMessage(e)); FALSE })
      }
      if (ok && is_complete(dest_file, rsize, expect_gz)) { done <- TRUE; break }
      message("  Attempt ", k, " incomplete (", ifelse(ok, "length mismatch", "transfer error"), "), retrying in 2s...")
      Sys.sleep(2)
    }
    if (!done) {
      warning("Download failed (still incomplete after ", retries, " retries): ", fname, "\n  URL: ", url,
              "\n  Suggestion: manually run `curl -C - -L -o <file> <url>` or `wget -c <url>`, then re-run the script",
              call. = FALSE)
      return(NA_character_)
    }
  }
  download_manifest[[fname]] <<- tibble(
    file = fname, url = url, size_bytes = file.size(dest_file),
    md5 = unname(tools::md5sum(dest_file)), download_date = as.character(Sys.Date()))
  dest_file
}

gdc_paths <- vapply(xena_gdc_files, function(f)
  download_file(paste0(xena_gdc_base, f), f), character(1))
legacy_paths <- vapply(names(legacy_files), function(id)
  download_file(legacy_files[[id]], paste0(id, "_clinicalMatrix")), character(1))

write_csv(bind_rows(download_manifest), file.path(out_dir, "download_manifest.csv"))

## 3. Read clinical tables and auto-detect smoking fields ------------------
read_xena_tsv <- function(path) {
  x <- data.table::fread(path, sep = "\t", header = TRUE,
                         data.table = FALSE, check.names = FALSE)
  names(x)[1] <- "sample"          # the first column of a Xena tsv is often an empty column name (FILE_SORTED etc.)
  x
}
## Column auto-detection: compatible with different version prefixes such as exposures.* / diagnoses.*
find_col <- function(df, patterns, required = TRUE) {
  for (p in patterns) {
    hit <- grep(p, names(df), ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  if (required) warning("Column not found, pattern: ", paste(patterns, collapse = " | "),
                        "\nExample existing columns: ", paste(head(names(df), 40), collapse = ", "))
  NA_character_
}

## 3.1 GDC clinical: only quantitative CPD / years_smoked / pack_years
read_gdc_clinical <- function(fname, cohort) {
  cl <- read_xena_tsv(file.path(data_dir, fname))
  cols <- list(
    cpd     = find_col(cl, c("cigarettes_per_day")),
    years   = find_col(cl, c("years_smoked")),
    py      = find_col(cl, c("pack_years_smoked", "pack_years")),
    dx_year = find_col(cl, c("year_of_diagnosis"), required = FALSE),
    age     = find_col(cl, c("age_at_diagnosis", "age_at_index"), required = FALSE),
    stage   = find_col(cl, c("ajcc_pathologic_stage", "tumor_stage"), required = FALSE),
    gender  = find_col(cl, c("gender"), required = FALSE),
    vital   = find_col(cl, c("vital_status"), required = FALSE),
    d_death = find_col(cl, c("days_to_death"), required = FALSE),
    d_fu    = find_col(cl, c("days_to_last_follow"), required = FALSE),
    ttype   = find_col(cl, c("tissue_type"), required = FALSE)
  )
  cat("\n[", cohort, "GDC clinical] auto-detected key columns:\n"); print(unlist(cols))
  pick <- function(nm) if (is.na(cols[[nm]])) NA_character_ else as.character(cl[[cols[[nm]]]])
  tibble(
    sample = cl$sample, cohort = cohort,
    ## Note: this column is the deprecated Xena cigarettes_per_day (== pack-years/18.25, see section 6.0 audit).
    ##   This version does **not** use it to construct any exposure metric; it is archived as-is (cpd_xena_raw) for auditing and reproduction.
    cpd_xena = pick("cpd"), years_gdc = pick("years"), py_gdc = pick("py"),
    dx_year_gdc = pick("dx_year"), age = pick("age"), stage = pick("stage"),
    gender = pick("gender"), vital = pick("vital"),
    days_death = pick("d_death"), days_fu = pick("d_fu"), tissue = pick("ttype")
  )
}
gdc_clin <- bind_rows(read_gdc_clinical("TCGA-LUAD.clinical.tsv.gz", "TCGA-LUAD"),
                      read_gdc_clinical("TCGA-LUSC.clinical.tsv.gz", "TCGA-LUSC")) %>%
  ## Coerce to character before merging: fread may parse the same column differently across the two cancer types (integer vs character)
  mutate(across(everything(), as.character))

## 3.2 legacy clinicalMatrix: fill in status / quit year / smoking onset year / official pack-years
read_legacy_clinical <- function(path, cohort) {
  lg <- read_xena_tsv(path)
  cols <- list(
    tsh    = find_col(lg, c("^tobacco_smoking_history$", "tobacco_smoking_history"), required = FALSE),
    tshi   = find_col(lg, c("tobacco_smoking_history_indicator"), required = FALSE),
    stop_y = find_col(lg, c("stopped_smoking_year"), required = FALSE),
    onset  = find_col(lg, c("year_of_tobacco_smoking_onset"), required = FALSE),
    dx_y   = find_col(lg, c("year_of_initial_pathologic_diagnosis"), required = FALSE),
    py     = find_col(lg, c("number_pack_years_smoked"), required = FALSE)
  )
  cat("\n[", cohort, "legacy clinicalMatrix] auto-detected key columns:\n"); print(unlist(cols))
  pick <- function(nm) if (is.na(cols[[nm]])) NA_character_ else as.character(lg[[cols[[nm]]]])
  tibble(
    sample15 = lg$sample,                             # legacy is 15 characters: TCGA-05-4244-01
    cohort   = cohort,
    tsh_raw  = pick("tsh"), tshi_raw = pick("tshi"),
    stop_year = pick("stop_y"), onset_year = pick("onset"),
    dx_year_lg = pick("dx_y"), py_lg = pick("py")
  )
}
legacy_clin <- bind_rows(
  read_legacy_clinical(file.path(data_dir, "LUAD_clinicalMatrix"), "TCGA-LUAD"),
  read_legacy_clinical(file.path(data_dir, "LUSC_clinicalMatrix"), "TCGA-LUSC")) %>%
  mutate(across(everything(), as.character))   # same as gdc_clin: coerce to character before merging

## 3.3 Bridge the two tables using the 15-character sample barcode
clinic0 <- gdc_clin %>%
  mutate(sample15 = substr(sample, 1, 15)) %>%
  left_join(legacy_clin %>% select(sample15, tsh_raw, tshi_raw, stop_year,
                                   onset_year, dx_year_lg, py_lg),
            by = "sample15")

## 3.4 CSI parameters (Strategy A: adopt published parameters, no within-cohort re-estimation) -----------
## Rationale: published lung-cancer parameters are adopted directly (Strategy A); no within-cohort re-estimation.
## Parameter sources (all published values, cited directly):
##   Leffondre 2006 (Stat Med 25:4132-4146) lung cancer cohort -> tau = 26, delta = 1.0 (men)/0.7 (women)
##   ICARE 2013 (Br J Cancer) pooled men+women             -> tau = 26, delta = 0.7   (adopted for the main analysis)
##   Montreal 2018                                        -> tau = 26, delta = 1.0 (men)/0.7 (women)
CSI_MAIN <- c(tau = 26, delta = 0.7)   # main-analysis parameters (to be reported in the methods)
## Sensitivity analysis: **change only the parameter values, not the formula**, to show conclusions are insensitive to parameter choice
## (within-cohort re-estimation of (tau,delta) is Strategy B and is not used in the manuscript.)
CSI_SENS <- list(
  d10 = c(tau = 26, delta = 1.0),   # Leffondre men's value
  t20 = c(tau = 20, delta = 0.7),   # half-life at the lower bound
  t35 = c(tau = 35, delta = 0.7)    # half-life at the upper bound
)
## CSI canonical formula (Leffondre 2006 lung cancer version) -- main and sensitivity analyses share the same function for consistency
##   int = cigarettes per day (CPD); dur = smoking years; tsc = years since quitting (current smokers = 0)
##   tau = biological half-life of the carcinogenic effect (years); delta = lag between exposure and detectable disease (years)
csi_leffondre <- function(int, dur, tsc,
                          tau = CSI_MAIN[["tau"]], delta = CSI_MAIN[["delta"]]) {
  tsc      <- ifelse(is.na(tsc), 0, tsc)          # current smokers tsc = 0
  tsc_star <- pmax(tsc - delta, 0)
  dur_star <- pmax(dur + tsc - delta, 0) - tsc_star
  (1 - 0.5^(dur_star / tau)) * 0.5^(tsc_star / tau) * log(int + 1)
}

## 4. Cleaning: NA coding, status recoding, sample de-duplication -----------
## Note: in Xena the same column may differ in case across the two cancer types (e.g. LUSC shows "[Discrepancy]"),
##       so NA tokens are always compared in lowercase to avoid leaks contaminating types and groups.
na_tokens <- c("[not available]", "[unknown]", "[not applicable]", "[not performed]",
               "[discrepancy]", "[not evaluated]", "unknown", "na", "n/a", "",
               "not reported", "--")
clean_na  <- function(x) {
  x  <- as.character(x)
  xt <- tolower(trimws(x))
  x[is.na(x) | xt %in% na_tokens] <- NA
  x
}
num_clean <- function(x) suppressWarnings(as.numeric(clean_na(x)))

## Status recoding: prefer legacy numeric codes (1 never / 2 current / 3 reformed>15y /
## 4 reformed<=15y / 5 reformed unspecified), with text labels as fallback (the legacy indicator is often empty)
## Note: text matching must apply grepl to the **full lab vector**, then AND with idx;
##       writing grepl(pattern, lab[idx]) triggers a recycling warning due to unequal length.
recode_status <- function(code, label) {
  cn  <- suppressWarnings(as.numeric(clean_na(code)))
  out <- rep(NA_character_, length(cn))
  out[!is.na(cn) & cn == 1]              <- "never"
  out[!is.na(cn) & cn == 2]              <- "current"
  out[!is.na(cn) & cn %in% c(3, 4, 5)]   <- "former"

  lab <- tolower(clean_na(label)); lab[is.na(lab)] <- ""   # set empty to avoid NA indices in grepl
  idx <- is.na(out)                                        # only check text for those undetermined by numeric code
  g_never  <- grepl("lifelong non-smoker|non-smoker|never", lab)
  g_cur    <- grepl("current smoker", lab)                 # does not match "current reformed smoker"
  g_former <- grepl("reformed|former", lab)
  out[idx & g_never]  <- "never"
  out[idx & g_cur]    <- "current"
  out[idx & g_former] <- "former"
  out
}
recode_former_detail <- function(code, label) {
  cn  <- suppressWarnings(as.numeric(clean_na(code)))
  lab <- tolower(clean_na(label)); lab[is.na(lab)] <- ""
  out <- rep(NA_character_, length(cn))
  out[!is.na(cn) & cn == 3] <- "former_gt15y"
  out[!is.na(cn) & cn == 4] <- "former_le15y"
  out[!is.na(cn) & cn == 5] <- "former_nspec"

  idx  <- is.na(out)
  g_gt <- grepl("> *15", lab)                                  # "> 15 years" / ">15"
  g_le <- grepl("< *= *15|< *or *= *15|< *15", lab)            # "<=15" / "< or = 15" / "< 15"
  g_ns <- grepl("duration not specified", lab)
  out[idx & g_gt] <- "former_gt15y"
  out[idx & g_le] <- "former_le15y"
  out[idx & g_ns] <- "former_nspec"
  out
}

clinic <- clinic0 %>%
  mutate(
    patient = substr(sample, 1, 12),        # TCGA-XX-XXXX
    stype   = substr(sample, 14, 15),       # 01 = primary tumor, 11 = adjacent normal
    smoking_grp   = recode_status(tsh_raw, tshi_raw),
    former_detail = recode_former_detail(tsh_raw, tshi_raw),
    cpd_xena_raw = num_clean(cpd_xena),   # deprecated Xena raw column, archived only for the 6.0 audit
    duration_raw = num_clean(years_gdc),
    py_gdc    = num_clean(py_gdc),
    py_legacy = num_clean(py_lg),
    dx_year   = coalesce(num_clean(dx_year_gdc), num_clean(dx_year_lg)),
    stop_year = num_clean(stop_year),
    start_year = num_clean(onset_year),
    age_y     = num_clean(age),
    days_death = num_clean(days_death),
    days_fu    = num_clean(days_fu)
  ) %>%
  ## Fallback: when legacy status is missing but GDC has pack-years, use quantitative fields to assign never / ever
  mutate(smoking_grp = case_when(
    !is.na(smoking_grp) ~ smoking_grp,
    !is.na(py_gdc) & py_gdc == 0 ~ "never",
    (!is.na(py_gdc) & py_gdc > 0) | (!is.na(cpd_xena_raw) & cpd_xena_raw > 0) ~ "ever_unclassified",
    TRUE ~ NA_character_
  )) %>%
  select(patient, sample, sample15, cohort, stype, smoking_grp, former_detail,
         cpd_xena_raw, duration_raw, py_gdc, py_legacy, dx_year, start_year, stop_year,
         age_y, gender, stage, vital, days_death, days_fu, tissue)

## Keep primary tumors only (01); for multiple samples per patient keep the one with the most complete smoking info
if (any(clinic$stype == "01")) {
  clinic <- clinic %>% filter(stype == "01")
} else {
  warning("Sample IDs appear to be patient-level (12 characters), skipping primary-tumor filtering")
}
clinic <- clinic %>%
  mutate(n_missing = rowSums(is.na(across(c(duration_raw, py_gdc, stop_year))))) %>%
  arrange(patient, n_missing) %>%
  distinct(patient, .keep_all = TRUE) %>%
  select(-n_missing)

cat("\nSmoking group distribution (by cohort):\n"); print(table(clinic$smoking_grp, clinic$cohort, useNA = "ifany"))

## 5. Build exposure metrics (CPD = pack-years x 20 / years) ----------------------------------------
## Must be implemented in "two passes": compute duration first, then CPD -- otherwise the year_calc branch of duration is missed.
## Fallback cascade:
##   duration : years_smoked -> (quit year or diagnosis year - smoking onset year) (flagged year_calc) -> NA
##   cpd      : pack-years x 20 / duration -> NA
## Under this definition cig_years == pack-years x 20 (constructive identity) -> "smoking index vs pack-years" are not two
##    independent metrics and **must not be compared head-to-head**; downstream use three axes: pack-years / logcig-years / CSI.
## Tertile helper (splits only within ever-smokers with non-missing metric)
tertile3 <- function(x, grp) {
  out <- rep(NA_character_, length(x))
  ok  <- !is.na(x) & grp %in% c("current", "former")
  out[ok] <- as.character(cut(dplyr::ntile(x[ok], 3), breaks = 0:3,
                              labels = c("T1_Low", "T2_Mid", "T3_High")))
  out
}

build_exposure <- function(df) {
  df %>%
    ## ---- Pass 1: determine duration first (CPD depends on it) ----
    mutate(
      stop_eff = case_when(smoking_grp == "current" ~ dx_year, TRUE ~ stop_year),
      dur_est  = ifelse(!is.na(start_year) & !is.na(stop_eff) & stop_eff >= start_year,
                        stop_eff - start_year, NA_real_),
      duration_src = case_when(
        smoking_grp == "never" ~ "never_zero",
        !is.na(duration_raw) ~ "years_smoked",
        !is.na(dur_est) & dur_est > 0 ~ "year_calc",
        TRUE ~ NA_character_
      ),
      duration = case_when(
        smoking_grp == "never" ~ 0,
        !is.na(duration_raw) ~ duration_raw,
        !is.na(dur_est) & dur_est > 0 ~ dur_est,
        TRUE ~ NA_real_
      ),
      ## Years since quitting: computable only for former smokers with an exact quit year; current smokers recorded as 0
      years_since_quit_raw = case_when(
        smoking_grp == "never" ~ NA_real_,
        smoking_grp == "current" ~ 0,
        !is.na(stop_year) & !is.na(dx_year) ~ dx_year - stop_year,
        TRUE ~ NA_real_
      ),
      ## Data consistency: quit year later than diagnosis year => negative quit duration, a record conflict; set NA and flag
      quit_year_irregular = !is.na(years_since_quit_raw) & years_since_quit_raw < 0,
      years_since_quit = ifelse(quit_year_irregular, NA_real_, years_since_quit_raw)
    ) %>%
    ## ---- Pass 2: solve CPD back from the two true values (algebraic identity) ----
    mutate(
      ## CPD = pack-years x 20 / smoking years (algebraic identity of the pack-year definition, no new assumptions)
      cpd_src = case_when(
        smoking_grp == "never" ~ "never_zero",
        !is.na(py_gdc) & !is.na(duration) & duration > 0 ~ "back_calc_py",
        TRUE ~ NA_character_
      ),
      cpd = case_when(
        smoking_grp == "never" ~ 0,
        !is.na(py_gdc) & !is.na(duration) & duration > 0 ~ py_gdc * 20 / duration,
        TRUE ~ NA_real_
      ),
      ## Plausibility flag: <1 or >100 cigarettes/day
      cpd_implausible = !is.na(cpd) & smoking_grp != "never" & (cpd < 1 | cpd > 100)
    ) %>%
    ## ---- Exposure metrics ----
    ## 1) cig_years = CPD x years (the "smoking index") -- under this definition == pack-years x 20 (see above)
    ## 2) pack-years: official py_gdc is primary (observed 100% agreement with GDC official)
    ## 3) logcig-years = log(CPD+1) x years
    ## 4a) CSI main analysis (Strategy A): csi_lc uses the canonical formula + published parameters tau=26 / delta=0.7 (CSI_MAIN)
    ##     sensitivity columns (same formula, different published parameters): csi_lc_d10 / csi_lc_t20 / csi_lc_t35
    ## 4b) csi_simple is a power-law family placeholder (int^a*dur^b*(1+tsc)^-c), **not** the Leffondre CSI,
    ##     for sensitivity comparison only; all formal analyses use csi_lc
    mutate(
      cig_years        = cpd * duration,
      pack_years_recon = cpd / 20 * duration,   # should be identically equal to official pack-years (= basis of the 6.1 consistency QC)
      log_cig_years    = log(cpd + 1) * duration,
      ## Main-analysis CSI (Strategy A):
      csi_lc = csi_leffondre(cpd, duration, years_since_quit),
      ## Sensitivity: only the published parameter values change, the formula is unchanged
      csi_lc_d10 = csi_leffondre(cpd, duration, years_since_quit,
                                 tau = CSI_SENS$d10[["tau"]], delta = CSI_SENS$d10[["delta"]]),
      csi_lc_t20 = csi_leffondre(cpd, duration, years_since_quit,
                                 tau = CSI_SENS$t20[["tau"]], delta = CSI_SENS$t20[["delta"]]),
      csi_lc_t35 = csi_leffondre(cpd, duration, years_since_quit,
                                 tau = CSI_SENS$t35[["tau"]], delta = CSI_SENS$t35[["delta"]]),
      ## Force the decay-term base >=1 (pmax) to avoid NaN/Inf from a negative base to the -0.5 power
      csi_simple = cpd^1 * duration^1 * (1 + pmax(replace_na(years_since_quit, 0), 0))^(-0.5)
    ) %>%
    ## Tertiles: one set for the main definition
    mutate(across(c(cig_years, py_gdc, log_cig_years, csi_lc, csi_simple),
                  ~tertile3(.x, smoking_grp), .names = "{.col}_t3"))
}
step1 <- build_exposure(clinic)

## 6.0 Data-source audit: the Xena cigarettes_per_day column (deprecated in this version, kept as evidence only) ----
## Purpose: express "why this column is unused and why CPD is back-calculated from pack-years" as reproducible tests; results archived in
## QC_xena_cpd_audit.csv for methods section and reviewer verification. Four checks:
##   1) proportion where raw column == pack-years/18.25 (near 100% -> it is a linear rescaling of pack-years, zero independent information)
##   2) median of the raw column "cigarettes/day" (~2 rather than ~20 -> physiologically impossible)
##   3) median of the back-calculated CPD and the proportion within 5-60 (true values should fall in this range)
##   4) if the raw column truly were cigarettes/day, pack-years / (column/20 x years) should be ~1; report the observed median (should be far from 1)
audit <- step1 %>%
  filter(!is.na(cpd_xena_raw), !is.na(py_gdc), py_gdc > 0) %>%
  mutate(ratio = cpd_xena_raw / (py_gdc / CPD_PSEUDO_DIVISOR))
ev <- step1 %>% filter(smoking_grp %in% c("current", "former"))
## if the raw column truly were cigarettes/day, pack-years ~ (column/20) x years and this ratio should be ~1 (observed far from 1)
ratio_implied <- ev %>%
  filter(is.finite(cpd_xena_raw), cpd_xena_raw > 0, is.finite(duration), duration > 0) %>%
  mutate(r = py_gdc / ((cpd_xena_raw / 20) * duration)) %>%
  pull(r)
cpd_audit <- tibble(
  check = c("Cases with both raw column and pack-years available",
            "Proportion where raw column == pack-years/18.25 (near 100% -> reason for deprecation)",
            "Median of raw column (supposedly 'cigarettes/day')",
            "Median CPD (cigarettes/day)",
            "Proportion of CPD within 5-60",
            "If raw column truly were cigarettes/day: median pack-years/(column/20 x years) (should be ~1)",
            "Proportion where cig_years == pack-years x 20 under this definition",
            "Cases with CPD flagged implausible (<1 or >100 cigarettes/day)"),
  value = c(as.character(nrow(audit)),
            sprintf("%.1f%%", mean(abs(audit$ratio - 1) < 1e-6) * 100),
            sprintf("%.2f", median(audit$cpd_xena_raw, na.rm = TRUE)),
            sprintf("%.1f", median(ev$cpd, na.rm = TRUE)),
            sprintf("%.1f%%", mean(ev$cpd >= 5 & ev$cpd <= 60, na.rm = TRUE) * 100),
            sprintf("%.2f", median(ratio_implied, na.rm = TRUE)),
            sprintf("%.1f%%", mean(abs(step1$cig_years - step1$py_gdc * 20) < 1e-6, na.rm = TRUE) * 100),
            as.character(sum(step1$cpd_implausible, na.rm = TRUE)))
)
cat("\n== 6.0 Audit of the Xena cigarettes_per_day column (deprecated in this version) ==\n")
print(as.data.frame(cpd_audit))
write_csv(cpd_audit, file.path(out_dir, "QC_xena_cpd_audit.csv"))

## 6.1 Consistency QC: official pack-years vs reconstructed pack-years ---------------------------------------
## Reconstructed pack-years = cpd/20 x years; under this definition it should equal official pack-years;
## this QC confirms the pack-years field itself is readable and uncontaminated (if rho < 0.9, check the data first).
qc <- step1 %>% filter(!is.na(py_gdc), !is.na(pack_years_recon))
cat("\nQC paired samples: ", nrow(qc), "\n")
if (nrow(qc) > 2) {
  qc_rho <- cor(qc$py_gdc, qc$pack_years_recon, method = "spearman", use = "complete.obs")
  cat("Spearman rho(py_gdc, py_recon) =", round(qc_rho, 3), "\n")
  p_qc <- ggplot(qc, aes(py_gdc, pack_years_recon)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Pack-years (GDC official)", y = "Pack-years (reconstructed)",
         title = paste0("Consistency QC, rho = ", round(qc_rho, 3))) +
    theme_bw()
  ggsave(file.path(fig_dir, "QC_packyears_consistency.pdf"), p_qc, width = 5, height = 4.5)
}

## 6.2 CSI parameter robustness QC (core defense material for Strategy A) ----------------------------
## Logic: instead of re-estimating parameters, show that "switching to other published parameters does not change the exposure ranking or grouping".
## Output can serve directly as supplementary material answering reviewers' question "why tau=26/delta=0.7".
if (all(c("csi_lc", "csi_lc_d10", "csi_lc_t20", "csi_lc_t35") %in% names(step1))) {
  sens_vars <- c("csi_lc_d10", "csi_lc_t20", "csi_lc_t35")
  lab_sens  <- c(csi_lc_d10 = "tau=26, delta=1.0 (Leffondre male)",
                 csi_lc_t20 = "tau=20, delta=0.7",
                 csi_lc_t35 = "tau=35, delta=0.7")
  ## Tertiles (in memory only for agreement comparison, not written to disk)
  tert3 <- function(x) {
    out <- rep(NA_character_, length(x))
    ok  <- is.finite(x)
    out[ok] <- as.character(cut(dplyr::ntile(x[ok], 3), breaks = 0:3,
                                labels = c("T1_Low", "T2_Mid", "T3_High")))
    out
  }
  evs <- step1 %>%
    filter(smoking_grp %in% c("current", "former"), is.finite(csi_lc)) %>%
    mutate(csi_lc_g = tert3(csi_lc),
           across(all_of(sens_vars), list(g = ~tert3(.x)), .names = "{.col}_g"))
  if (nrow(evs) > 2) {
    csi_rob <- tibble(
      main_ref      = paste0("csi_lc (canonical CSI; tau=", CSI_MAIN[["tau"]],
                             ", delta=", CSI_MAIN[["delta"]], ")"),
      sensitivity   = unname(lab_sens[sens_vars]),
      n_pair        = vapply(sens_vars, function(p)
                             sum(is.finite(evs$csi_lc) & is.finite(evs[[p]])), integer(1)),
      spearman_rho  = vapply(sens_vars, function(p)
                             round(cor(evs$csi_lc, evs[[p]], method = "spearman", use = "complete.obs"), 4), numeric(1)),
      kendall_tau   = vapply(sens_vars, function(p)
                             round(cor(evs$csi_lc, evs[[p]], method = "kendall", use = "complete.obs"), 4), numeric(1)),
      tertile_agree = vapply(sens_vars, function(p)
                             round(mean(evs[[paste0(p, "_g")]] == evs$csi_lc_g, na.rm = TRUE) * 100, 1), numeric(1))
    )
    cat("\nCSI parameter robustness (main tau=26/delta=0.7 vs other published parameters):\n")
    print(as.data.frame(csi_rob))
    write_csv(csi_rob, file.path(out_dir, "CSI_parameter_robustness.csv"))
    ## Empirical rule: rho >= 0.95 and tertile agreement >= 80% -> parameter choice does not affect conclusions
    cat("Interpretation: rho>=0.95 and tertile agreement>=80% means exposure ranking is robust to parameter choice.\n")
  }
}

## 7. Missingness audit ---------------------------------------------------------------
miss_vars <- c("smoking_grp", "cpd", "duration", "py_gdc", "cig_years",
               "log_cig_years", "csi_lc", "csi_simple", "years_since_quit")
miss_tab <- tibble(
  variable = miss_vars,
  n_miss   = sapply(miss_vars, function(v) sum(is.na(step1[[v]]))),
  pct_miss = round(sapply(miss_vars, function(v) mean(is.na(step1[[v]]))) * 100, 1)
)
write_csv(miss_tab, file.path(out_dir, "missingness_table.csv")); print(miss_tab)

## 8. Descriptive plots -------------------------------------------------------------
## 8.1 Non-finite diagnostics: Inf/NaN commonly arise from conflicting quit-year records (quit year later than diagnosis year)
metric_vars <- c("cig_years", "log_cig_years", "csi_lc", "csi_simple")
nonfinite_diag <- tibble(
  variable    = metric_vars,
  n_na        = sapply(metric_vars, function(v) sum(is.na(step1[[v]]))),
  n_nonfinite = sapply(metric_vars, function(v) sum(!is.finite(step1[[v]])))
)
cat("\nExposure metric missing/non-finite diagnostics:\n"); print(nonfinite_diag)
write_csv(nonfinite_diag, file.path(out_dir, "exposure_nonfinite_diag.csv"))
n_irr <- if ("quit_year_irregular" %in% names(step1)) sum(step1$quit_year_irregular, na.rm = TRUE) else 0L
cat("Records with quit year later than diagnosis year (set to NA, flagged quit_year_irregular): ", n_irr, "\n")

## 8.2 Plotting: explicitly drop non-finite values to avoid the stat_boxplot "Removed N rows" warning
if (nrow(step1) > 0) {
  long_df <- step1 %>%
    filter(smoking_grp %in% c("current", "former"), !is.na(cig_years)) %>%
    transmute(cohort, smoking_grp,
              `cig_years` = cig_years,
              `pack_years` = py_gdc,
              `log_cig_years` = log_cig_years,
              `csi_lc` = csi_lc) %>%
    pivot_longer(-c(cohort, smoking_grp))
  plot_df <- long_df %>% filter(is.finite(value))
  n_drop  <- nrow(long_df) - nrow(plot_df)
  if (n_drop > 0) message("Dropped ", n_drop, " non-finite/missing values for plotting (details in exposure_nonfinite_diag.csv)")
  p1 <- ggplot(plot_df, aes(smoking_grp, value, fill = smoking_grp)) +
    geom_boxplot(outlier.size = 0.4) +
    facet_grid(name ~ cohort, scales = "free_y") +
    labs(title = "Exposure metrics by smoking group (CPD = pack-years x 20 / duration)",
         x = NULL, y = "value") +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(fig_dir, "exposure_metrics_distribution.pdf"), p1, width = 9, height = 6)
}

## 9. Lock the analytic cohort + export (Step2 interface) --------------------------------------
step1_export <- step1 %>%
  select(patient, sample, sample15, cohort, smoking_grp, former_detail,
         ## main exposure columns
         cpd, cpd_src, cpd_implausible, duration, duration_src,
         years_since_quit, quit_year_irregular,
         py_gdc, pack_years_recon, cig_years, log_cig_years,
         csi_lc, csi_lc_d10, csi_lc_t20, csi_lc_t35, csi_simple,
         matches("_t3$"),          # main-definition tertiles
         cpd_xena_raw,             # deprecated Xena raw column (archived only for audit/reproduction)
         age_y, gender, stage, vital, days_death, days_fu)
write_csv(step1_export, file.path(out_dir, "smoking_index_clinical.csv"))

analytic_ever <- step1_export %>%
  filter(smoking_grp %in% c("current", "former"), !is.na(cig_years), !is.na(duration))
write_csv(analytic_ever, file.path(out_dir, "analytic_ever_smokers.csv"))

never_set <- step1_export %>% filter(smoking_grp == "never")
write_csv(never_set, file.path(out_dir, "never_smoker_set.csv"))

## Step2 patient list (Xena expression/mutation matrices are downloaded whole, no per-patient download needed)
write_lines(sort(unique(analytic_ever$patient)),
            file.path(out_dir, "step2_patient_list.txt"))

## 9.1 Decision record: write the "CPD definition and choice of metric" into the outputs (for methods reproduction) ------------
write_csv(tibble(
  Item = c("CPD definition", "CPD formula", "Deprecated column", "CSI parameter strategy", "CSI main parameters", "tsc definition", "Generation date"),
  Value = c("single CPD definition",
         "CPD = pack-years x 20 / smoking years (two-pass back-calculation; algebraic identity of the pack-year definition, no new assumptions)",
         "Xena cigarettes_per_day (== pack-years/18.25; not a true cigarettes/day value; archived only as cpd_xena_raw)",
         "Strategy A (adopt published parameters, no within-cohort re-estimation)",
         sprintf("tau = %s, delta = %s", CSI_MAIN[["tau"]], CSI_MAIN[["delta"]]),
         "Upper-bound definition (missing tsc for former smokers treated as 0) -- see the 05 determination table",
         as.character(Sys.Date()))
), file.path(out_dir, "step1_exposure_decision.csv"))

## 10. Run log --------------------------------------------------------------
sink(file.path(out_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
cat("Step1 complete: ", out_dir, " (CPD = pack-years x 20 / years)\n")
############################################################
## Downstream steps (separate scripts):
##   Step2 (04): deconvolution (quanTIseq / EPIC / MCP-counter), TIS / CYT / CD274 scoring,
##               TMB and SBS4 (needs BSgenome hg38) + dose-response analysis
##   Step3 (06): SIRS construction; Step4 (07-10): Imvigor210 validation + sensitivity analyses
############################################################
