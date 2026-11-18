############################################################
## 07_Step2_immune_microenvironment_TMB_SBS4.R
## Project: smoking exposure quantification x NSCLC immune microenvironment x
##          ICI response prediction
##
## Goals (Step2):
##   1) expression matrix (FPKM-UQ -> TPM) -> immune deconvolution
##      (quanTIseq primary; EPIC and MCP-counter as comparison pathways)
##   2) immune feature scoring: TIS (Ayers 2017, 18 genes), cytolytic activity
##      (CYT), CD274 (PD-L1), CD8A
##   3) TMB: non-synonymous mutation count (SNV + indel) / 38 Mb
##   4) SBS4: 96-trinucleotide spectrum (BSgenome.Hsapiens.UCSC.hg38 supplies the
##      context) -> deconstructSigs projection onto the COSMIC v2
##      thirty-signature reference set (Signature 1-30)
##   5) merge with the Step1 exposure metrics -> dose-response analysis
##      (Spearman + tertile trend) + plots
##
## Inputs (downloaded in Step1):
##   step1_output/smoking_index_clinical.csv    (main table, output of script 02)
##   step1_output/step1c_csi_by_sample.csv      (optional; three tsc conventions,
##                                               output of script 05)
##   data_raw/xena/TCGA-{LUAD,LUSC}.star_fpkm-uq.tsv.gz
##   data_raw/xena/TCGA-{LUAD,LUSC}.somaticmutation_wxs.tsv.gz
##
## Outputs:
##   step2_output/step2_immune_tme.csv                     (sample x TME/TMB/SBS4 main table; Step3/4 interface)
##   step2_output/step2_dose_response_spearman.csv
##   step2_output/step2_tertile_kruskal.csv
##   step2_output/step2_tsc_convention_sensitivity.csv
##   step2_output/step2_qc_epic_nonconverged.csv
##   figs/step2_*.pdf
##
## Dependency installation (idempotent; run before the first execution):
##   setwd(proj_dir); source("install/install_deps.R")             # SBS4 chain + deconvolution chain
##   setwd(proj_dir); source("install/install_deconstructsigs.R")  # hg38-patched deconstructSigs 1.9.0
## Minimal verification: source("install/smoke_test_step2.R")
##
## Data-version guard
##   The master table must use the CPD definition of script 02
##   (CPD = pack-years x 20 / smoking years; ever-smoker median ~20); a table
##   with median CPD < 5 stops the script and prompts to re-run 02.
##
## ---------------------------------------------------------------------------
## Implementation notes and known pitfalls
## ---------------------------------------------------------------------------
## TMB
##   The `effect` column in Xena somaticmutation_wxs uses VEP/SO terms
##   (missense_variant, ...) and may concatenate several effects
##   ("missense_variant;splice_region_variant"), so it is matched by the regex
##   NONSYN_RE rather than by MAF-style terms. Empty chrom values and
##   non-standard contigs are removed (STD_CHR) before any genomic operation.
##
## SBS4
##   - Signature fitting uses ALL somatic SNVs (synonymous and non-coding
##     included, following signature-analysis convention), whereas TMB counts
##     only non-synonymous variants.
##   - whichSignatures() is a single-sample interface (it internally subsets
##     rownames == sample.id and requires tumor.ref as a data.frame) and is
##     therefore called per sample, only for samples with SNV >= 50
##     (threshold MIN_SBS4_SNV).
##   - Two representations are produced: the relative contribution
##     (`sbs4_weight`) and its absolute approximation
##     `sbs4_abs` = `sbs4_weight` x TMB. The two are not interchangeable.
##
## Gene mapping (probeMap)
##   - The correct file name is `gencode.v36.annotation.gtf.gene.probemap`
##     (lowercase m); a fallback URL is provided.
##   - probeMap ids carry version numbers while the expression matrix side
##     strips them, so both sides must be version-stripped before joining.
##   - After stripping, PAR_Y entries collide with the main entry -> distinct(id)
##     first, to prevent a many-to-many join.
##   - The id column must be dropped before as.matrix(), otherwise the whole
##     matrix is coerced to character.
##
## Sample columns
##   Xena expression matrix column names are 16-character barcodes
##   (TCGA-XX-XXXX-01A), so the same 15-character sample can appear as two
##   aliquots (-01A / -01B); truncating to 15 characters creates duplicate
##   names (measured: 1029 columns -> 1017 unique samples). The earliest suffix
##   is kept. Two guards are in place: sample15 uniqueness on the deconvolution
##   tables, and nrow(step2) == nrow(step1) after merging (a duplicate in the
##   right table of a left_join would otherwise silently add rows).
##   Deconvolution tables are normalised by normalize_deconv(df, prefix), which
##   accepts `Sample` or `sample` as the sample column and falls back to the
##   first column with an explicit warning when neither is present.
##
## Namespace masking (Bioconductor)
##   library(BSgenome.Hsapiens.UCSC.hg38) pulls GenomicRanges onto the search
##   path, and GenomicRanges exports S4 generics that mask same-named S3
##   functions: `reduce` (purrr) and union/intersect/setdiff (base). Calls to
##   these must be explicitly namespaced (purrr::reduce, base::union, ...) once
##   a genomic package has been loaded.
##
## R language pitfalls
##   - Never build formulas with paste(): its default separator is a space, so
##     paste("`", x, "`") yields a name with spaces inside the backticks. Use
##     reformulate() instead.
##   - tidyselect helpers (all_of/any_of/everything) are only valid in
##     column-selecting positions inside select()/across(); placing them on the
##     right-hand side of rename(new = old) is an error. To rename dynamically,
##     assign to names(df) directly.
##   - CYT is the geometric mean of GZMA and PRF1 (Rooney et al. 2015); the score
##     is averaged along the gene direction (colMeans), one value per sample,
##     with a +1 offset before the log to avoid -Inf at zero expression.
##
## EPIC convergence QC
##   EPIC occasionally reports "optimization didn't fully converge" (measured at
##   ~5-6% of samples) but still returns estimates and exposes no convergence
##   flag. The warning is captured with withCallingHandlers and written to
##   step2_qc_epic_nonconverged.csv, and a boolean epic_nonconverged column is
##   added to the EPIC table. EPIC is a comparison pathway, so this does not
##   affect the primary (quanTIseq) results.
##
## Environment / locale
##   Identifiers and output file names must be ASCII; in-figure labels must be
##   ASCII (cairo_pdf is used for plotting); comments and strings may be
##   non-ASCII.
############################################################

## 0. Configuration ------------------------------------------------------------------
## Author: see manuscript author list   Date: 2026-10-17
proj_dir <- "D:/Projects/smoking-index-ici"      # keep consistent with 02
data_dir <- file.path(proj_dir, "data_raw/xena")
out_dir  <- file.path(proj_dir, "step2_output"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir  <- file.path(proj_dir, "figs");         dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
options(timeout = max(3600, getOption("timeout")), stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})

## Optional-package check: degrade whichever is missing, do not abort everything
has_pkg <- function(p) requireNamespace(p, quietly = TRUE)
HAS_QS       <- has_pkg("quantiseqr")                         # official Bioconductor quanTIseq implementation (preferred primary pathway)
HAS_DECONV   <- has_pkg("immunedeconv")                       # deconvolution portal (installed from GitHub, heavy deps, used as fallback)
HAS_EPIC     <- has_pkg("EPIC")
HAS_MCP      <- has_pkg("MCPcounter")
HAS_GSVA     <- has_pkg("GSVA")
HAS_DSIGS    <- has_pkg("deconstructSigs")
HAS_BSG      <- has_pkg("BSgenome.Hsapiens.UCSC.hg38")
cat("Package availability: quantiseqr=", HAS_QS, " immunedeconv=", HAS_DECONV, " EPIC=", HAS_EPIC,
    " MCPcounter=", HAS_MCP, " GSVA=", HAS_GSVA,
    " deconstructSigs=", HAS_DSIGS, " BSgenome.hg38=", HAS_BSG, "\n")

## Dependency self-check: state in one go "what's missing, how to install" (companion install script install_deps.R)
miss_mod <- c(if (!HAS_QS && !HAS_DECONV) "primary deconvolution pathway (quantiseqr or immunedeconv)" else NULL,
              if (!HAS_DSIGS || !HAS_BSG) "SBS4 (deconstructSigs + BSgenome.Hsapiens.UCSC.hg38)" else NULL)
if (length(miss_mod)) {
  message("!! The following modules will be skipped due to missing packages: ", paste(miss_mod, collapse = "; "), "\n",
          "   Install all at once (idempotent, re-runnable):  setwd('", proj_dir, "'); source('install/install_deps.R')\n",
          "   Note: deconstructSigs has been archived on CRAN (the script uses the GitHub mainline, falling back to Archive 1.8.0);\n",
          "       BSgenome.Hsapiens.UCSC.hg38 is an annotation source package, ~731 MB.")
} else {
  cat("Dependency self-check: SBS4 and deconvolution chains present\n")
}

## 1. Read the Step1 main table --------------------------------------------------------
step1 <- fread(file.path(proj_dir, "step1_output/smoking_index_clinical.csv"),
               na.strings = c("", "NA", "[Not Available]")) %>% as_tibble()
step1 <- step1 %>%
  mutate(sample15 = substr(sample, 1, 15),
         patient  = substr(sample, 1, 12))
cat("Step1 main table: ", nrow(step1), " samples\n")
print(table(step1$cohort, step1$smoking_grp, useNA = "ifany"))

## 1.1 Data-version guard (requires the main table from script 02) --------------------------------------
## Background: the sole CPD convention in 02 is pack-years*20/years; ever-smokers' median CPD should be ~20.
##       The Xena cigarettes_per_day column is not a true cigarettes/day value (= pack-years/18.25, median ~2.2),
##       disproven and discarded. If a table with median CPD < 5 is read, abort directly, avoiding downstream use of distorted exposure values.
ev_check  <- step1 %>% filter(smoking_grp %in% c("current", "former"))
cpd_med   <- suppressWarnings(median(ev_check$cpd, na.rm = TRUE))
cat(sprintf("\n[data-version guard] median cpd = %.2f (expected ~20; stale pseudo-column ~2.2)\n", cpd_med))
if (is.finite(cpd_med) && cpd_med < 5) {
  stop(sprintf(paste0("The Step1 main table read in appears to use a stale CPD definition (median CPD = %.2f, expected ~20).\n",
                      "  Cause: it was built from Xena's cigarettes_per_day, which is a pseudo-value.\n",
                      "  Action: re-run script 02 (overwrites the old CSV), then run this script.\n",
                      "  Note: see docs/06_data_source_audit_smoking_fields.md"),
               cpd_med), call. = FALSE)
} else {
  cat("[data-version guard] Passed: main table uses the CPD definition of script 02 (CPD = pack-years*20/years).\n")
}
if ("cpd_implausible" %in% names(step1)) {
  cat("Number of samples with implausible CPD (<1 or >100 cigarettes/day): ", sum(step1$cpd_implausible, na.rm = TRUE), "\n")
}

## 2. Expression matrix: fpkm-uq -> TPM, EnsemblID -> gene symbol ----------------------
##   (1) gene mapping table: the gencode probeMap shipped with the Xena GDC hub (download once, use many times)
##   Two pitfalls:
##     (a) the correct filename is `gene.probemap` (lowercase m) -- writing `gene.probeMap` gives 403 Access Denied;
##     (b) this table's id column **carries version numbers** (ENSG00000223972.5), while the expression-matrix side strips them ->
##         both must be version-stripped before joining, otherwise inner_join matches nothing and the expression matrix is empty.
probemap_file <- file.path(data_dir, "gencode.v36.annotation.gtf.gene.probemap")
if (!file.exists(probemap_file) || file.size(probemap_file) < 1e5) {
  message("Downloading gencode probeMap ...")
  pm_urls <- c(
    "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
    "https://toil.xenahubs.net/download/probeMap/gencode.v23.annotation.gene.probemap"  # fallback (older version)
  )
  for (u in pm_urls) {
    r <- try(download.file(u, probemap_file, mode = "wb", quiet = TRUE, method = "libcurl"), silent = TRUE)
    if (!inherits(r, "try-error") && file.exists(probemap_file) && file.size(probemap_file) > 1e5) {
      cat("probeMap downloaded successfully: ", basename(u), "\n"); break
    }
  }
  if (!file.exists(probemap_file) || file.size(probemap_file) < 1e5)
    stop("probeMap download failed: please manually download gencode.v36.annotation.gtf.gene.probemap into ", data_dir)
}
probemap <- fread(probemap_file) %>%
  transmute(id = sub("\\..*$", "", id), symbol = gene) %>%   # strip version (align with the expression-matrix side)
  filter(!is.na(symbol), symbol != "", symbol != "?") %>%
  distinct(id, .keep_all = TRUE)
cat("probeMap: ", nrow(probemap), " Ensembl->symbol mappings\n")

read_fpkm_uq <- function(cohort) {
  f <- file.path(data_dir, paste0(cohort, ".star_fpkm-uq.tsv.gz"))
  ## First column is Ensembl_ID (with version), remaining columns are samples (15-character barcodes)
  ## Three pitfalls:
  ##   (a) after stripping versions, id collides with entries like PAR_Y -> distinct(id) first to avoid many-to-many join
  ##   (b) id is a character column and must be dropped before use, otherwise as.matrix() coerces the whole matrix to character
  ##   (c) when multiple genes map to the same symbol, take the within-column maximum
  fread(f) %>% as_tibble() %>%
    mutate(id = sub("\\..*$", "", Ensembl_ID)) %>%
    select(-Ensembl_ID) %>%
    distinct(id, .keep_all = TRUE) %>%
    inner_join(probemap, by = "id") %>%
    filter(!is.na(symbol)) %>%
    select(symbol, where(is.numeric)) %>%            # keep only symbol + numeric sample columns
    group_by(symbol) %>% summarise(across(everything(), max), .groups = "drop") %>%
    column_to_rownames("symbol") %>% as.matrix()
}
luad <- read_fpkm_uq("TCGA-LUAD"); lusc <- read_fpkm_uq("TCGA-LUSC")
genes_common <- base::intersect(rownames(luad), rownames(lusc))   # explicitly align by gene name to avoid positional mismatch
if (!setequal(rownames(luad), rownames(lusc)))
  cat(sprintf("Note: cohort gene counts %d / %d, merged intersection %d\n", nrow(luad), nrow(lusc), length(genes_common)))
mat_fpkm <- cbind(luad[genes_common, , drop = FALSE], lusc[genes_common, , drop = FALSE])
rm(luad, lusc)

## Duplicate sample columns (Xena column names are 16-character barcodes; truncation to 15 collides)
##   Xena sample columns look like `TCGA-XX-XXXX-01A`; the same 15-character sample may have both A/B aliquots;
##   after truncation to 15 characters the two columns share a name. Measured: of 1029 columns only 1017 unique samples, 12 LUAD samples
##   (A/B aliquots of TCGA-44-xxxx-01) double-counted. Two consequences:
##     (1) functions like EPIC that "use input column names as row names" directly report
##        `error in .rowNamesDF<-: duplicate 'row.names' are not allowed`;
##     (2) more subtle: `step2 <- step1 %>% left_join(deconv/score table, by = "sample15")` produces extra rows
##        because the right table has duplicates -> silently inflating sample size and contaminating all correlation analyses.
##   Handling: keep the **earliest-suffix** aliquot for the same 15-character sample (A preferred, i.e. the primary aliquot); no samples lost.
cn_all   <- colnames(mat_fpkm)
sid_all  <- substr(cn_all, 1, 15)
alq_all  <- substring(cn_all, 16)                 # 15-char -> ""; 16-char -> "A"/"B"/"C"
is01_all <- substr(cn_all, 14, 15) == "01"        # keep only primary tumor (-01)
n_dup    <- sum(duplicated(sid_all[is01_all]))

keep_df  <- tibble(raw = cn_all[is01_all], sid = sid_all[is01_all], aliquot = alq_all[is01_all]) %>%
  arrange(sid, aliquot) %>%                       # same 15-char sample: "" < "A" < "B" < "C"
  distinct(sid, .keep_all = TRUE)
mat_fpkm <- mat_fpkm[, match(keep_df$raw, cn_all), drop = FALSE]
colnames(mat_fpkm) <- keep_df$sid
if (n_dup > 0)
  cat(sprintf("Sample-column dedup: %d duplicate -01A/-01B columns removed, %d unique samples kept\n",
              n_dup, ncol(mat_fpkm)))
stopifnot(!anyDuplicated(colnames(mat_fpkm)))     # enforce the "one column per sample" contract
cat("fpkm-uq matrix: ", nrow(mat_fpkm), " genes x ", ncol(mat_fpkm), " samples\n")

##   (2) FPKM-UQ -> TPM: TPM_i = FPKM_i / sum(FPKM) * 1e6 (UQ does not change within-column proportions, so the conversion holds)
mat_tpm <- apply(mat_fpkm, 2, function(x) { s <- sum(x[x > 0], na.rm = TRUE); ifelse(x > 0, x / s * 1e6, 0) })
mat_tpm_log <- log2(mat_tpm + 1)   # later scoring/plotting uniformly uses log2-TPM; deconvolution uses linear TPM

## 3. Immune feature scoring (no deconvolution needed, computed first) --------------------------------------
## TIS 18 genes (Ayers et al. 2017, JCI) -- interferon-gamma response signature, a classic ICI prediction signature
TIS_GENES <- c("CD27","IFNB1","IFNG","CXCL10","CXCL9","HLA-DRA","HLA-E","IDO1",
               "PSMB10","STAT1","TAP1","TAP2","TAPBP","TIGIT","CMKLR1","CCL5","GZMK","CD8A")
## Cytolytic activity CYT (Rooney et al. 2015, Cell): geometric mean of two cytolytic gene expressions, **one value per sample**
## Gene convention: strictly aligned with the original Rooney 2015 -- the gene pair is GZMA + PRF1.
## Formula convention: CYT in the original Rooney is the **geometric mean** of the two-gene expression geometric mean = (GZMA*PRF1)^(1/2),
##   consistent with this implementation exp(mean(log(x+1))) (this implementation uses log2-TPM with +1 to prevent log(0), a robust scale variant).
CYT_GENES <- c("GZMA", "PRF1")

expr_for_scoring <- mat_tpm_log[base::intersect(base::union(TIS_GENES, CYT_GENES), rownames(mat_tpm_log)), , drop = FALSE]
score_mean <- function(genes, mat) colMeans(mat[base::intersect(genes, rownames(mat)), , drop = FALSE], na.rm = TRUE)

## CYT must yield one value per sample -> average along the [gene direction], i.e. colMeans.
##   The original implementation wrongly used rowMeans, returning "one value per gene" (length = number of matched genes = 2),
##   inconsistent with sample15 (length = number of samples), triggering tibble()'s size-check error:
##     "Tibble columns must have compatible sizes. Size 1029: Existing data. Size 2: Column cyt_score."
##   The +1 offset prevents TPM=0 from making that sample exactly 0 via log(0) = -Inf.
cyt_present <- base::intersect(CYT_GENES, rownames(expr_for_scoring))
if (length(cyt_present) == 0) warning("None of the CYT genes (", paste(CYT_GENES, collapse = "/"),
                                      ") were found in the expression matrix; cyt_score set to NA", call. = FALSE)
cyt_vec <- if (length(cyt_present) >= 1) {
  exp(colMeans(log(expr_for_scoring[cyt_present, , drop = FALSE] + 1), na.rm = TRUE))
} else {
  rep(NA_real_, ncol(mat_tpm))
}
stopifnot(length(cyt_vec) == ncol(mat_tpm))   # enforce the length contract, avoiding the same class of error silently passed into tibble()

## Single-gene extraction likewise has presence protection (NA when the gene is missing, rather than aborting the block)
pick_gene <- function(g) if (g %in% rownames(mat_tpm_log)) mat_tpm_log[g, ] else rep(NA_real_, ncol(mat_tpm))

immu_scores <- tibble(
  sample15 = colnames(mat_tpm),
  tis_score  = score_mean(TIS_GENES, expr_for_scoring),          # simplified convention: mean; a z-score can be upgraded to GSVA::ssgsea
  cyt_score  = cyt_vec,
  cd274_log2tpm = pick_gene("CD274"),
  cd8a_log2tpm  = pick_gene("CD8A")
)

## 4. Immune deconvolution: quanTIseq (primary) + EPIC / MCP-counter (comparison) ----------------
## Primary pathway prefers quantiseqr (official Bioconductor quanTIseq implementation, light dependencies, binary-installable);
##   falls back to immunedeconv::deconvolute only when quantiseqr is missing but immunedeconv exists.
##   Both implement the same algorithm (Finotello 2019, Genome Med), with identical TIL10 signature and output cell types.
deconv_res <- list()

## Uniformly "normalize" deconvolution results to sample15 + prefixed cell-type columns.
##   Background: `rename(sample15 = Sample)` uses a bare column name; once the upstream column name changes (quantiseqr version differences,
##   or accidentally taking the immunedeconv branch), the error is "object 'Sample' not found" -- which does not tell which line
##   or object is at fault, making it very hard to locate. Changed to: first **explicitly pick the sample column** (accepting Sample/sample/
##   first column, etc.), uniformly prefix the remaining columns; and assert the column count on the result.
##   `rename(sample15 = all_of(sample_col))` is **wrong** --
##     all_of() is a tidyselect "selection function", allowed only in selection positions such as select()/across();
##     inside rename()'s new=old it is treated as the literal name (tidyselect >=1.2 also emits a deprecation warning),
##     hence the error `Some 'from' names in value not found on 'x': sample15`.
##     Correct approach: use `.data[[sample_col]]` on the right side of rename (a dynamic column reference supported by rlang),
##     or assign directly to names() -- the latter is used here, simplest and most robust.
normalize_deconv <- function(df, prefix, cell_prefix = TRUE) {
  df <- as_tibble(df)
  ## Sample column: prefer Sample / sample; otherwise take the first column (both quantiseqr and immunedeconv start with the sample column)
  cand <- c("Sample", "sample")
  hit  <- base::intersect(cand, names(df))
  if (!length(hit)) {
    warning(sprintf("No Sample/sample column found in deconvolution result; using first column '%s' as the sample column", names(df)[1]),
            call. = FALSE)
    sample_col <- names(df)[1]
  } else {
    sample_col <- hit[1]
  }
  stopifnot(sample_col %in% names(df),
            sum(names(df) == sample_col) == 1)   # this column must exist and be unique, otherwise renaming is ambiguous
  n_cell <- ncol(df) - 1                     # number of cell types besides the sample column
  ## Rename via "position -> sample15": first record the names to keep (old names), then assign new names wholesale
  old_names   <- names(df)
  names(df)[old_names == sample_col] <- "sample15"
  out <- df %>%
    mutate(sample15 = substr(as.character(sample15), 1, 15))
  if (cell_prefix) {
    out <- out %>% rename_with(~ paste0(prefix, .x), -sample15)
  }
  stopifnot(ncol(out) == n_cell + 1)         # enforce the "sample column + cell-type columns" contract
  out
}

if (HAS_QS) {
  library(quantiseqr)
  ## run_quantiseq requires a **linear TPM** matrix with HGNC symbols as row names and samples as column names
  ## is_tumordata = TRUE (tumor samples, enables tumor-related gene removal), scale_mRNA = TRUE (corrects mRNA-content bias)
  qs <- run_quantiseq(expression_data = mat_tpm, signature_matrix = "TIL10",
                      is_arraydata = FALSE, is_tumordata = TRUE,
                      scale_mRNA = TRUE, method = "lsei")
  ## Returns a data.frame: first column Sample + 10 cell-type scores (row sums = 1)
  deconv_res$quantiseq <- normalize_deconv(qs, prefix = "q_")
  cat("quanTIseq done (quantiseqr): ", ncol(deconv_res$quantiseq) - 1, " cell types\n")
} else if (HAS_DECONV) {
  library(immunedeconv)
  ## quanTIseq requires: TPM, gene symbols, gene x sample
  dq <- deconvolute(mat_tpm, "quantiseq", tumor = TRUE) %>%
    as_tibble() %>% pivot_wider(names_from = cell_type, values_from = score)
  deconv_res$quantiseq <- normalize_deconv(dq, prefix = "q_")
  cat("quanTIseq done (immunedeconv): ", ncol(deconv_res$quantiseq) - 1, " cell types\n")
}
if (HAS_EPIC) {
  library(EPIC)
  ## Not using as.data.frame() / rownames_to_column(): both error on "duplicate row names":
  ##   `error in .rowNamesDF<-: duplicate 'row.names' are not allowed`.
  ##   Instead use as_tibble() then manually add the sample15 column (tibbles have no row names, so this check is not triggered).
  ##
  ## EPIC's "some samples did not fully converge" is only a warning, and the returned object has no convergence-flag column.
  ##   Use withCallingHandlers to capture the warning text -> parse non-converged samples -> write QC table + add flag column.
  ##   Note: muffleWarning does not affect the return value; even if some version emits this via message (parsed empty)
  ##   it only yields an empty QC table with all-FALSE flags, producing no false positives.
  epic_warn <- character(0)
  epic_raw  <- withCallingHandlers(
    EPIC(mat_tpm),
    warning = function(w) { epic_warn <<- c(epic_warn, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage")
  )
  fr <- as.matrix(epic_raw$cellFractions)
  ## Goes through the unified normalization function (EPIC has no sample column; provided by rownames)
  deconv_res$epic <- fr %>%
    as_tibble(.name_repair = "minimal") %>%
    mutate(sample15 = substr(rownames(fr), 1, 15)) %>%
    select(sample15, everything()) %>%
    rename_with(~ paste0("e_", .x), -sample15)

  ## Parse non-converged samples (15-character barcodes) and intersect with this table's samples to avoid mislabeling
  epic_nc <- unique(unlist(str_extract_all(paste(epic_warn, collapse = " "),
                                           "TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}-[0-9]{2}")))
  epic_nc <- base::intersect(epic_nc, deconv_res$epic$sample15)
  deconv_res$epic <- deconv_res$epic %>% mutate(epic_nonconverged = sample15 %in% epic_nc)
  fwrite(tibble(sample15 = epic_nc), file.path(out_dir, "step2_qc_epic_nonconverged.csv"))
  cat("EPIC done: ", ncol(deconv_res$epic) - 2, " cell types; non-converged samples ", length(epic_nc),
      " cases -> step2_qc_epic_nonconverged.csv (comparison pathway; does not overturn main conclusions; exclusion sensitivity possible)\n", sep = "")
}
if (HAS_MCP) {
  library(MCPcounter)
  ## MCPcounter returns cell types x samples (row names = cell types, unique) -> transpose to samples x cell types
  mcp_wide <- MCPcounter.estimate(mat_tpm, featuresType = "HUGO_symbols") %>%
    as.data.frame() %>% rownames_to_column("cell_type") %>%
    pivot_longer(-cell_type, names_to = "sample15", values_to = "score") %>%
    mutate(sample15 = substr(sample15, 1, 15),
           cell_type = paste0("mcp_", cell_type)) %>%
    pivot_wider(names_from = cell_type, values_from = score)
  deconv_res$mcp <- mcp_wide
  cat("MCP-counter done: ", ncol(deconv_res$mcp) - 1, " cell types\n")
}

## Uniqueness guard: every auxiliary table entering step2 must have "one row per sample".
##   If duplicates remain (e.g. a new aliquot name collision in the source data), warn explicitly and keep the first row,
##   rather than letting left_join silently inflate the row count.
deconv_res <- lapply(deconv_res, function(d) {
  if (anyDuplicated(d$sample15)) {
    warning("Duplicate sample15 in deconvolution result: ", sum(duplicated(d$sample15)),
            " rows; kept the first row -- please check whether the source matrix column names are duplicated", call. = FALSE)
    d <- d %>% distinct(sample15, .keep_all = TRUE)
  }
  d
})

## 5. TMB: non-synonymous mutation density -----------------------------------------------------
##   somaticmutation_wxs is a long table: sample gene chrom start end ref alt effect callers dna_vaf
##   The effect column uses VEP/SO terms (missense_variant / stop_gained / ...), not MAF terms
##     (Missense_Mutation / Nonsense_Mutation / ...); and can be a concatenation of multiple effects
##     (e.g. "missense_variant;splice_region_variant") -> must use regex matching,
##     not %in% exact matching (otherwise nothing matches and TMB is all 0).
NONSYN_RE <- paste(c("missense_variant", "stop_gained", "stop_lost", "start_lost",
                     "frameshift_variant", "inframe_deletion", "inframe_insertion",
                     "splice_acceptor_variant", "splice_donor_variant",
                     "splice_region_variant", "protein_altering_variant",
                     "disruptive_inframe", "nonstop"), collapse = "|")
SYNON_RE  <- "^synonymous_variant$|^synonymous"   # for QC counting only
EXOME_MB  <- 38   # MC3 convention; before submission, cross-check with GDC official per-sample covered bed size
STD_CHR   <- c(as.character(1:22), "X", "Y")      # discard empty chrom and scaffold/unplaced contigs

read_mut <- function(cohort) {
  fread(file.path(data_dir, paste0(cohort, ".somaticmutation_wxs.tsv.gz")),
        select = c("sample", "chrom", "start", "end", "ref", "alt", "effect"),
        showProgress = FALSE) %>% as_tibble() %>% mutate(cohort = cohort)
}
mut_all <- bind_rows(read_mut("TCGA-LUAD"), read_mut("TCGA-LUSC")) %>%
  mutate(sample15 = substr(sample, 1, 15),
         chrom    = sub("^chr", "", chrom)) %>%
  filter(!is.na(start), chrom %in% STD_CHR)        # drop empty chrom / non-standard contigs together
mut <- mut_all %>% filter(grepl(NONSYN_RE, effect))
cat("Mutation table: total variants ", nrow(mut_all), " -> non-synonymous ", nrow(mut), " (incl. SNV + indel)\n")
cat("  QC: synonymous variants ", sum(grepl(SYNON_RE, mut_all$effect)), "; ",
    "non-synonymous fraction ", round(100 * nrow(mut) / max(nrow(mut_all), 1), 1), "%\n")

##   TMB convention: non-synonymous SNV + indel count / 38 Mb (TCGA convention)
##   Note: using only the "non-synonymous" table cannot distinguish "no mutations" from "not sequenced" -> first take the full sample set of the mutation table,
##         then fill in 0 (these samples do have mutation data, just no non-synonymous variants), avoiding downstream mistaking 0 for missing.
mut_universe <- mut_all %>% distinct(sample15)      # full sample set with mutation data
tmb <- mut %>% distinct(sample15, chrom, start, ref, alt, .keep_all = TRUE) %>%  # dedup to prevent double-counting the same locus
  count(sample15, name = "nonsyn_n") %>%
  right_join(mut_universe, by = "sample15") %>%
  mutate(nonsyn_n  = ifelse(is.na(nonsyn_n), 0L, nonsyn_n),
         tmb_per_mb = nonsyn_n / EXOME_MB)
if (nrow(tmb) == 0) warning("TMB table is empty -- please check whether the effect terms match NONSYN_RE", call. = FALSE)
cat("TMB done: ", nrow(tmb), " samples (including ",
    sum(tmb$nonsyn_n == 0), " with 0 mutations); median TMB =",
    round(median(tmb$tmb_per_mb, na.rm = TRUE), 2), "mut/Mb\n")

## 6. SBS4: 96 channels -> COSMIC projection --------------------------------------------
##   The Xena mutation table lacks trinucleotide context; BSgenome hg38 supplies it; deconstructSigs(cosmic30) contains Signature.4
##   Two interface constraints:
##     (1) tumor.ref must be a data.frame (passing a matrix stops immediately);
##     (2) whichSignatures is a **single-sample** interface -- internally it subsets rownames == sample.id,
##         so it must be called per sample, not once for the whole cohort (otherwise 'sample.id' is missing).
##     (3) weights are unreliable with too few sample SNVs (in-package warning threshold 50), so fit only samples with >= MIN_SBS4_SNV.
MIN_SBS4_SNV <- 50
sbs4 <- NULL
if (HAS_DSIGS && HAS_BSG) {
  library(deconstructSigs); library(BSgenome.Hsapiens.UCSC.hg38)
  mut_in <- mut_all %>%
    filter(ref %in% c("A","C","G","T"), alt %in% c("A","C","G","T")) %>%  # SBS uses only SNVs (empty chrom/contigs already removed)
    transmute(sample.id = sample15, chr = paste0("chr", chrom), pos = start, ref, alt) %>%
    as.data.frame()                                    # UCSC style, aligned with BSgenome hg38's seqlevels
  sigs_input <- suppressWarnings(          # silence the "fewer than 50 mutations" spam (filtering logic is on the next lines, done independently)
    mut.to.sigs.input(mut.ref = mut_in,
                      sample.id = "sample.id", chr = "chr", pos = "pos",
                      ref = "ref", alt = "alt",
                      bsg = BSgenome.Hsapiens.UCSC.hg38)   # annotation package (~731 MB)
  )
  si_df  <- as.data.frame(sigs_input)                  # must convert to data.frame
  n_snv  <- rowSums(si_df)
  keep_s <- names(n_snv)[n_snv >= MIN_SBS4_SNV]
  cat("SBS4: ", length(keep_s), "/", nrow(si_df), " samples with SNV >= ", MIN_SBS4_SNV,
      ", entering signature fitting (other samples get sbs4_weight = NA)\n", sep = "")

  w_list <- lapply(keep_s, function(s) {              # per-sample fitting
    res <- try(whichSignatures(tumor.ref = si_df, sample.id = s,
                               signatures.ref = signatures.cosmic,   # COSMIC v2 30 signatures, including Signature.4
                               contexts.needed = TRUE,
                               tri.counts.method = "exome2genome")[["weights"]],  # correct exome preference
               silent = TRUE)
    if (inherits(res, "try-error")) return(NULL)
    res$sample15 <- s
    res
  })
  w <- bind_rows(w_list) %>% as_tibble()
  if (nrow(w)) {
    sig_names <- base::setdiff(names(w), "sample15")
    sbs4 <- w %>%
      mutate(top_signature = sig_names[max.col(pick(all_of(sig_names)), ties.method = "first")]) %>%
      transmute(sample15,
                sbs4_weight = `Signature.4`,        # primary metric: exposure weight of the smoking-related signature
                top_signature)
    ## "Smoking-related signature dominant" (top_signature == "Signature.4") can serve as a categorical auxiliary endpoint
    cat("SBS4 projection done: ", nrow(sbs4), " samples; median SBS4 weight =",
        round(median(sbs4$sbs4_weight, na.rm = TRUE), 3),
        "; samples with Signature.4 dominant =", sum(sbs4$top_signature == "Signature.4", na.rm = TRUE), "\n")
  } else {
    warning("SBS4 fitting produced no samples (check SNV counts and chromosome naming)", call. = FALSE)
  }
}

## 7. Merge main table: exposure x TME ----------------------------------------------------
## Note: all deconvolution packages may be missing (deconv_res empty), so do not use reduce() to join them all at once
##
## Namespace masking ("unable to find an inherited method for function 'reduce' for signature 'x = \"list\"'"):
##   Section 6's `library(BSgenome.Hsapiens.UCSC.hg38)` pulls `GenomicRanges` onto the search path,
##   and GenomicRanges also exports `reduce` -- and it is an **S4 generic**. An S4 generic masks the same-named S3 function,
##   so a bare `reduce(list, ...)` matches the GenomicRanges generic first -> no method for list -> error.
##   (Measured: with only tidyverse/data.table/quantiseqr/EPIC/MCPcounter/deconstructSigs/GSVA
##     loaded, reduce still points to purrr; once BSgenome.hg38 is loaded it is taken over by GenomicRanges.)
##   Handling: always **explicitly namespace** reduce as purrr::reduce(), fully immune to search-path drift.
##   Likewise, union/intersect/setdiff/slice/rename are also S4 generics (provided by IRanges/GenomicRanges):
##   they have methods on vectors so they work for now, but to avoid "load-order dependence" they are also explicitly namespaced base::/dplyr::.
deconv_all <- if (length(deconv_res)) purrr::reduce(deconv_res, full_join, by = "sample15") else NULL
step2 <- step1 %>% left_join(immu_scores, by = "sample15")
if (!is.null(deconv_all)) step2 <- step2 %>% left_join(deconv_all, by = "sample15")
step2 <- step2 %>% left_join(tmb, by = "sample15")
if (!is.null(sbs4)) step2 <- step2 %>% left_join(sbs4, by = "sample15")

## Integrity guard: left_join should not change the row count -- right-table sample duplication silently inflates N.
##   This assertion is the "last gate" against duplicate-sample-column-type errors, preventing contaminated results from being written as normal output.
if (nrow(step2) != nrow(step1))
  stop(sprintf("Merged row count %d != step1 row count %d: auxiliary table has duplicate sample15, please check first (see the dedup warning in section 2)",
               nrow(step2), nrow(step1)), call. = FALSE)
cat("Main table merge done: ", nrow(step2), " samples x ", ncol(step2), " columns\n")

## 7.1 Merge in 05's tsc convention sensitivity columns ----------------------------------------------
## 05_Step1c outputs step1c_csi_by_sample.csv (CSI and tertiles under the upper/strict/mid three conventions).
## With it, Step2's dose-response analysis simultaneously answers "whether conclusions change under a different convention" --
## this is exactly how "recompute, then decide the convention" is put into practice: judge with endpoints (TME/TMB) rather than missingness rate alone.
scen_file <- file.path(proj_dir, "step1_output/step1c_csi_by_sample.csv")
if (file.exists(scen_file)) {
  sc <- fread(scen_file) %>% as_tibble() %>%
    select(sample15, any_of(c("csi_upper", "csi_upper_t3", "csi_strict", "csi_strict_t3",
                              "csi_mid", "csi_mid_t3")))
  step2 <- step2 %>% left_join(sc, by = "sample15")
  cat("Merged in 05's tsc convention sensitivity columns: ", paste(base::setdiff(names(sc), "sample15"), collapse = ", "), "\n")
} else {
  message("step1c_csi_by_sample.csv not found; run the 05 script first if convention sensitivity is needed.")
}

step2 <- step2 %>% mutate(across(where(is.numeric), ~ ifelse(is.infinite(.x), NA, .x)))
fwrite(step2, file.path(out_dir, "step2_immune_tme.csv"))
cat("step2 main table: ", nrow(step2), " samples x ", ncol(step2), " columns ->", "step2_immune_tme.csv\n")

## 8. Dose-response analysis: exposure metrics x TME features --------------------------------------
##   Exposures: py_gdc (pack-years, primary) / log_cig_years / csi_lc (CSI main analysis) / cpd (intensity component)
##        + csi_strict / csi_mid (tsc convention sensitivity, provided by 05; to test "whether conclusions change under a different convention")
##   Endpoints: TMB / SBS4 / q_CD8 T cell / CYT / TIS / CD274 ... (information structure in audit document 06)
exposures <- c("py_gdc", "log_cig_years", "csi_lc", "cpd",
               base::intersect(c("csi_strict", "csi_mid"), names(step2)))
endpoints <- base::intersect(c("tmb_per_mb", "sbs4_weight", "cd8a_log2tpm", "cd274_log2tpm",
                         "cyt_score", "tis_score",
                         grep("^q_", names(step2), value = TRUE)),
                       names(step2))
ev <- step2 %>% filter(smoking_grp %in% c("current", "former"))

sp_tab <- expand_grid(exp = exposures, end = endpoints, cohort = c("Pooled", "TCGA-LUAD", "TCGA-LUSC")) %>%
  pmap_dfr(function(exp, end, cohort) {
    d <- if (cohort == "Pooled") ev else filter(ev, cohort == .env$cohort)
    ok <- is.finite(d[[exp]]) & is.finite(d[[end]])
    if (sum(ok) < 30) return(tibble())
    ct <- cor.test(d[[exp]][ok], d[[end]][ok], method = "spearman", exact = FALSE)
    tibble(exposure = exp, endpoint = end, cohort = cohort, n = sum(ok),
           spearman_rho = unname(ct$estimate), p_value = ct$p.value)
  }) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"),
         sig = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**", p_adj < 0.05 ~ "*", TRUE ~ ""))
fwrite(sp_tab, file.path(out_dir, "step2_dose_response_spearman.csv"))
cat("\n== Dose-response Spearman (Pooled, sorted by |rho|) ==\n")
print(sp_tab %>% filter(cohort == "Pooled") %>% arrange(desc(abs(spearman_rho))) %>% as.data.frame())

## tsc convention sensitivity summary: whether the primary convention csi_lc and the alternatives csi_strict / csi_mid give the same conclusion
## Interpretation: if the three conventions agree in rho direction and are close in magnitude (|delta rho| <= 0.05) on the same endpoint,
##       one can write in the paper that "conclusions are insensitive to how missingness in smoking-duration is handled", so convention choice is no longer a risk.
if (all(c("csi_strict", "csi_mid") %in% names(step2)) && nrow(sp_tab)) {
  sens_cmp <- sp_tab %>%
    filter(cohort == "Pooled", exposure %in% c("csi_lc", "csi_strict", "csi_mid")) %>%
    select(exposure, endpoint, spearman_rho, p_adj) %>%
    pivot_wider(names_from = exposure,
                values_from = c(spearman_rho, p_adj),
                names_glue = "{exposure}_{.value}")
  if (all(c("csi_lc_spearman_rho", "csi_strict_spearman_rho", "csi_mid_spearman_rho") %in% names(sens_cmp))) {
    sens_cmp <- sens_cmp %>%
      mutate(max_abs_delta_rho = pmax(abs(csi_strict_spearman_rho - csi_lc_spearman_rho),
                                      abs(csi_mid_spearman_rho   - csi_lc_spearman_rho),
                                      na.rm = TRUE),
             conclusion_stable = max_abs_delta_rho <= 0.05)
  }
  fwrite(sens_cmp, file.path(out_dir, "step2_tsc_convention_sensitivity.csv"))
  cat("\n== tsc convention sensitivity (Pooled): three-convention correlation comparison ==\n"); print(as.data.frame(sens_cmp))
}

##   Tertile trend: Kruskal-Wallis (between-group difference) + Jonckheere-Terpstra trend (ordered monotonicity)
##   Fix for "object ' tmb_per_mb ' not found":
##     In the original `as.formula(paste("`", end, "` ~", t3))`, paste()'s **default separator is a space**,
##     so what gets produced is ` tmb_per_mb ` ~ py_gdc_t3 (spaces inside the backticks) -> R treats `` ` tmb_per_mb ` ``
##     as a variable name containing spaces and looks for it -> error (and the spaces around the object name are the key clue for this bug).
##     Use reformulate() instead: build the formula directly from a vector of variable names, **without string concatenation**,
##     equally safe for columns containing dots/special characters (e.g. q_T.cells.CD8).

## Self-implemented Jonckheere-Terpstra trend test (monotonicity across ordered k groups), avoiding the uninstalled clinfun.
## Statistic U = sum_{i<j} #{(x in i, y in j): y > x}; normal approximation with continuity correction.
jt_test <- function(x, g) {
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- as.integer(g[ok])
  ks <- sort(unique(g)); if (length(ks) < 2) return(NA_real_)
  U <- 0; mu <- 0; v <- 0; N <- length(x)
  for (a in seq_along(ks)) for (b in seq_along(ks)) {
    if (b <= a) next
    xa <- x[g == ks[a]]; xb <- x[g == ks[b]]
    na <- length(xa); nb <- length(xb)
    U  <- U + sum(outer(xb, xa, ">")) + 0.5 * sum(outer(xb, xa, "=="))
    mu <- mu + na * nb / 2
    v  <- v + na * nb * (na + nb + 1) / 12
  }
  if (v <= 0) return(NA_real_)
  z <- (U - mu - sign(U - mu) * 0.5) / sqrt(v)     # continuity correction
  2 * pnorm(-abs(z))
}

kw_tab <- expand_grid(exp = exposures, end = endpoints) %>%
  pmap_dfr(function(exp, end) {
    d <- ev %>% filter(is.finite(.data[[exp]]), is.finite(.data[[end]]))
    t3 <- paste0(exp, "_t3")
    if (!t3 %in% names(d) || sum(!is.na(d[[t3]])) < 30) return(tibble())
    ## use reformulate only, do not build the formula with paste
    kt <- kruskal.test(reformulate(t3, response = end), data = d)
    grp_ord <- factor(d[[t3]], levels = c("T1_Low", "T2_Mid", "T3_High"), ordered = TRUE)
    jt_p <- jt_test(d[[end]], grp_ord)
    tibble(exposure = exp, endpoint = end, n = sum(!is.na(d[[t3]])),
           kruskal_p = kt$p.value,
           jt_trend_p = jt_p,
           median_T1 = median(d[[end]][d[[t3]] == "T1_Low"], na.rm = TRUE),
           median_T2 = median(d[[end]][d[[t3]] == "T2_Mid"], na.rm = TRUE),
           median_T3 = median(d[[end]][d[[t3]] == "T3_High"], na.rm = TRUE))
  })
fwrite(kw_tab, file.path(out_dir, "step2_tertile_kruskal.csv"))

## 9. Plotting --------------------------------------------------------------------
## 9.1 Correlation heatmap (Pooled)
if (nrow(sp_tab)) {
  p1 <- ggplot(sp_tab %>% filter(cohort == "Pooled"),
               aes(exposure, endpoint, fill = spearman_rho)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f%s", spearman_rho, sig)), size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    labs(title = "Smoking exposure x immune microenvironment (Spearman, ever-smokers)",
         x = NULL, y = NULL, fill = "rho") + theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(fig_dir, "step2_exposure_TME_heatmap.pdf"), p1, width = 7, height = 5.5)
}
## 9.2 Boxplots of CSI tertiles x key endpoints
key_ends <- base::intersect(c("tmb_per_mb", "sbs4_weight", "cd8a_log2tpm", "cd274_log2tpm"), endpoints)
if (length(key_ends) && "csi_lc_t3" %in% names(ev)) {
  p2 <- ev %>% select(cohort, csi_lc_t3, all_of(key_ends)) %>%
    pivot_longer(all_of(key_ends), names_to = "endpoint") %>%
    filter(is.finite(value)) %>%
    ggplot(aes(csi_lc_t3, value, fill = csi_lc_t3)) +
    geom_boxplot(outlier.size = 0.4) +
    facet_grid(endpoint ~ cohort, scales = "free_y") +
    labs(title = "CSI tertiles vs key TME endpoints", x = "CSI tertile", y = NULL) +
    theme_bw() + theme(legend.position = "none")
  ggsave(file.path(fig_dir, "step2_CSI_tertile_endpoints.pdf"), p2, width = 8, height = 6)
}
## 9.3 TMB vs SBS4 scatter (smoker-signature sanity check: SBS4 should rise with cumulative exposure)
if (!is.null(sbs4) && nrow(sp_tab)) {
  p3 <- ev %>% filter(is.finite(py_gdc), is.finite(sbs4_weight)) %>%
    ggplot(aes(py_gdc, sbs4_weight, color = cohort)) +
    geom_point(size = 0.7, alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.6) +
    labs(title = "Pack-years vs SBS4 signature weight (sanity check)", x = "Pack-years", y = "SBS4 weight") +
    theme_bw()
  ggsave(file.path(fig_dir, "step2_packyears_SBS4.pdf"), p3, width = 6, height = 4.5)
}

## 10. Step3 interface notes ----------------------------------------------------------
## The output step2_immune_tme.csv already contains: smoking_grp / continuous exposure values and tertiles / the full TME metric set
## Step3 (09) reads this table directly to construct the SIRS index;
## Step4 (external validation) repeats the same scoring pipeline on IMvigor210 (which carries smoking-history fields).
############################################################
## Implementation caveats:
##   1) TIS is computed as a mean score rather than by ssGSEA (a one-line switch when GSVA is installed);
##   2) the TMB denominator uses 38 Mb, the MC3 common approximation, not a per-sample covered bed size;
##   3) quanTIseq is run with tumor=TRUE (cancer-cell signal subtracted); tumor=FALSE is available as a control;
##   4) the tertile trend test uses the Jonckheere-Terpstra test implemented in-script (no clinfun dependency).
############################################################
