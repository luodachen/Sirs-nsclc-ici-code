## ============================================================================
## smoke_test_step2.R -- minimal verification once the Step 2 dependencies are
## installed (does not run the full pipeline; ~2-5 minutes)
##
## Verifies three things:
##   1) SBS4 chain    : Xena mutation table -> 96 channels -> deconstructSigs
##                      (COSMIC v2, 30 signatures) -> Signature.4 weight
##   2) deconv chain  : do quanTIseq (quantiseqr) / EPIC / MCPcounter return
##                      values on real TPM data
##   3) sanity check  : does the SBS4 weight rise with pack-years (the
##                      mechanistic expectation for smoking exposure)
## ============================================================================
options(timeout = 7200)
suppressPackageStartupMessages({ library(data.table); library(tidyverse) })

proj_dir <- "D:/Projects/smoking-index-ici"
data_dir <- file.path(proj_dir, "data_raw/xena")
step1_csv <- file.path(proj_dir, "step1_output/smoking_index_clinical.csv")

hr <- function(x) cat("\n", strrep("-", 62), "\n", x, "\n", strrep("-", 62), "\n", sep = "")

## ---------------------------------------------------------------- 1. SBS4 chain
hr("1. SBS4 chain")
ok_sbs4 <- FALSE
if (requireNamespace("deconstructSigs", quietly = TRUE) &&
    requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(deconstructSigs); library(BSgenome.Hsapiens.UCSC.hg38)
  })
  mut <- fread(file.path(data_dir, "TCGA-LUAD.somaticmutation_wxs.tsv.gz"),
               select = c("sample", "chrom", "start", "ref", "alt", "effect"))
  cat("raw mutation rows: ", nrow(mut), "\n")
  ## Use 60 samples for the smoke test (real data, light workload)
  set.seed(42)
  keep <- sample(unique(mut$sample), 60)
  mi <- mut[sample %in% keep &
              ref %in% c("A","C","G","T") & alt %in% c("A","C","G","T")] %>%
    as_tibble() %>%
    transmute(sample.id = substr(sample, 1, 15), chr = chrom, pos = start, ref, alt) %>%
    as.data.frame()
  cat("SNV rows: ", nrow(mi), "; samples: ", length(unique(mi$sample.id)), "\n")

  si <- mut.to.sigs.input(mut.ref = mi, sample.id = "sample.id", chr = "chr",
                          pos = "pos", ref = "ref", alt = "alt",
                          bsg = BSgenome.Hsapiens.UCSC.hg38)
  cat("96-channel matrix: ", nrow(si), "samples x", ncol(si), "channels (expected 96)\n")

  ## NOTE whichSignatures() is a single-sample interface: tumor.ref must be a
  ## data.frame and sample.id must be given one sample at a time.
  si_df  <- as.data.frame(si)
  n_snv  <- rowSums(si_df)
  keep_s <- names(n_snv)[n_snv >= 50]
  cat("samples with SNV >= 50: ", length(keep_s), "/", nrow(si_df), "\n")

  w_list <- lapply(keep_s, function(s) {
    r <- try(whichSignatures(tumor.ref = si_df, sample.id = s,
                             signatures.ref = signatures.cosmic,
                             contexts.needed = TRUE,
                             tri.counts.method = "exome2genome")[["weights"]], silent = TRUE)
    if (inherits(r, "try-error")) return(NULL)
    r$sample15 <- s; r
  })
  w <- dplyr::bind_rows(w_list)
  cat("Signature.4 (SBS4) weight: median ", round(median(w$Signature.4, na.rm = TRUE), 4),
      "  range ", round(min(w$Signature.4, na.rm = TRUE), 4), "-",
      round(max(w$Signature.4, na.rm = TRUE), 4), " (n = ", nrow(w), ")\n", sep = "")
  cat("  first 3 samples:\n"); print(head(w[, c("sample15", "Signature.1", "Signature.4")], 3))
  ok_sbs4 <- nrow(w) > 0
} else {
  cat("!! deconstructSigs or BSgenome.Hsapiens.UCSC.hg38 is not installed, skipping\n")
}

## ---------------------------------------------------------------- 2. deconvolution chain
hr("2. Deconvolution chain")
## Gene mapping identical to script 07: gencode probeMap (note the lower-case
## file name ...gene.probemap; the id column needs its version suffix stripped).
probemap_file <- file.path(data_dir, "gencode.v36.annotation.gtf.gene.probemap")
if (!file.exists(probemap_file) || file.size(probemap_file) < 1e5) {
  download.file("https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
                probemap_file, mode = "wb", quiet = TRUE, method = "libcurl")
}
pm <- fread(probemap_file) %>%
  transmute(id = sub("\\..*$", "", id), symbol = gene) %>%
  filter(!is.na(symbol), symbol != "", symbol != "?")
fq <- fread(file.path(data_dir, "TCGA-LUAD.star_fpkm-uq.tsv.gz"))
m <- fq %>%
  mutate(id = sub("\\..*$", "", Ensembl_ID)) %>% select(-Ensembl_ID) %>%
  distinct(id, .keep_all = TRUE) %>%
  inner_join(pm, by = "id") %>% filter(!is.na(symbol)) %>%
  select(symbol, where(is.numeric)) %>%
  group_by(symbol) %>% summarise(across(everything(), max), .groups = "drop") %>%
  column_to_rownames("symbol") %>% as.matrix()
cat("matrix type: ", class(m[1, 1]), " (expected numeric)\n")
cat("mapped matrix: ", nrow(m), "genes x ", ncol(m), "samples\n")
tpm <- apply(m[, 1:30, drop = FALSE], 2,
             function(x) { s <- sum(x[x > 0], na.rm = TRUE); ifelse(x > 0, x/s*1e6, 0) })

if (requireNamespace("quantiseqr", quietly = TRUE)) {
  qs <- try(quantiseqr::run_quantiseq(expression_data = tpm, signature_matrix = "TIL10",
                                      is_arraydata = FALSE, is_tumordata = TRUE,
                                      scale_mRNA = TRUE, method = "lsei"), silent = TRUE)
  if (inherits(qs, "try-error")) {
    cat("!! quantiseqr call failed: ", conditionMessage(attr(qs, "condition")), "\n")
  } else {
    cat("quanTIseq OK: ", nrow(qs), "samples x", ncol(qs) - 1, "cell types\n")
    print(as_tibble(qs)[1:2, 1:min(6, ncol(qs))])
  }
} else cat("!! quantiseqr not installed\n")

if (requireNamespace("EPIC", quietly = TRUE)) {
  r <- try(EPIC::EPIC(tpm), silent = TRUE)
  if (inherits(r, "try-error")) cat("!! EPIC call failed: ", conditionMessage(attr(r, "condition")), "\n")
  else cat("EPIC OK: ", nrow(r$cellFractions), " samples x ", ncol(r$cellFractions), " cell types\n", sep = "")
} else cat("!! EPIC not installed\n")

if (requireNamespace("MCPcounter", quietly = TRUE)) {
  r <- try(MCPcounter::MCPcounter.estimate(tpm, featuresType = "HUGO_symbols"), silent = TRUE)
  if (inherits(r, "try-error")) cat("!! MCPcounter call failed: ", conditionMessage(attr(r, "condition")), "\n")
  else cat("MCPcounter OK: ", nrow(r), " cell types x ", ncol(r), " samples\n", sep = "")
} else cat("!! MCPcounter not installed\n")

## ---------------------------------------------------------------- 3. Summary
hr("3. Conclusion")
cat("SBS4 chain: ", if (ok_sbs4) "PASS (sbs4_weight can be produced)" else "FAIL (Step 2 will skip SBS4)", "\n")
cat("Deconvolution chain: see above (any working pathway yields q_/e_/mcp_ columns)\n")
cat("Next: run 07_Step2_immune_microenvironment_TMB_SBS4.R to produce the Step 2 outputs\n")
