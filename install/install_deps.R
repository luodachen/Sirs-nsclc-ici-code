## ============================================================================
## install_deps.R -- install the Step 2 dependencies (idempotent; safe to re-run)
##
## Two dependency chains are installed:
##   [SBS4 chain]    BSgenome + BSgenome.Hsapiens.UCSC.hg38 + deconstructSigs
##   [deconv chain]  quantiseqr (Bioconductor quanTIseq) + EPIC + MCPcounter
##
## Usage:  Rscript install_deps.R
## Env:    R 4.5.1 / Bioconductor 3.21 / no Rtools, so binary packages are
##         preferred and only pure-R packages are built from source.
## ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 7200, warn = 1)

log <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush.console() }
need <- function(p) !requireNamespace(p, quietly = TRUE)

log("R:", R.version.string)
log("libPaths:", paste(.libPaths(), collapse = " | "))

## ---------------------------------------------------------------- 1. BiocManager
if (need("BiocManager")) install.packages("BiocManager")
suppressMessages(library(BiocManager))
log("BiocManager", as.character(packageVersion("BiocManager")),
    "-> Bioconductor", as.character(BiocManager::version()))

if (need("remotes")) install.packages("remotes")

## ---------------------------------------------------------------- 2. SBS4 chain
## 2.1 BSgenome (binary, no compilation)
if (need("BSgenome")) {
  log("installing BSgenome ...")
  BiocManager::install("BSgenome", ask = FALSE, update = FALSE)
}
if (need("GenomeInfoDb")) {
  log("installing GenomeInfoDb ...")
  BiocManager::install("GenomeInfoDb", ask = FALSE, update = FALSE)
}

## 2.2 hg38 reference genome (annotation source package, ~731 MB -- the slowest step)
if (need("BSgenome.Hsapiens.UCSC.hg38")) {
  log("installing BSgenome.Hsapiens.UCSC.hg38 (~731 MB, this takes a while) ...")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg38",
                       ask = FALSE, update = FALSE, type = "source")
} else {
  log("BSgenome.Hsapiens.UCSC.hg38 already present, skipping")
}

## 2.3 deconstructSigs (archived on CRAN since 2024)
##     GitHub master (= 1.9.x) is preferred; on failure, fall back to CRAN Archive 1.8.0.
##     NOTE: 1.8.0 hard-codes hg19 and cannot be used with hg38 coordinates.
##     The hg38-capable build must come from install_deconstructsigs.R (patched source).
if (need("reshape2")) install.packages("reshape2")
if (need("plyr"))     install.packages("plyr")
if (need("deconstructSigs")) {
  log("installing deconstructSigs (GitHub raerose01/master) ...")
  try(remotes::install_github("raerose01/deconstructSigs",
                              upgrade = "never", quiet = FALSE), silent = TRUE)
  if (need("deconstructSigs")) {
    log("GitHub install failed -> falling back to CRAN Archive 1.8.0 ...")
    install.packages("https://cran.r-project.org/src/contrib/Archive/deconstructSigs/deconstructSigs_1.8.0.tar.gz",
                     repos = NULL, type = "source")
  }
} else {
  log("deconstructSigs already present, skipping")
}

## ---------------------------------------------------------------- 3. Deconvolution chain
## 3.1 quantiseqr: the Bioconductor implementation of quanTIseq (avoids the heavy
##     immunedeconv dependency tree)
if (need("quantiseqr")) {
  log("installing quantiseqr (Bioc, official quanTIseq implementation) ...")
  BiocManager::install("quantiseqr", ask = FALSE, update = FALSE)
}
## 3.2 EPIC (GitHub)
if (need("EPIC")) {
  log("installing EPIC (GitHub GfellerLab/EPIC) ...")
  try(remotes::install_github("GfellerLab/EPIC", upgrade = "never", quiet = FALSE), silent = TRUE)
}
## 3.3 MCP-counter (GitHub; the package body lives in the Source/ subdirectory)
if (need("MCPcounter")) {
  log("installing MCPcounter (GitHub ebecht/MCPcounter, subdir = Source) ...")
  try(remotes::install_github("ebecht/MCPcounter", subdir = "Source",
                              upgrade = "never", quiet = FALSE), silent = TRUE)
}

## ---------------------------------------------------------------- 4. Summary
cat("\n================ install summary ================\n")
for (p in c("BSgenome", "BSgenome.Hsapiens.UCSC.hg38", "deconstructSigs",
            "quantiseqr", "EPIC", "MCPcounter", "GSVA", "immunedeconv")) {
  cat(sprintf("%-32s %s\n", p,
              if (requireNamespace(p, quietly = TRUE))
                paste0("OK  ", as.character(packageVersion(p))) else "missing"))
}
cat("=================================================\n")
