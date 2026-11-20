## ============================================================================
## install_deconstructsigs.R -- install an hg38-capable build of deconstructSigs
##
## Why this script is needed:
##   * deconstructSigs 1.8.0 (last CRAN release, now archived) hard-codes hg19 in
##     the source (no `bsg` argument) and therefore cannot be used with hg38.
##   * deconstructSigs 1.9.0 (GitHub master) adds the `bsg` argument, but its
##     DESCRIPTION still hard-depends on BSgenome.Hsapiens.UCSC.hg19 (another
##     ~700 MB package) and CRAN has archived it, so a plain install always fails.
##
## Strategy (in priority order):
##   [1] already installed            -> skip
##   [2] reuse the retained patch source (vendor/deconstructSigs-1.9.0-hg38patched,
##       see its PATCH_NOTE.md)
##   [3] no local copy                -> download through a mirror, apply the patch
##       (drop hg19 from DESCRIPTION + hg19 -> hg38 in R/*.R), install from source
##   [4] everything failed            -> fall back: install BSgenome.Hsapiens.UCSC.hg19
##       and install the package unmodified
##
## Usage: Rscript install_deconstructsigs.R
## ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 7200)
log <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n"); flush.console() }

proj_dir <- "D:/Projects/smoking-index-ici"
## Local patched source, redistributed under vendor/.
patched <- file.path(proj_dir, "vendor/deconstructSigs-1.9.0-hg38patched")

## Result report
report <- function(ok) {
  cat("\n================ result ================\n")
  cat("deconstructSigs: ", if (ok) paste0("OK ", as.character(packageVersion("deconstructSigs")))
      else "not installed", "\n")
  cat("BSgenome.Hsapiens.UCSC.hg38: ",
      requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE), "\n")
  cat("========================================\n")
}

## ---------------------------------------------------------------- [1] installed?
if (requireNamespace("deconstructSigs", quietly = TRUE)) {
  log("deconstructSigs already installed (", as.character(packageVersion("deconstructSigs")), "), skipping")
  report(TRUE); quit(save = "no")
}
if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE))
  stop("Please run install_deps.R first to install BSgenome.Hsapiens.UCSC.hg38")

## ---------------------------------------------------------------- patch function
patch_pkg <- function(pkg) {
  desc_f <- file.path(pkg, "DESCRIPTION")
  d <- readLines(desc_f, warn = FALSE)
  d <- gsub("\\s*BSgenome\\.Hsapiens\\.UCSC\\.hg19\\s*,?", "", d)
  d <- gsub(",\\s*,", ",", d); d <- gsub(",\\s*$", "", d)
  writeLines(d, desc_f)
  n <- 0
  for (f in list.files(file.path(pkg, "R"), pattern = "\\.R$", full.names = TRUE)) {
    x <- readLines(f, warn = FALSE)
    y <- gsub("BSgenome\\.Hsapiens\\.UCSC\\.hg19", "BSgenome.Hsapiens.UCSC.hg38", x)
    if (!identical(x, y)) { writeLines(y, f); n <- n + sum(x != y) }
  }
  log("patch applied: hg19 dependency removed from DESCRIPTION; ",
      n, " line(s) changed hg19 -> hg38 in R/")
}

## ---------------------------------------------------------------- [2] local patched copy
ok <- FALSE
if (dir.exists(patched)) {
  log("reusing local patched source: ", basename(patched))
  r <- try(install.packages(patched, repos = NULL, type = "source"), silent = TRUE)
  ok <- requireNamespace("deconstructSigs", quietly = TRUE)
}

## ---------------------------------------------------------------- [3] mirror download + patch
if (!ok) {
  log("no usable local copy -> downloading the source through a mirror ...")
  tmp <- file.path(tempdir(), "dsig_src"); unlink(tmp, recursive = TRUE)
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  urls <- c(
    "https://gh-proxy.com/https://github.com/raerose01/deconstructSigs/archive/refs/heads/master.tar.gz",
    "https://ghfast.top/https://github.com/raerose01/deconstructSigs/archive/refs/heads/master.tar.gz",
    "https://codeload.github.com/raerose01/deconstructSigs/tar.gz/refs/heads/master"
  )
  for (u in urls) {
    dest <- file.path(tmp, "src.tar.gz")
    r <- try(download.file(u, dest, mode = "wb", quiet = TRUE, method = "libcurl"), silent = TRUE)
    if (inherits(r, "try-error") || !file.exists(dest) || file.size(dest) < 1e5) next
    if (!inherits(try(untar(dest, exdir = tmp), silent = TRUE), "try-error")) {
      pkg <- list.files(tmp, pattern = "^deconstructSigs-", full.names = TRUE)[1]
      if (!is.na(pkg)) {
        log("download succeeded (", basename(u), "), applying the patch ...")
        patch_pkg(pkg)
        try(install.packages(pkg, repos = NULL, type = "source"), silent = TRUE)
        ok <- requireNamespace("deconstructSigs", quietly = TRUE)
        if (ok) break
      }
    }
  }
}

## ---------------------------------------------------------------- [4] fallback: install hg19
if (!ok) {
  log("patch route failed -> fallback: install BSgenome.Hsapiens.UCSC.hg19 (~700 MB), then install unmodified")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", ask = FALSE, update = FALSE, type = "source")
  if (requireNamespace("remotes", quietly = TRUE))
    try(remotes::install_github("raerose01/deconstructSigs", upgrade = "never", quiet = FALSE))
  ok <- requireNamespace("deconstructSigs", quietly = TRUE)
}

report(ok)
