## ============================================================================
## 12_Step4c_single_algorithm_sensitivity.R
## Algorithm-sensitivity analysis
## ============================================================================
## Purpose: test whether the SIRS construction, and therefore the IC0 finding,
##   is an artefact of MIXING three deconvolution algorithms.
##
##   The main SIRS combines eight features drawn from three algorithms
##   (MCP-counter 3, EPIC 3, quanTIseq 2). Here SIRS is rebuilt from ONE
##   algorithm at a time under the identical construction rule:
##       partial Spearman rho between feature and anchor (csi_lc / sbs4_abs),
##       adjusting for age/gender/stage/cohort; BH FDR < 0.10; weight = rho;
##       each feature keeps the anchor with the larger |rho|.
##
## Design notes
##   (a) The BH adjustment in 09_Step3 is applied WITHIN each anchor
##       (m = number of features per anchor), not across anchors x features.
##       This is replicated exactly: the full-pool re-derivation reproduces the
##       Step3 weights (max |delta| = 4.5e-07).
##   (b) Under the same rule the EPIC pool alone selects 3 features
##       (e_CD4_Tcells, e_NKcells, e_Endothelial) - exactly EPIC's share of the
##       main SIRS - so a genuine EPIC-only index exists and can be tested.
##   (c) The EPIC-only feature set is identical under three independent
##       definitions (pool-FDR / top-3 / full-pool sub-score), and no pool
##       feature is degenerate (max zero fraction 0.00), so a
##       degeneracy-filtered variant would be identical to it. Both checks are
##       reported as guards.
##
##   Indices: main (mixed) | MCP-only | quanTIseq-only | EPIC-only.
##
## Input:
##   step2_output/step2_immune_tme.csv       (TCGA-NSCLC; all three feature blocks)
##   data/IMvigor210/deconv_imvigor.rds      (cached: EPIC 348x8, MCP 10x348, quanTIseq 348x12)
##   step3_output/step3_SIRS_all_samples.csv (main SIRS, TCGA side)
##   step4_output/step4_analysis_table.csv   (348 rows; response, IC level, main SIRS)
##
## Output (all under step4_output/):
##   step4c_weights_by_algorithm.csv
##   step4c_tcga_exposure_validity.csv
##   step4c_imvig_auc_overall.csv
##   step4c_imvig_ic_stratified_auc.csv
##   step4c_ic0_headtohead_delong.csv
##   step4c_complementarity.csv
##   figs/step4c_algorithm_ic_stratified_auc.pdf
##   12_Step4c_summary.md
##
## Conventions honoured
##   - ASCII identifiers only; English comments; no non-ASCII inside figures.
##   - Offline environment: no new packages. DeLong (1988), bootstrap and
##     permutation tests are implemented from their original definitions.
##   - Rscript cannot read Chinese file names -> ASCII names only.
## ============================================================================

## 0. Environment -------------------------------------------------------------
## Author: see manuscript author list   Date: 2026-10-17
proj_dir <- "D:/Projects/smoking-index-ici"
setwd(proj_dir)
out_dir  <- file.path(proj_dir, "step4_output")
fig_dir  <- file.path(out_dir, "figs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({ library(data.table) })
options(stringsAsFactors = FALSE)
set.seed(20261010)

FDR_MAIN <- 0.10
ANCHORS  <- c("csi_lc", "sbs4_abs")
CONF     <- c("age", "gender_m", "stage_n", "cohort_n")
QS_EXCLUDE <- c("q_T.cells.CD4", "q_Monocytes")   # zero-inflated, dropped in Step3 (C2)

cat("========== Step4c: single-algorithm sensitivity of SIRS ==========\n\n")

## ============================================================================
## 1. Shared analytics (identical estimators to 11_Step4b) -------------------
## ============================================================================
## 1.1 Confounder-adjusted partial Spearman (numerically equivalent to
##     ppcor::pcor.test; validated in Step2b). Verbatim from 09_Step3.
partial_spearman <- function(x, y, Z) {
  ok <- is.finite(x) & is.finite(y) & stats::complete.cases(Z)
  x <- x[ok]; y <- y[ok]; Z <- as.data.frame(Z)[ok, , drop = FALSE]
  rx <- rank(x); ry <- rank(y)
  Z  <- as.data.frame(lapply(Z, function(z) rank(z)))
  res <- stats::residuals(stats::lm(cbind(rx, ry) ~ as.matrix(Z)))
  r <- stats::cor(res[, 1], res[, 2])
  n <- length(rx); k <- ncol(Z); df <- n - k - 2
  tt <- r * sqrt(df / (1 - r^2))
  list(rho = r, n = n, p_value = 2 * stats::pt(-abs(tt), df = df))
}

## 1.2 Mann-Whitney AUC (ties = 0.5)
auc_mw <- function(x, y) {
  ok <- is.finite(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2 || length(unique(x)) < 2) return(NA_real_)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(rank(x)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## 1.3 Bootstrap AUC interval (stratified resampling)
boot_auc <- function(x, y, B = 1000, seed = 20261010) {
  ok <- is.finite(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2 || length(unique(x)) < 2)
    return(c(auc = NA_real_, lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  i1 <- which(y == 1); i0 <- which(y == 0)
  bs <- replicate(B, {
    s <- c(sample(i1, length(i1), TRUE), sample(i0, length(i0), TRUE))
    if (length(unique(y[s])) < 2) NA_real_ else auc_mw(x[s], y[s])
  })
  bs <- bs[is.finite(bs)]
  c(auc = auc_mw(x, y), lo = unname(stats::quantile(bs, .025)),
    hi = unname(stats::quantile(bs, .975)))
}

## 1.4 Permutation test for AUC != 0.5
perm_auc <- function(x, y, B = 1000, seed = 20261010) {
  ok <- is.finite(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  if (length(unique(y)) < 2 || length(unique(x)) < 2) return(c(auc = NA_real_, p_emp = NA_real_))
  stat <- auc_mw(x, y); set.seed(seed)
  null <- replicate(B, auc_mw(x, sample(y)))
  null <- null[is.finite(null)]
  c(auc = stat, p_emp = (sum(abs(null - 0.5) >= abs(stat - 0.5)) + 1) / (length(null) + 1))
}

## 1.5 DeLong (1988) correlated-ROC components on pairwise-complete observations
delong_parts <- function(x, y) {
  ok <- is.finite(x) & !is.na(y); x <- x[ok]; y <- y[ok]
  X <- x[y == 1]; Y <- x[y == 0]; m <- length(X); n <- length(Y)
  psi <- outer(X, Y, function(a, b) (a > b) + 0.5 * (a == b))
  list(auc = mean(psi), V10 = rowMeans(psi), V01 = colMeans(psi), m = m, n = n)
}
delong_var <- function(p) stats::var(p$V10) / p$m + stats::var(p$V01) / p$n
delong_cov <- function(p1, p2) stats::cov(p1$V10, p2$V10) / p1$m + stats::cov(p1$V01, p2$V01) / p1$n
delong_p_paired <- function(x1, x2, y) {
  ok <- is.finite(x1) & is.finite(x2) & !is.na(y)
  if (sum(y[ok] == 1) < 2 || sum(y[ok] == 0) < 2)
    return(c(p = NA_real_, n = sum(ok), d = NA_real_))
  p1 <- delong_parts(x1[ok], y[ok]); p2 <- delong_parts(x2[ok], y[ok])
  v <- delong_var(p1) + delong_var(p2) - 2 * delong_cov(p1, p2)
  c(p = if (!is.finite(v) || v <= 0) NA_real_ else
      2 * stats::pnorm(-abs((p1$auc - p2$auc) / sqrt(v))),
    n = sum(ok), d = p1$auc - p2$auc)
}

## ============================================================================
## 2. TCGA side: analysis set, screening pools, index construction -----------
## ============================================================================
cat("================ 2. TCGA side: analysis set + per-algorithm weights ================\n")
step2 <- data.table::fread(file.path(proj_dir, "step2_output/step2_immune_tme.csv"),
                           na.strings = c("", "NA"))
stopifnot(length(grep("^e_",   names(step2))) == 8L,
          length(grep("^mcp_", names(step2))) == 10L,
          length(grep("^q_",   names(step2))) == 11L)   # Guard 0: all three blocks present

st <- step2
st[, os_event := vital == "Dead"]
st[, os_time  := ifelse(vital == "Dead", days_death, days_fu)]
st[!is.finite(os_time) & is.finite(days_fu),    os_time := days_fu]
st[!is.finite(os_time) & is.finite(days_death), os_time := days_death]
st[, age      := age_y / 365.25]
## Stage collapsing with first-match-wins semantics (equivalent to case_when);
## sequential grepl overwriting would let "^Stage I" consume "Stage III/IV"
st[, stage_grp := ifelse(grepl("^Stage\\s*IV",  stage), "IV",
                  ifelse(grepl("^Stage\\s*III", stage), "III",
                  ifelse(grepl("^Stage\\s*II",  stage), "II",
                  ifelse(grepl("^Stage\\s*I",   stage), "I", NA_character_))))]
st[, sbs4_abs := sbs4_weight * tmb_per_mb]
st[, os_time  := pmax(os_time, 0.5)]
st[, gender_m := as.integer(factor(gender)) - 1L]
st[, stage_n  := as.integer(factor(stage_grp, levels = c("I", "II", "III", "IV")))]
st[, cohort_n := as.integer(factor(cohort)) - 1L]

ev <- st[smoking_grp %in% c("current", "former")]
dm <- as.data.frame(ev[vital %in% c("Alive", "Dead") & is.finite(os_time) &
                         !is.na(age) & !is.na(gender_m) & !is.na(stage_n) & !is.na(cohort_n)])
cat(sprintf("TCGA modelling set: n = %d, events = %d\n", nrow(dm), sum(dm$os_event)))
## Guard 1: reproduce Step3's modelling set. Verified by an explicit side-by-side
## against the verbatim dplyr construction in 09_Step3 (identical sample set).
## NB: n = 855 / events = 341 in step3_SIRS_formula.md refers to `use` (modelling
## set further restricted to complete SIRS features), NOT to `dm`.
stopifnot(nrow(dm) == 862L, sum(dm$os_event) == 345L)

## --- screening engines ---
feature_pool <- function(dat, algo) {
  f <- switch(algo,
              EPIC = grep("^e_",   names(dat), value = TRUE),
              MCP  = grep("^mcp_", names(dat), value = TRUE),
              QS   = base::setdiff(grep("^q_", names(dat), value = TRUE), QS_EXCLUDE),
              FULL = c("tis_score", "cyt_score", "cd274_log2tpm", "cd8a_log2tpm",
                       grep("^e_", names(dat), value = TRUE),
                       grep("^mcp_", names(dat), value = TRUE),
                       base::setdiff(grep("^q_", names(dat), value = TRUE), QS_EXCLUDE)))
  f[sapply(f, function(v) sum(is.finite(dat[[v]])) >= 0.8 * nrow(dat))]
}
screen_pool <- function(dat, fp) {
  ## IMPORTANT: in 09_Step3 the BH adjustment is applied WITHIN each anchor
  ## (inside the per-anchor map_dfr), i.e. m = number of features per anchor,
  ## NOT m = features x anchors. Reproducing this exactly is essential: pooling
  ## the two anchors into one BH family silently changes the selection (it drops
  ## the 8-feature main SIRS to 3 features). Verified against Step3 below.
  do.call(rbind, lapply(ANCHORS, function(a) {
    d0 <- dat[is.finite(dat[[a]]), ]
    r <- do.call(rbind, lapply(fp, function(f) {
      d <- d0[is.finite(d0[[f]]), ]
      if (nrow(d) < 30) return(NULL)
      pr <- partial_spearman(d[[a]], d[[f]], d[, CONF, drop = FALSE])
      data.frame(feature = f, anchor = a, n = pr$n, rho = pr$rho,
                 p_value = pr$p_value, stringsAsFactors = FALSE)
    }))
    r$p_adj <- stats::p.adjust(r$p_value, method = "BH")
    r
  }))
}
best_anchor <- function(scr) {   # one row per feature: the anchor with larger |rho|
  do.call(rbind, lapply(unique(scr$feature), function(f) {
    s <- scr[scr$feature == f, ]; s[which.max(abs(s$rho)), ]
  }))
}
w_fdr  <- function(scr, fdr = FDR_MAIN) { s <- best_anchor(scr); s <- s[s$p_adj < fdr, ]; s[order(-abs(s$rho)), ] }
w_topk <- function(scr, k) { s <- best_anchor(scr); s <- s[order(-abs(s$rho)), ]; head(s, k) }

## --- Step3-identical full pool: re-derive the main weights and guard them ---
scr_full <- screen_pool(dm, feature_pool(dm, "FULL"))
w_full   <- w_fdr(scr_full)
W_MAIN   <- c("mcp_Myeloid dendritic cells", "mcp_Neutrophils", "mcp_Endothelial cells",
              "q_Macrophages.M2", "q_T.cells.CD8", "e_CD4_Tcells", "e_NKcells", "e_Endothelial")
W_MAIN_W <- c(-0.179523, -0.138195, -0.083697, -0.103068, 0.079900, -0.091021, 0.082585, -0.080418)
cat(sprintf("\nFull-pool screening (Step3 rule): %d features in pool -> %d selected at FDR < %.2f\n",
            length(feature_pool(dm, "FULL")), nrow(w_full), FDR_MAIN))
dmax <- max(abs(w_full$rho[match(W_MAIN, w_full$feature)] - W_MAIN_W))
cat(sprintf("Guard 2: full-pool re-derivation reproduces the Step3 weights (max |delta| = %.2e)\n", dmax))
stopifnot(setequal(w_full$feature, W_MAIN), dmax < 1e-4)

## --- per-algorithm screens and index definitions ---
scr_by_algo <- lapply(c("EPIC", "MCP", "QS"),
                      function(a) screen_pool(dm, feature_pool(dm, a)))
names(scr_by_algo) <- c("EPIC", "MCP", "QS")

## Degeneracy check (a near-constant feature would be amplified by z-standardisation
## and could make an index fail for the wrong reason). Reported explicitly, whether
## or not it changes anything.
zf <- lapply(c("EPIC", "MCP", "QS"), function(a) {
  f <- feature_pool(dm, a)
  data.frame(algorithm = a, feature = f,
             zero_frac = sapply(f, function(v) mean(dm[[v]] == 0, na.rm = TRUE)),
             stringsAsFactors = FALSE)
})
zf <- do.call(rbind, zf)
cat(sprintf("\n-- degeneracy check: max zero fraction across all pool features = %.3f --\n", max(zf$zero_frac)))
cat(sprintf("   (features with >50%% zeros: %d -> a degeneracy-filtered EPIC variant would be %s)\n",
            sum(zf$zero_frac > 0.5),
            if (sum(zf$zero_frac[zf$algorithm == "EPIC"] > 0.5)) "different" else "IDENTICAL"))

idx_make <- function(nm, w) {
  if (is.null(w) || !nrow(w)) return(NULL)
  list(name = nm, feats = w$feature, weights = stats::setNames(w$rho, w$feature))
}
IDX <- list()
IDX[["main_mixed"]]  <- list(name = "Main SIRS (mixed)",
                             feats = W_MAIN, weights = stats::setNames(W_MAIN_W, W_MAIN))
if (!is.null(i <- idx_make("MCP-only", w_fdr(scr_by_algo$MCP))))      IDX[["MCP"]]  <- i
if (!is.null(i <- idx_make("quanTIseq-only", w_fdr(scr_by_algo$QS)))) IDX[["QS"]]   <- i
IDX[["EPIC"]] <- idx_make("EPIC-only", w_fdr(scr_by_algo$EPIC))

## The EPIC-only definition is derived three independent ways; if they coincide the
## index is not an artefact of one particular selection rule. Reported as a guard.
epic_defs <- list(pool_FDR  = IDX$EPIC$feats,
                  top3      = w_topk(scr_by_algo$EPIC, 3)$feature,
                  sub_score = w_full$feature[grepl("^e_", w_full$feature)])
epic_same <- length(unique(lapply(epic_defs, sort))) == 1L
cat(sprintf("\nEPIC-only defined three independent ways (pool-FDR / top-3 / full-pool sub-score);\n  identical feature set: %s -> %s\n",
            epic_same, paste(sort(epic_defs$pool_FDR), collapse = ", ")))
stopifnot(epic_same)

## Symmetric follow-up: use EVERY feature of each algorithm (no FDR filter at all,
## weight = the feature's own adjusted partial rho). This pre-empts the objection
## that a single-algorithm index failed because the "wrong" features were selected.
for (a in c("EPIC", "MCP", "QS")) {
  ba <- best_anchor(scr_by_algo[[a]])
  IDX[[paste0(a, "_all")]] <- list(
    name   = paste0(c(EPIC = "EPIC", MCP = "MCP-counter", QS = "quanTIseq")[a],
                    "-only (all features)"),
    feats  = ba$feature,
    weights = stats::setNames(ba$rho, ba$feature))
}

cat("\n-- indices constructed (identical rule, one algorithm at a time) --\n")
w_tab <- do.call(rbind, lapply(names(IDX), function(k) {
  x <- IDX[[k]]; if (is.null(x)) return(NULL)
  data.frame(index = x$name, algorithm = k, n_features = length(x$feats),
             features = paste(x$feats, collapse = "; "),
             weights = paste(sprintf("%+.4f", x$weights), collapse = "; "),
             stringsAsFactors = FALSE)
}))
print(w_tab[, c("index", "n_features", "features")], row.names = FALSE)
n_epic <- if (is.null(IDX$EPIC)) 0L else length(IDX$EPIC$feats)
if (n_epic > 0) {
  cat(sprintf("\nEPIC-only under the SAME rule selects %d feature(s): the index IS formable,\n", n_epic))
  cat("  so the IC0 result can be tested on a genuine EPIC-only index.\n")
} else {
  cat("\nEPIC-only under the SAME rule selects 0 features: the index is NOT formable.\n")
}
data.table::fwrite(w_tab, file.path(out_dir, "step4c_weights_by_algorithm.csv"))

## --- TCGA-side scoring + exposure validity + correlation with the main SIRS ---
scorer_tcga <- function(dat, feats, w) {          # train-parameter z on the complete-case subset
  ok <- rep(TRUE, nrow(dat)); for (f in feats) ok <- ok & is.finite(dat[[f]])
  X  <- as.matrix(dat[ok, feats, drop = FALSE])
  sdv <- apply(X, 2, stats::sd)
  if (any(!is.finite(sdv) | sdv <= 0)) stop("Zero-variance TCGA feature")
  Z <- scale(X, center = colMeans(X), scale = sdv)
  out <- rep(NA_real_, nrow(dat)); out[ok] <- as.numeric(Z %*% w); out
}
main_s <- data.table::fread(file.path(proj_dir, "step3_output/step3_SIRS_all_samples.csv"))
tc <- data.frame(sample15 = dm$sample15, main = main_s$sirs[match(dm$sample15, main_s$sample15)])
## Guard 3: the TCGA-side main SIRS recomputed here must match the Step3 export.
## Tolerance 1e-4: the export CSV stores `sirs` as rounded text. Verified that the
## residual is a writing-precision artefact, not a computational one -- recomputing
## with the published 6-dp standardisation constants gives max|delta| = 4.8e-04,
## i.e. WORSE than this full-precision path (5.1e-06).
tc_main_re <- scorer_tcga(dm, W_MAIN, W_MAIN_W)
okm <- is.finite(tc$main) & is.finite(tc_main_re)
cat(sprintf("\nGuard 3: TCGA main SIRS vs Step3 export: n = %d, max |delta| = %.2e, rho = %.4f\n",
            sum(okm), max(abs(tc$main[okm] - tc_main_re[okm])),
            stats::cor(tc$main[okm], tc_main_re[okm], method = "spearman")))
stopifnot(max(abs(tc$main[okm] - tc_main_re[okm])) < 1e-4)

val_rows <- list()
for (k in names(IDX)) {
  x <- IDX[[k]]; if (is.null(x)) next
  sc <- scorer_tcga(dm, x$feats, unname(x$weights))
  tc[[k]] <- sc
  for (v in c("csi_lc", "sbs4_abs", "tmb_per_mb")) {
    o <- is.finite(dm[[v]]) & is.finite(sc); if (sum(o) < 30) next
    ct <- suppressWarnings(stats::cor.test(dm[[v]][o], sc[o], method = "spearman"))
    val_rows[[length(val_rows) + 1]] <- data.frame(index = x$name, algorithm = k,
      variable = v, n = sum(o), rho = unname(ct$estimate), p_value = ct$p.value,
      stringsAsFactors = FALSE)
  }
}
val_tab <- do.call(rbind, val_rows)     # raw p only: these are descriptive validity
                                        # checks across different indices, not one family
cat("\n-- TCGA exposure validity of each index --\n")
print(val_tab, row.names = FALSE, digits = 4)
data.table::fwrite(val_tab, file.path(out_dir, "step4c_tcga_exposure_validity.csv"))

corr_tcga <- data.frame()
for (k in names(IDX)) {
  x <- IDX[[k]]; if (is.null(x) || k == "main_mixed") next
  o <- is.finite(tc$main) & is.finite(tc[[k]])
  ct <- suppressWarnings(stats::cor.test(tc$main[o], tc[[k]][o], method = "spearman"))
  corr_tcga <- rbind(corr_tcga, data.frame(index = x$name, algorithm = k,
    n = sum(o), rho_vs_main = unname(ct$estimate), stringsAsFactors = FALSE))
}
cat("\n-- TCGA: index vs main SIRS --\n"); print(corr_tcga, row.names = FALSE, digits = 4)

## EPIC non-convergence sensitivity
n_nc <- sum(dm$epic_nonconverged %in% TRUE)
w_epic_nc <- w_fdr(screen_pool(dm[!(dm$epic_nonconverged %in% TRUE), ],
                               feature_pool(dm, "EPIC")))
cat(sprintf("\n-- EPIC non-convergence: %d of %d modelling-set samples flagged; EPIC-only top3 unchanged: %s --\n",
            n_nc, nrow(dm),
            setequal(w_topk(scr_by_algo$EPIC, 3)$feature, w_topk(
              screen_pool(dm[!(dm$epic_nonconverged %in% TRUE), ], feature_pool(dm, "EPIC")), 3)$feature)))

## ============================================================================
## 3. IMvigor210 side: per-algorithm feature matrices ------------------------
## ============================================================================
cat("\n================ 3. IMvigor210: per-algorithm feature matrices ================\n")
dc <- readRDS(file.path(proj_dir, "data/IMvigor210/deconv_imvigor.rds"))
D  <- as.data.frame(data.table::fread(file.path(out_dir, "step4_analysis_table.csv")),
                    check.names = FALSE)

## Orientation verified by probe (not assumed):
##   dc$epic      : samples x cell types (348 x 8), rownames = sample IDs
##   dc$mcp       : cell types x samples (10 x 348)          -> transpose
##   dc$quantiseq : data.frame samples x cell types (348 x 12), first column "Sample"
ep <- dc$epic;                        stopifnot(nrow(ep) == 348L, ncol(ep) == 8L)
mc <- t(dc$mcp);                      stopifnot(nrow(mc) == 348L, ncol(mc) == 10L)
qs <- as.data.frame(dc$quantiseq); rownames(qs) <- qs$Sample; qs$Sample <- NULL
qs <- as.matrix(qs);                  stopifnot(nrow(qs) == 348L, ncol(qs) == 11L)

colnames(ep) <- paste0("e_",   colnames(ep))
colnames(mc) <- paste0("mcp_", colnames(mc))
colnames(qs) <- paste0("q_",   colnames(qs))

## Assemble ONE sample x feature matrix so every index is scored identically
Mall <- cbind(ep, mc, qs)
stopifnot(!any(duplicated(colnames(Mall))))
idx_row <- match(D$sample, rownames(Mall)); stopifnot(!anyNA(idx_row))
Mall <- Mall[idx_row, , drop = FALSE]
cat(sprintf("IMvigor210 combined feature matrix: %d samples x %d features\n", nrow(Mall), ncol(Mall)))

## Guard 4: every feature used by any index must exist here
for (k in names(IDX)) {
  x <- IDX[[k]]; if (is.null(x)) next
  miss <- base::setdiff(x$feats, colnames(Mall))
  if (length(miss)) stop(sprintf("%s: missing features in IMvigor210: %s",
                                 x$name, paste(miss, collapse = ", ")))
}
cat("Guard 4 passed: all index features present in IMvigor210\n")

## ============================================================================
## 4. Transport each index and evaluate --------------------------------------
## ============================================================================
cat("\n================ 4. Transport + evaluation in IMvigor210 ================\n")
S0 <- D[!is.na(D$resp), ]; S0$ic <- factor(S0$ic, levels = c("IC0", "IC1", "IC2+"))
cat(sprintf("Response-evaluable set: n = %d (CR/PR %d)\n", nrow(S0), sum(S0$resp)))
stopifnot(nrow(S0) == 298L, sum(S0$resp) == 68L)                # Guard 5
## Guard 5b: select the IC0 subset with which(); `S0$ic == "IC0"` returns NA for
## the single unknown IC level, and logical indexing with NA inserts NA rows
ic0_0 <- S0[which(S0$ic == "IC0"), ]
stopifnot(nrow(ic0_0) == 83L, sum(ic0_0$resp) == 13L)

score_cz <- function(M, feats, w) {                             # cohort-internal z (primary convention)
  X <- as.matrix(M[, feats, drop = FALSE])
  sdv <- apply(X, 2, stats::sd)
  if (any(!is.finite(sdv) | sdv <= 0))
    stop("Zero-variance feature in IMvigor210: ",
         paste(feats[!is.finite(sdv) | sdv <= 0], collapse = ", "))
  as.numeric(scale(X) %*% w)
}
## NB: assign through an explicit vector. `D["col"][idx] <- val` would NOT work:
## single-bracket indexing on a data.frame selects COLUMNS, so a numeric idx would
## silently create new columns instead of filling rows.
row_of <- match(rownames(Mall), D$sample)
v_main <- rep(NA_real_, nrow(D)); v_main[row_of] <- score_cz(Mall, W_MAIN, W_MAIN_W)
D$sirs_main_re <- v_main

for (k in names(IDX)) {
  x <- IDX[[k]]; if (is.null(x)) next
  v <- rep(NA_real_, nrow(D)); v[row_of] <- score_cz(Mall, x$feats, unname(x$weights))
  D[[paste0("idx_", k)]] <- v
}
## Guard 6: the recomputed main SIRS must equal the stored sirs_cz
o <- is.finite(D$sirs_cz) & is.finite(D$sirs_main_re)
cat(sprintf("Guard 6: recomputed main SIRS vs analysis_table sirs_cz: max |delta| = %.2e (n = %d)\n",
            max(abs(D$sirs_cz[o] - D$sirs_main_re[o])), sum(o)))
stopifnot(max(abs(D$sirs_cz[o] - D$sirs_main_re[o])) < 1e-6)

S <- D[!is.na(D$resp), ]; S$ic <- factor(S$ic, levels = c("IC0", "IC1", "IC2+"))
ic0 <- S[which(S$ic == "IC0"), ]
markers <- c("sirs_cz", paste0("idx_", intersect(names(IDX),
             c("MCP", "QS", "EPIC", "MCP_all", "QS_all", "EPIC_all"))))
markers <- markers[markers %in% names(S)]
m_lab <- c(sirs_cz = "Main SIRS (mixed)", idx_MCP = "MCP-counter (FDR-selected)",
           idx_QS = "quanTIseq (FDR-selected)", idx_EPIC = "EPIC (FDR-selected)",
           idx_MCP_all = "MCP-counter (all features)", idx_QS_all = "quanTIseq (all features)",
           idx_EPIC_all = "EPIC (all features)")
markers <- markers[!vapply(markers, function(v) length(unique(S[[v]])) < 2, TRUE)]  # drop constant scores

## 4.1 consistency with the main SIRS inside IMvigor210
corr_imvig <- do.call(rbind, lapply(markers[markers != "sirs_cz"], function(v) {
  ct <- suppressWarnings(stats::cor.test(D$sirs_cz, D[[v]], method = "spearman"))
  data.frame(marker = v, n = nrow(D), rho_vs_main = unname(ct$estimate), stringsAsFactors = FALSE)
}))
cat("\n-- IMvigor210: index vs main SIRS (cohort z) --\n"); print(corr_imvig, row.names = FALSE, digits = 4)

## 4.2 overall AUC
ov <- do.call(rbind, lapply(markers, function(v) {
  b <- boot_auc(S[[v]], S$resp); pr <- perm_auc(S[[v]], S$resp)
  data.frame(marker = m_lab[v], auc = b[["auc"]], ci_lo = b[["lo"]], ci_hi = b[["hi"]],
             perm_p = pr[["p_emp"]], n = nrow(S), n_resp = sum(S$resp), stringsAsFactors = FALSE)
}))
cat("\n-- overall AUC (all patients) --\n"); print(ov, row.names = FALSE, digits = 4)
data.table::fwrite(ov, file.path(out_dir, "step4c_imvig_auc_overall.csv"))

## 4.3 pre-specified stratification by PD-L1 IC level + interaction test
strat <- data.frame(); inter <- data.frame()
for (v in markers) {
  for (lv in levels(S$ic)) {
    s <- S[which(S$ic == lv), ]
    b <- boot_auc(s[[v]], s$resp); pr <- perm_auc(s[[v]], s$resp)
    strat <- rbind(strat, data.frame(marker = m_lab[v], ic_level = lv, auc = b[["auc"]],
      ci_lo = b[["lo"]], ci_hi = b[["hi"]], perm_p = pr[["p_emp"]],
      n = nrow(s), n_resp = sum(s$resp), stringsAsFactors = FALSE))
  }
  S$.mk <- S[[v]]
  mi <- glm(resp ~ scale(.mk) * ic01, data = S, family = binomial)
  ci <- summary(mi)$coef
  inter <- rbind(inter, data.frame(marker = m_lab[v], beta_interaction = ci[4, 1],
                                   p_interaction = ci[4, 4], aic = stats::AIC(mi),
                                   stringsAsFactors = FALSE))
}
cat("\n-- stratified AUC by PD-L1 IC level --\n"); print(strat, row.names = FALSE, digits = 4)
cat("\n-- interaction test: response ~ marker * IC0 --\n"); print(inter, row.names = FALSE, digits = 4)
data.table::fwrite(strat, file.path(out_dir, "step4c_imvig_ic_stratified_auc.csv"))

## 4.4 IC0 head-to-head (DeLong, correlated ROC)
dl <- data.frame()
for (v in markers[markers != "sirs_cz"]) {
  dp <- delong_p_paired(ic0[["sirs_cz"]], ic0[[v]], ic0$resp)
  dl <- rbind(dl, data.frame(comparison = sprintf("main SIRS vs %s", m_lab[v]),
    auc_left = auc_mw(ic0[["sirs_cz"]], ic0$resp), auc_right = auc_mw(ic0[[v]], ic0$resp),
    delta = dp[["d"]], n = dp[["n"]], delong_p = dp[["p"]], stringsAsFactors = FALSE))
  dp2 <- delong_p_paired(ic0[[v]], ic0$tmb, ic0$resp)
  dl <- rbind(dl, data.frame(comparison = sprintf("%s vs TMB", m_lab[v]),
    auc_left = auc_mw(ic0[[v]], ic0$resp), auc_right = auc_mw(ic0$tmb, ic0$resp),
    delta = dp2[["d"]], n = dp2[["n"]], delong_p = dp2[["p"]], stringsAsFactors = FALSE))
}
cat("\n-- IC0 head-to-head (DeLong) --\n"); print(dl, row.names = FALSE, digits = 4)
data.table::fwrite(dl, file.path(out_dir, "step4c_ic0_headtohead_delong.csv"))

## 4.5 complementarity of every index with TMB
comp <- do.call(rbind, lapply(markers, function(v) {
  do.call(rbind, lapply(levels(S$ic), function(lv) {
    s <- S[which(S$ic == lv), ]
    data.frame(marker = m_lab[v], ic_level = lv, n = nrow(s),
               auc_index = auc_mw(s[[v]], s$resp), auc_tmb = auc_mw(s$tmb, s$resp),
               stringsAsFactors = FALSE)
  }))
}))
comp$stronger <- ifelse(comp$auc_index > comp$auc_tmb, "index", "TMB")
cat("\n-- complementarity: index vs TMB within each PD-L1 stratum --\n")
print(comp, row.names = FALSE, digits = 4)
data.table::fwrite(comp, file.path(out_dir, "step4c_complementarity.csv"))

## 4.6 WHY the IC0 AUCs coincide: rank agreement within the IC0 stratum.
##     An identical AUC with DeLong p = 1 can only arise from an identical ordering,
##     so this table is what makes the coincidence interpretable rather than mysterious.
ic0_rank <- do.call(rbind, lapply(markers[markers != "sirs_cz"], function(v) {
  ok <- is.finite(ic0[["sirs_cz"]]) & is.finite(ic0[[v]])
  ct <- suppressWarnings(stats::cor.test(ic0[["sirs_cz"]][ok], ic0[[v]][ok], method = "spearman"))
  data.frame(marker = m_lab[v], n = sum(ok), rho_within_IC0 = unname(ct$estimate),
             auc_IC0 = auc_mw(ic0[[v]], ic0$resp), stringsAsFactors = FALSE)
}))
cat("\n-- rank agreement with the main SIRS WITHIN IC0 (explains identical AUCs) --\n")
print(ic0_rank, row.names = FALSE, digits = 4)

## 4.7 export the per-sample scores (reproducibility of every index)
score_out <- data.frame(sample = D$sample, resp = D$resp, ic = as.character(D$ic),
                        tmb = D$tmb, sirs_cz = D$sirs_cz)
for (k in names(IDX)) if (!is.null(D[[paste0("idx_", k)]]))
  score_out[[paste0("idx_", k)]] <- D[[paste0("idx_", k)]]
data.table::fwrite(score_out, file.path(out_dir, "step4c_index_scores_imvig.csv"))
cat(sprintf("\nPer-sample scores exported: step4c_index_scores_imvig.csv (%d rows, %d columns)\n",
            nrow(score_out), ncol(score_out)))

## ============================================================================
## 5. Figure ------------------------------------------------------------------
## ============================================================================
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  if (isTRUE(capabilities("cairo"))) options(ggplot2.device = grDevices::cairo_pdf)
  plt <- strat
  lv  <- unname(m_lab[markers])
  plt$marker <- factor(plt$marker, levels = lv)
  plt$grp <- factor(ifelse(grepl("^Main", lv), "Main SIRS (mixed)",
                    ifelse(grepl("^MCP", lv), "MCP-counter",
                    ifelse(grepl("^quanTIseq", lv), "quanTIseq", "EPIC")))[match(plt$marker, lv)],
                    levels = c("Main SIRS (mixed)", "MCP-counter", "quanTIseq", "EPIC"))
  plt$variant <- ifelse(grepl("all features", lv), "all features",
                 ifelse(grepl("FDR-selected", lv), "FDR-selected", "mixed"))[match(plt$marker, lv)]
  p <- ggplot(plt, aes(x = ic_level, y = auc, colour = variant, group = variant)) +
    geom_hline(yintercept = 0.5, linetype = 2, colour = "grey45") +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = .16) +
    facet_wrap(~ grp, nrow = 1) +
    scale_colour_manual(values = c(mixed = "#111111", `FDR-selected` = "#C0392B",
                                   `all features` = "#2471A3")) +
    labs(x = "PD-L1 immune-cell level", y = "AUC (objective response)", colour = NULL,
         title = "Algorithm sensitivity of the SIRS-response association",
         subtitle = "Bootstrap 95% CI; cohort-internal z standardisation") +
    theme_bw(base_size = 10.5) + theme(legend.position = "top")
  ggsave(file.path(fig_dir, "step4c_algorithm_ic_stratified_auc.pdf"), p,
         width = 8.4, height = 4.0, device = grDevices::cairo_pdf)
  cat("\nFigure written: figs/step4c_algorithm_ic_stratified_auc.pdf\n")
}

## ============================================================================
## 6. Summary -----------------------------------------------------------------
## ============================================================================
f <- function(x, d = 3) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")
sm <- c("# Step4c summary - single-algorithm sensitivity of SIRS", "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Purpose: rebuild SIRS from ONE deconvolution algorithm at a time under the",
  "identical construction rule, and re-test the IC0 result in IMvigor210.", "",
  sprintf("TCGA modelling set (identical to Step3): n = %d, events = %d.", nrow(dm), sum(dm$os_event)),
  sprintf("Full-pool re-derivation reproduces the Step3 weights exactly (max |delta| = %.1e).", dmax),
  "", "## 1. Indices constructed", "", "| index | n features | features |", "|---|---|---|")
for (i in seq_len(nrow(w_tab))) {
  sm <- c(sm, sprintf("| %s | %d | %s |", w_tab$index[i], w_tab$n_features[i], w_tab$features[i]))
}
sm <- c(sm, "", sprintf("Note: the EPIC-only feature set is identical across three independent definitions (pool-FDR / top-3 / full-pool sub-score): %s.",
  paste(sort(IDX$EPIC$feats), collapse = ", ")),
  sprintf("Degeneracy check: max zero fraction over all pool features = %.3f, so no feature is degenerate; a degeneracy-filtered EPIC variant is IDENTICAL to EPIC-only.",
          max(zf$zero_frac)), "",
  "## 2. Overall AUC in IMvigor210", "", "| index | AUC | 95% CI | permutation p |", "|---|---|---|---|")
for (i in seq_len(nrow(ov))) sm <- c(sm, sprintf("| %s | %s | %s-%s | %s |",
  ov$marker[i], f(ov$auc[i]), f(ov$ci_lo[i]), f(ov$ci_hi[i]), f(ov$perm_p[i], 4)))
sm <- c(sm, "", "## 3. Stratified AUC by PD-L1 IC level", "",
  "| index | IC0 | IC1 | IC2+ | interaction p |", "|---|---|---|---|---|")
for (v in unname(m_lab[markers])) {
  g <- function(lv) f(strat$auc[strat$marker == v & strat$ic_level == lv])
  sm <- c(sm, sprintf("| %s | %s | %s | %s | %s |", v, g("IC0"), g("IC1"), g("IC2+"),
                      f(inter$p_interaction[inter$marker == v], 4)))
}
sm <- c(sm, "", "## 3b. Within-IC0 rank agreement with the main SIRS", "",
  "| index | rho within IC0 | AUC in IC0 |", "|---|---|---|")
for (i in seq_len(nrow(ic0_rank)))
  sm <- c(sm, sprintf("| %s | %s | %s |", ic0_rank$marker[i],
                      f(ic0_rank$rho_within_IC0[i]), f(ic0_rank$auc_IC0[i])))
sm <- c(sm, "", "## 4. Consistency with the main SIRS", "",
  "| index | rho vs main (TCGA) | rho vs main (IMvigor210) |", "|---|---|---|")
for (k in names(IDX)) {
  x <- IDX[[k]]; if (is.null(x) || k == "main_mixed") next
  tc_r <- corr_tcga$rho_vs_main[corr_tcga$algorithm == k]
  im_r <- corr_imvig$rho_vs_main[corr_imvig$marker == paste0("idx_", k)]
  sm <- c(sm, sprintf("| %s | %s | %s |", x$name,
                      f(if (length(tc_r)) tc_r else NA_real_),
                      f(if (length(im_r)) im_r else NA_real_)))
}
sm <- c(sm, "", "## 5. Files", paste0("- step4_output/", c(
  "step4c_weights_by_algorithm.csv", "step4c_tcga_exposure_validity.csv",
  "step4c_imvig_auc_overall.csv", "step4c_imvig_ic_stratified_auc.csv",
  "step4c_ic0_headtohead_delong.csv", "step4c_complementarity.csv",
  "step4c_index_scores_imvig.csv", "figs/step4c_algorithm_ic_stratified_auc.pdf")))
writeLines(sm, file.path(out_dir, "12_Step4c_summary.md"))
cat("\nSummary written: step4_output/12_Step4c_summary.md\n")
cat("\n========== Step4c complete ==========\n")
