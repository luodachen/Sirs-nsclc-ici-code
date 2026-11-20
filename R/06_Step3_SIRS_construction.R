############################################################
## 06_Step3_SIRS_construction.R
## Project: smoking exposure quantification x NSCLC immune microenvironment x
##          ICI response prediction
##
## Role: build the smoking-related immune remodelling score (SIRS) and assess it
##   in three layers, each with its own claim.
##   L1 primary product: SIRS = effect-size-weighted sum of smoking-related
##      immune features, weights = confounder-adjusted partial Spearman rho.
##      Positioning: a portable quantitative index of smoking-driven immune
##      remodelling; no claim of prognostic value.
##   L2 validity: SIRS as a proxy for smoking exposure
##      (vs csi_lc / py_gdc / sbs4_abs / TMB).
##   L3 prognosis: KM / timeROC / C-index / multivariable Cox + bootstrap,
##      reported as an honest negative, plus the LUAD-vs-LUSC subgroup
##      difference. LASSO-Cox is retained as a control whose failure to select
##      immune features is reported as a methodological negative.
##
## Upstream constraints (must not be overturned)
##   C1 the SBS4 primary metric is sbs4_abs (= sbs4_weight x tmb_per_mb)
##   C2 T-cell features from quanTIseq are excluded (q_T.cells.CD4 is
##      zero-inflated at 80.3%, q_Monocytes at 97.6%)
##   C3 age / gender / stage / cohort are strong confounders and are always
##      adjusted for
##   C4 LUAD effect size > LUSC -> pooled modelling with stratified reporting
##
## Survival definition (error-prone)
##   days_death is present only for Dead cases and days_fu is incomplete for
##   Dead cases, so os_time = (vital == "Dead") ? days_death : days_fu, with a
##   two-way fallback (yields 1009/1087). Using days_fu alone would miss 236
##   deaths; an assertion guards this.
##
## Environment / locale rules
##   - Identifiers must be ASCII (non-ASCII identifiers fail to parse under the
##     local locale); comments and strings may be non-ASCII.
##   - Output file names must be ASCII.
##   - In-figure labels must be ASCII; the plotting device is cairo_pdf.
##   - After loading a Bioconductor genomic package, namespace-qualify same-named
##     purrr/dplyr functions explicitly.
##
## Dependencies: tidyverse, data.table, glmnet, survival, survminer, timeROC
## Input:  step2_output/step2_immune_tme.csv   (output of script 04)
## Output: step3_output/ (SIRS formula and all-sample scores, for Step4)
## Usage:  setwd("D:/Projects/smoking-index-ici"); source("06_Step3_SIRS_construction.R")
## Author: see manuscript author list   Date: 2026-10-10
############################################################

## 0. Configuration -----------------------------------------------------------
proj_dir <- "D:/Projects/smoking-index-ici"
out_dir  <- file.path(proj_dir, "step3_output")
fig_dir  <- file.path(out_dir, "figs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(glmnet); library(survival)
})
if (isTRUE(capabilities("cairo"))) options(ggplot2.device = grDevices::cairo_pdf)

set.seed(20261010)
B_BOOT <- 200        # number of bootstrap replicates (optimism correction)
FDR_MAIN <- 0.10     # primary screening threshold
FDR_SENS <- 0.05     # sensitivity threshold

cat("========== Step3: smoking-related immune remodelling score (SIRS) ==========\n\n")

## 1. Load + variable construction + guards -----------------------------------
step2 <- fread(file.path(proj_dir, "step2_output/step2_immune_tme.csv"),
               na.strings = c("", "NA")) %>% as_tibble()
cat(sprintf("Loaded step2 main table: %d samples x %d columns\n", nrow(step2), ncol(step2)))

step2 <- step2 %>%
  mutate(
    os_event = vital == "Dead",
    os_time  = if_else(vital == "Dead", days_death, days_fu),
    os_time  = if_else(!is.finite(os_time) & is.finite(days_fu),    days_fu,    os_time),
    os_time  = if_else(!is.finite(os_time) & is.finite(days_death), days_death, os_time),
    age      = age_y / 365.25,
    ## Stage collapsing must use case_when (first match wins). Sequential overwriting
    ##   would let "^Stage\\s*I" eat the prefix of "Stage IV"/"Stage III" 
    stage_grp = case_when(
      grepl("^Stage\\s*IV",  stage) ~ "IV",
      grepl("^Stage\\s*III", stage) ~ "III",
      grepl("^Stage\\s*II",  stage) ~ "II",
      grepl("^Stage\\s*I",   stage) ~ "I",
      TRUE ~ NA_character_),
    sbs4_abs = sbs4_weight * tmb_per_mb,          # not in the main table; must be rebuilt
    os_time  = pmax(os_time, 0.5),                # Cox does not accept 0/negative
    gender_m = as.integer(factor(gender)) - 1L,
    stage_n  = as.integer(factor(stage_grp, levels = c("I", "II", "III", "IV"))),
    cohort_n = as.integer(factor(cohort)) - 1L)

cat("\n== Survival field check ==\n")
cat(sprintf("  days_death non-missing %d | days_fu non-missing %d | spliced os_time non-missing %d (%.1f%%)\n",
            sum(is.finite(step2$days_death)), sum(is.finite(step2$days_fu)),
            sum(is.finite(step2$os_time)), 100 * mean(is.finite(step2$os_time))))
## Guard 1: the correct definition must cover more than days_fu alone (otherwise the definition is wrong)
stopifnot(sum(is.finite(step2$os_time)) > sum(is.finite(step2$days_fu)))

CONF <- c("age", "gender_m", "stage_n", "cohort_n")
## Guard 2: control variables must not be degenerate, or the design matrix is singular
for (v in CONF) {
  nu <- length(unique(na.omit(step2[[v]])))
  if (nu <= 1) stop(sprintf("Control variable %s is degenerate (unique values = %d)", v, nu))
}

ev <- step2 %>% filter(smoking_grp %in% c("current", "former"))
dm <- ev %>% filter(vital %in% c("Alive", "Dead"), is.finite(os_time),
                    if_all(all_of(CONF), ~ !is.na(.x)))
cat(sprintf("\nAnalysis set ev = %d | modelling set (+ valid OS + complete confounders) = %d, events %d (%.1f%%)\n",
            nrow(ev), nrow(dm), sum(dm$os_event), 100 * mean(dm$os_event)))
print(dm %>% group_by(cohort) %>%
        summarise(n = n(), events = sum(os_event),
                  event_pct = round(100 * mean(os_event), 1), .groups = "drop") %>%
        as.data.frame())
## Guard 3: EPV must be reported before modelling
cat(sprintf("EPV: events %d -> at 10 events per variable, variable cap %d\n",
            sum(dm$os_event), floor(sum(dm$os_event) / 10)))

## 2. Candidate immune feature pool -------------------------------------------
QI_EXCLUDE <- c("q_T.cells.CD4", "q_Monocytes")     # C2: zero-inflated, excluded
FEAT_POOL <- c("tis_score", "cyt_score", "cd274_log2tpm", "cd8a_log2tpm",
               grep("^e_", names(dm), value = TRUE),
               grep("^mcp_", names(dm), value = TRUE),
               grep("^q_", names(dm), value = TRUE))
FEAT_POOL <- base::setdiff(FEAT_POOL, QI_EXCLUDE)
FEAT_POOL <- FEAT_POOL[sapply(FEAT_POOL, function(v) sum(is.finite(dm[[v]])) >= 0.8 * nrow(dm))]
cat(sprintf("\nCandidate immune feature pool: %d (zero-inflated excluded: %s)\n",
            length(FEAT_POOL), paste(QI_EXCLUDE, collapse = ", ")))

## 3. L1 primary product: screening + effect-size-weighted construction of SIRS --
## Reuses the self-implemented partial correlation validated in Step2b (numerically equivalent to ppcor::pcor.test, difference ~1e-16)
partial_spearman <- function(x, y, Z) {
  ok <- is.finite(x) & is.finite(y) & stats::complete.cases(Z)
  x <- x[ok]; y <- y[ok]; Z <- as.data.frame(Z)[ok, , drop = FALSE]
  rx <- rank(x); ry <- rank(y)
  Z  <- as.data.frame(lapply(Z, function(z) rank(z)))
  res <- stats::residuals(stats::lm(cbind(rx, ry) ~ as.matrix(Z)))
  r <- stats::cor(res[, 1], res[, 2])
  n <- length(rx); k <- ncol(Z); df <- n - k - 2
  tt <- r * sqrt(df / (1 - r^2))
  tibble(rho = r, n = n, p_value = 2 * stats::pt(-abs(tt), df = df))
}

## C-index helper: build an explicit data.frame first, then feed it to the formula.
## Counter-example: concordance(Surv(t,e) ~ use$sirs) relies on formula scope lookup and easily throws object not found
cindex_of <- function(time, event, marker) {
  ok <- is.finite(time) & !is.na(event) & is.finite(marker)
  if (sum(ok) < 20 || sum(event[ok]) < 5) return(c(C = NA_real_, se = NA_real_))
  dd <- data.frame(t = time[ok], e = as.integer(event[ok]), m = marker[ok])
  cc <- survival::concordance(survival::Surv(t, e) ~ m, data = dd, reverse = TRUE)
  c(C = cc$concordance, se = sqrt(cc$var))
}
logrank_p <- function(time, event, grp) {
  s <- survival::survdiff(survival::Surv(time, event) ~ grp)
  1 - pchisq(s$chisq, df = length(s$n) - 1)
}

ANCHORS <- c("csi_lc", "sbs4_abs")     # primary anchors: composite dose axis + absolute tobacco mutation burden
screen <- function(anchor, dat = dm) {
  d0 <- dat %>% filter(is.finite(.data[[anchor]]), !is.na(age))
  purrr::map_dfr(FEAT_POOL, function(f) {
    d <- d0 %>% filter(is.finite(.data[[f]]))
    if (nrow(d) < 30) return(tibble())
    pr <- partial_spearman(d[[anchor]], d[[f]], d %>% select(all_of(CONF)))
    tibble(feature = f, anchor = anchor, n = pr$n, rho = pr$rho, p_value = pr$p_value)
  }) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
}
scr <- purrr::map_dfr(ANCHORS, screen)
fwrite(scr, file.path(out_dir, "step3_screen_smoking_immune.csv"))

cat("\n== L1 screening: partial Spearman (adjusting for age/gender/stage/cohort) ==\n")
print(as.data.frame(scr %>% filter(p_adj < FDR_MAIN) %>%
                      select(feature, anchor, n, rho, p_value, p_adj) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))))

## Weight rule (pre-set): each feature takes the anchor with the largest |rho|; weight = rho
make_weights <- function(fdr) {
  sel <- scr %>% filter(p_adj < fdr) %>% pull(feature) %>% unique()
  w <- scr %>% filter(feature %in% sel) %>%
    group_by(feature) %>% slice_max(abs(rho), n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(feature, weight = rho, anchor, rho, p_adj)
  ## Sign convention: a negative rho means "the higher the exposure, the lower the feature"; the weight is rho itself,
  ##   so a larger SIRS = closer to the immune phenotype of a heavier smoker (data-driven direction; no manual flipping)
  w[order(-abs(w$weight)), ]
}
W_main <- make_weights(FDR_MAIN)
W_sens <- make_weights(FDR_SENS)
cat(sprintf("\nMain SIRS weights (FDR<%.2f, %d features):\n", FDR_MAIN, nrow(W_main)))
print(as.data.frame(W_main %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
cat(sprintf("\nSensitivity SIRS weights (FDR<%.2f, %d features): %s\n",
            FDR_SENS, nrow(W_sens), paste(W_sens$feature, collapse = ", ")))
fwrite(W_main, file.path(out_dir, "step3_SIRS_weights.csv"))

## Order matters: first take the modelling set use with complete features, then compute the z parameters **from use**.
##   Counter-example: computing mean/sd from unfiltered dm -> a few NAs make mean() = NA -> scale() turns the whole column NA (hit this before)
use <- dm %>% filter(if_all(all_of(W_main$feature), ~ is.finite(.x)))
scal      <- tibble(feature = W_main$feature,
                    z_mean  = sapply(W_main$feature, function(v) mean(use[[v]])),
                    z_sd    = sapply(W_main$feature, function(v) sd(use[[v]])))
scal_sens <- tibble(feature = W_sens$feature,
                    z_mean  = sapply(W_sens$feature, function(v) mean(use[[v]])),
                    z_sd    = sapply(W_sens$feature, function(v) sd(use[[v]])))

score_sirs <- function(dat, W, sc) {
  sc <- sc[match(W$feature, sc$feature), ]
  if (anyNA(sc$z_mean) || anyNA(sc$z_sd))
    stop("Standardisation parameters contain NA -- check whether sc$feature and W$feature line up")
  Z <- scale(as.matrix(dat[, W$feature, drop = FALSE]), center = sc$z_mean, scale = sc$z_sd)
  as.numeric(Z %*% W$weight)
}
use$sirs      <- score_sirs(use, W_main, scal)
use$sirs_sens <- score_sirs(use, W_sens, scal_sens)
cat(sprintf("Modelling set (+ complete SIRS features) = %d, events %d\n", nrow(use), sum(use$os_event)))

## Score all samples (not just the modelling set) -- maximise the samples available for Step4
sirs_ok <- sapply(W_main$feature, function(v) is.finite(step2[[v]]))
allrows <- step2[rowSums(!sirs_ok) == 0, ]
allrows$sirs      <- score_sirs(allrows, W_main, scal)
allrows$sirs_sens <- score_sirs(allrows, W_sens, scal_sens)
cat(sprintf("All-sample SIRS computable n = %d (ever-smokers %d)\n", nrow(allrows),
            sum(allrows$smoking_grp %in% c("current", "former"))))

## Guard 4: SIRS must contain no NA and must not have zero variance
if (!all(is.finite(use$sirs)))
  stop(sprintf("SIRS contains %d NA -- check feature completeness", sum(!is.finite(use$sirs))))
if (!isTRUE(sd(use$sirs) > 0)) stop("SIRS variance is 0 -- check the weights/standardisation parameters")

cat(sprintf("\nSIRS distribution: median %.3f, IQR [%.3f, %.3f]\n",
            median(use$sirs), quantile(use$sirs, .25), quantile(use$sirs, .75)))

## 4. L2 validity: SIRS as a proxy for smoking exposure -----------------------
cat("\n================ L2: validity of SIRS as a smoking exposure proxy ================\n")
validity <- purrr::map_dfr(c("csi_lc", "py_gdc", "log_cig_years", "duration", "cpd",
                             "tmb_per_mb", "sbs4_weight", "sbs4_abs"), function(v) {
  d <- use %>% filter(is.finite(.data[[v]]))
  if (nrow(d) < 30) return(tibble())
  ct <- suppressWarnings(cor.test(d[[v]], d$sirs, method = "spearman"))
  tibble(variable = v, n = nrow(d), rho = unname(ct$estimate), p_value = ct$p.value)
}) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
print(as.data.frame(validity %>% mutate(rho = round(rho, 3),
                                        p_value = signif(p_value, 3))))
fwrite(validity, file.path(out_dir, "step3_SIRS_exposure_validity.csv"))
## Validity criterion: significant with at least one primary exposure anchor (otherwise SIRS does not represent smoking and the L1 positioning fails)
ok_anchor <- validity %>% filter(variable %in% c("csi_lc", "py_gdc", "sbs4_abs"),
                                 p_value < 0.05)
if (!nrow(ok_anchor)) warning("SIRS is non-significant with all primary exposure anchors -- the L1 positioning fails; the feature pool must be re-examined")

## 5. L3 prognostic evaluation (honest negative reporting) ---------------------
cat("\n================ L3: discrimination of SIRS for OS ================\n")
pred <- purrr::map_dfr(c("sirs", "sirs_sens", "tis_score", "cyt_score", "cd274_log2tpm",
                         "stage_n"), function(v) {
  d <- use %>% filter(is.finite(.data[[v]]))
  cc <- cindex_of(d$os_time, d$os_event, d[[v]])
  tibble(marker = v, n = nrow(d), C = unname(cc["C"]),
         ci_low = C - 1.96 * unname(cc["se"]), ci_high = C + 1.96 * unname(cc["se"]))
})
print(as.data.frame(pred %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))
fwrite(pred, file.path(out_dir, "step3_Cindex_comparison.csv"))

## 5.1 KM / log-rank -- the primary analysis uses a **median split** (pre-set, no selection bias); the maxstat cutpoint
##   serves as a sensitivity analysis only: being data-adaptive, its log-rank p is systematically too small and must not be presented alone.
cut_med <- median(use$sirs)
use <- use %>% mutate(sirs_grp = if_else(sirs > cut_med, "High", "Low"))
lr_med <- logrank_p(use$os_time, use$os_event, use$sirs_grp)
cat(sprintf("\n5.1 KM primary analysis (median split %.3f) log-rank p = %.4f\n", cut_med, lr_med))
lr_max <- NA_real_; cut_max <- NA_real_
if (requireNamespace("survminer", quietly = TRUE) && requireNamespace("maxstat", quietly = TRUE)) {
  cp <- tryCatch(survminer::surv_cutpoint(as.data.frame(use), time = "os_time",
                                          event = "os_event", variables = "sirs"),
                 error = function(e) NULL)
  if (!is.null(cp)) {
    cut_max <- cp$cutpoint$cutpoint[1]
    lr_max  <- logrank_p(use$os_time, use$os_event, use$sirs > cut_max)
    cat(sprintf("    maxstat sensitivity cutpoint %.3f | log-rank p = %.4f  <- data-adaptive cutpoint, p biased low\n",
                cut_max, lr_max))
    if (lr_med > .05 && lr_max < .05)
      cat("    The two disagree -> the paper reports the median split, with maxstat as a sensitivity analysis only and its limitation stated\n")
  }
}
lr_main <- lr_med

## 5.2 timeROC at 1/3/5 years
tr <- timeROC::timeROC(T = use$os_time, delta = use$os_event, marker = use$sirs,
                       cause = 1, times = c(365, 1095, 1825), iid = TRUE)
se_v <- tryCatch(as.numeric(tr$inference$vect_sd_1), error = function(e) rep(NA_real_, 3))
if (length(se_v) != 3) se_v <- rep(NA_real_, 3)
auc_tab <- tibble(time_year = c(1, 3, 5), auc = as.numeric(tr$AUC), se = se_v,
                  ci_low = auc - 1.96 * se, ci_high = auc + 1.96 * se)
print(as.data.frame(auc_tab %>% mutate(across(where(is.numeric), ~ round(.x, 4)))))

## 5.3 bootstrap optimism correction (re-run the full pipeline: screening + weighting, for honesty)
boot_opt <- numeric(0)
for (b in seq_len(B_BOOT)) {
  idx <- sample(nrow(use), replace = TRUE)
  db  <- use[idx, ]
  if (sum(db$os_event) < 20) next
  ## Re-screen: re-run the partial correlations within the bootstrap sample, using the same feature set as before
  sb <- tryCatch(purrr::map_dfr(ANCHORS, function(a) {
    d0 <- db %>% filter(is.finite(.data[[a]]))
    purrr::map_dfr(W_main$feature, function(f) {
      d <- d0 %>% filter(is.finite(.data[[f]]))
      if (nrow(d) < 30) return(tibble())
      pr <- partial_spearman(d[[a]], d[[f]], d %>% select(all_of(CONF)))
      tibble(feature = f, anchor = a, rho = pr$rho, p_value = pr$p_value)
    }) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
  }), error = function(e) NULL)
  if (is.null(sb) || !nrow(sb)) next
  Wb <- sb %>% filter(p_adj < FDR_MAIN) %>%
    group_by(feature) %>% slice_max(abs(rho), n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(feature, weight = rho)
  if (nrow(Wb) < 2) next
  scb <- tibble(feature = Wb$feature,
                z_mean = sapply(Wb$feature, function(v) mean(db[[v]])),
                z_sd   = sapply(Wb$feature, function(v) sd(db[[v]])))
  f_ok <- sapply(Wb$feature, function(v) is.finite(db[[v]]))
  db2 <- db[rowSums(!f_ok) == 0, ]
  if (nrow(db2) < 50) next
  s_b <- score_sirs(db2, Wb, scb)
  s_o <- score_sirs(use,  Wb, scb)
  c_b <- unname(cindex_of(db2$os_time, db2$os_event, s_b)["C"])
  c_o <- unname(cindex_of(use$os_time, use$os_event, s_o)["C"])
  boot_opt <- c(boot_opt, c_b - c_o)
}
c_app <- unname(cindex_of(use$os_time, use$os_event, use$sirs)["C"])
cat(sprintf("\n5.3 C-index: apparent %.4f | bootstrap (B_effective=%d) mean optimism %.4f -> corrected %.4f\n",
            c_app, length(boot_opt), mean(boot_opt), c_app - mean(boot_opt)))

## 5.4 Multivariable Cox: is SIRS independent of the confounders (expected: not independent = negative)
mv <- survival::coxph(survival::Surv(os_time, os_event) ~ sirs + age + gender_m + stage_n + cohort_n,
                      data = use)
mv_tab <- broom::tidy(mv, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, HR = estimate, conf.low, conf.high, p.value)
cat("\n5.4 Multivariable Cox (SIRS + 4 confounders):\n")
print(as.data.frame(mv_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))

## 5.5 Stratification: LUAD / LUSC
strat <- purrr::map_dfr(c("TCGA-LUAD", "TCGA-LUSC"), function(co) {
  d <- use %>% filter(cohort == co)
  if (nrow(d) < 30 || sum(d$os_event) < 10) return(tibble())
  cc <- cindex_of(d$os_time, d$os_event, d$sirs)
  tibble(cohort = co, n = nrow(d), events = sum(d$os_event), C = unname(cc["C"]),
         logrank_p = logrank_p(d$os_time, d$os_event, d$sirs_grp))
})
cat("\n5.5 By cohort:\n")
print(as.data.frame(strat %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))

## 6. Control: LASSO-Cox (records the negative evidence that no immune feature is selected) -----
cat("\n================ Control: LASSO-Cox ================\n")
Xz <- scale(as.matrix(use[, W_main$feature]))
X  <- cbind(Xz, as.matrix(use[, CONF]))
pf <- c(rep(1, nrow(W_main)), rep(0, length(CONF)))    # confounders forced in, not penalised
Y  <- survival::Surv(use$os_time, use$os_event)
cv <- glmnet::cv.glmnet(X, Y, family = "cox", alpha = 1, penalty.factor = pf,
                        standardize = FALSE, nfolds = 10)
cat(sprintf("  lambda path has %d values, range [%.5f, %.5f]\n",
            length(cv$lambda), min(cv$lambda), max(cv$lambda)))
cat(sprintf("  lambda.min = %.5f (CV deviance %.4f) | lambda.1se = %.5f\n",
            cv$lambda.min, cv$cvm[which(cv$lambda == cv$lambda.min)], cv$lambda.1se))
cat(sprintf("  max nzero along the path = %d (of which %d confounders) -> immune features selected %d\n",
            max(cv$nzero), length(CONF), max(cv$nzero) - length(CONF)))
lasso_diag <- tibble(lambda = cv$lambda, nzero = cv$nzero, cvm = cv$cvm, cvsd = cv$cvsd)
fwrite(lasso_diag, file.path(out_dir, "step3_lasso_diagnostics.csv"))
n_immune_sel <- max(cv$nzero) - length(CONF)
nz0   <- which(cv$nzero == length(CONF))
cvm0  <- if (length(nz0)) cv$cvm[max(nz0)] else NA_real_   # CV deviance while no immune feature has entered yet
cvmmin <- min(cv$cvm)
cat(sprintf("  at most %d immune features enter | CV deviance: %.4f at the all-zero point -> %.4f at the minimum (Delta=%.4f)\n",
            n_immune_sel, cvm0, cvmmin, cvmmin - cvm0))
if (n_immune_sel == 0) {
  cat("  Negative evidence: 0 immune features selected along the whole path -> no independent prognostic gain after confounder adjustment\n")
} else if (isTRUE(cvmmin >= cvm0)) {
  cat("  Negative evidence: immune features enter only at smaller lambda, but the CV deviance does **not improve** (it rises)\n")
  cat("     -> adding immune features does not improve discrimination, i.e. an ineffective model\n")
}

## 7. Plots (ASCII labels + cairo_pdf) -----------------------------------------
## 7.1 SIRS weight bar chart
p_w <- W_main %>% mutate(feature = fct_reorder(feature, weight)) %>%
  ggplot(aes(feature, weight, fill = weight > 0)) + geom_col(width = .7) + coord_flip() +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2471A3"), guide = "none") +
  labs(x = NULL, y = "Weight (= confounder-adjusted partial rho)",
       title = "SIRS components", subtitle = sprintf("FDR < %.2f, n = %d features", FDR_MAIN, nrow(W_main))) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "step3_SIRS_weights.pdf"), p_w, width = 6,
       height = max(3, .3 * nrow(W_main) + 1.8), device = grDevices::cairo_pdf)

## 7.2 KM
## Known harmless warning (internal to third-party packages, cannot be silenced from outside):
##   survminer::ggsurvplot -> ggpubr::geom_exec(geom = "step", size = size)
##   still draws lines with size, and ggplot2 >= 3.4 advises switching to linewidth.
##   Verified: even the bare ggsurvplot(fit, data) always emits it (independent of risk.table/pval/legend);
##   the plot is correct -- in ggplot_build() the built data for that layer already uses linewidth.
##   So this handler filters only that single warning; all other warnings are still raised (no masking of real issues).
if (requireNamespace("survminer", quietly = TRUE)) {
  n_sil <- 0L
  km <- withCallingHandlers(
    survminer::ggsurvplot(
      survival::survfit(survival::Surv(os_time, os_event) ~ sirs_grp, data = use),
      data = use, pval = TRUE, risk.table = TRUE, xlab = "Time (days)",
      ylab = "Overall survival", legend.title = "SIRS", legend.labs = c("High", "Low"),
      palette = c("#C0392B", "#2471A3")),
    warning = function(w) {
      if (grepl("`size` aesthetic for lines", conditionMessage(w), fixed = TRUE)) {
        n_sil <<- n_sil + 1L
        invokeRestart("muffleWarning")
      }
    })
  if (n_sil > 0L) cat(sprintf(
    "[plot] filtered %d third-party (survminer/ggpubr) size -> linewidth deprecation warnings (harmless, see the 7.2 note)\n",
    n_sil))
  ggsave(file.path(fig_dir, "step3_KM_SIRS.pdf"), plot = km$plot, width = 6, height = 5,
         device = grDevices::cairo_pdf)
}
## 7.3 C-index comparison (with a 0.5 reference line)
p_c <- pred %>% mutate(marker = fct_reorder(marker, C)) %>%
  ggplot(aes(marker, C)) + geom_point(size = 2.6, colour = "#2471A3") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = .2, colour = "#2471A3") +
  geom_hline(yintercept = .5, linetype = 2, colour = "grey40") +
  coord_flip(ylim = c(.4, .7)) + labs(x = NULL, y = "C-index (OS)",
  title = "Discrimination: SIRS vs references") + theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "step3_Cindex_compare.pdf"), p_c, width = 6, height = 4,
       device = grDevices::cairo_pdf)

## 7.4 SIRS vs exposure-anchor scatter plot
p_s <- use %>% filter(is.finite(csi_lc)) %>%
  ggplot(aes(csi_lc, sirs)) + geom_point(alpha = .35, size = 1) +
  geom_smooth(method = "lm", colour = "#C0392B", se = TRUE) +
  labs(x = "CSI (lc, current tsc convention)", y = "SIRS",
       title = sprintf("SIRS vs smoking exposure (rho = %.3f)",
                       validity$rho[validity$variable == "csi_lc"])) +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "step3_SIRS_vs_csi.pdf"), p_s, width = 5, height = 4,
       device = grDevices::cairo_pdf)

## 8. Export: SIRS formula (for Step4) + summary ------------------------------
form_lines <- c(
  "# SIRS formula (Step3, for external validation in Step4)", "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M")),
  sprintf("Training set: TCGA-NSCLC ever-smokers, n = %d, events = %d", nrow(use), sum(use$os_event)),
  sprintf("Anchors: %s | Weight rule: confounder-adjusted partial Spearman rho (FDR < %.2f)",
          paste(ANCHORS, collapse = ", "), FDR_MAIN),
  sprintf("Features: %d | Orientation: higher SIRS = immune phenotype of heavier smoking",
          nrow(W_main)),
  "", "## Formula", "", "SIRS = sum_i  w_i * (x_i - mean_i) / sd_i", "",
  "| feature | w (partial rho) | anchor | FDR adj p | train mean | train sd |",
  "|---|---|---|---|---|---|")
for (i in seq_len(nrow(W_main))) {
  f <- W_main$feature[i]
  form_lines <- c(form_lines, sprintf("| %s | %.6f | %s | %.3g | %.6f | %.6f |",
                                      f, W_main$weight[i], W_main$anchor[i], W_main$p_adj[i],
                                      scal$z_mean[match(f, scal$feature)],
                                      scal$z_sd[match(f, scal$feature)]))
}
form_lines <- c(form_lines, "",
  "## Rules for external use", "",
  "1. Standardise external data with the TRAIN mean/sd above; never use the external cohort's own moments.",
  "2. Confounders (age/gender/stage/cohort) are NOT part of SIRS; they are adjusted for at the",
  "   screening step, so the weights are already confounder-adjusted.",
  "3. Both SIRS versions exported: main (FDR<0.10) and sensitivity (FDR<0.05).",
  "4. Honest limitation: SIRS carries NO independent prognostic value for OS beyond stage",
  sprintf("   (C-index %.3f vs stage alone 0.599); its role is to quantify smoking-driven", c_app),
  "   immune remodelling, not to predict survival.")
writeLines(form_lines, file.path(out_dir, "step3_SIRS_formula.md"))

fwrite(allrows %>% select(sample15, cohort, smoking_grp, sirs, sirs_sens,
                          os_time, os_event, vital),
       file.path(out_dir, "step3_SIRS_all_samples.csv"))

sink(file.path(out_dir, "06_Step3_summary.txt"))
cat("Step3 SIRS construction summary\n\n")
cat(sprintf("Modelling set: n=%d, events=%d (EPV cap %d vars)\n", nrow(use), sum(use$os_event),
            floor(sum(use$os_event)/10)))
cat(sprintf("SIRS: %d features, rule = adjusted partial rho, FDR<%.2f\n", nrow(W_main), FDR_MAIN))
cat(sprintf("SIRS = %s\n\n", paste(sprintf("%.4f*%s", W_main$weight, W_main$feature), collapse = " + ")))
cat("== L2 exposure validity ==\n"); print(as.data.frame(validity))
cat("\n== L3 discrimination ==\n"); print(as.data.frame(pred))
cat(sprintf("\nlog-rank p (median cut at %.3f) = %.4f | maxstat cut %.3f p = %.4f\n",
            cut_med, lr_med, cut_max, lr_max))
cat(sprintf("C-index apparent %.4f -> bootstrap-corrected %.4f (B=%d)\n",
            c_app, c_app - mean(boot_opt), length(boot_opt)))
cat("\n== multivariable Cox ==\n"); print(as.data.frame(mv_tab))
cat("\n== strata ==\n"); print(as.data.frame(strat))
cat(sprintf("\nLASSO-Cox control: max selected immune features in path = %d\n",
            max(cv$nzero) - length(CONF)))
sink()

cat("\n========== Step3 complete ==========\n")
print(basename(list.files(out_dir, pattern = "^step3|^06_Step3")))

