#!/usr/bin/env Rscript
# Step 3: Find optimal BAI percentile cutoff via Leave-One-Cohort-Out CV (LOCO-CV).
#
# In each fold one cohort is held out. The remaining cohorts are used to sweep
# candidate percentile cutoffs and select the one whose binarised-BAI HR has the
# smallest p-value. That cutoff is then applied to the held-out cohort.
# The cross-validated HR is obtained by pooling predictions across all
# held-out cohorts and fitting the chosen model on the pooled data.
#
# Model options (--model=):
#   finegray  Fine-Gray competing-risk model, death as competing event (default)
#   coxph     Cause-specific Cox PH model (death treated as censored)
#
# Usage:
#   Rscript step3_find_optimal_cutoff.R                   # noAPOE4, Fine-Gray
#   Rscript step3_find_optimal_cutoff.R --apoe4           # withAPOE4
#   Rscript step3_find_optimal_cutoff.R --model=coxph     # cause-specific Cox PH
#   Rscript step3_find_optimal_cutoff.R --cutoffs=50:80   # custom range (default: 50,55,…,80)
#   Rscript step3_find_optimal_cutoff.R --horizons=5,7    # horizons (years) for Se/PPV/AUC (default: 5)
#   Rscript step3_find_optimal_cutoff.R --n-boot=1000     # bootstrap replicates for Se/PPV/AUC CIs
#   Rscript step3_find_optimal_cutoff.R --n-boot-select=1000 --cores=8
#                                       # bootstrap replicates of the cutoff selection (0 = skip)
#
#   --suffix= selects the step-2 output directly, overriding --apoe4. Use it for the
#   APOE4-matching comparison, where both arms share one reference pool and one set
#   of individuals (see step2 --without-apoe4-same-dataset):
#     Rscript step3_find_optimal_cutoff.R --suffix=_noAPOE4_sameDataset
#     Rscript step3_find_optimal_cutoff.R --apoe4

library(survival)

# ---- Parse CLI arguments ------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

apoe4      <- '--apoe4' %in% args

model_arg  <- grep('^--model=', args, value = TRUE)
model_type <- if (length(model_arg) > 0) sub('--model=', '', model_arg[1]) else 'finegray'
if (!model_type %in% c('finegray', 'coxph'))
  stop('--model must be "finegray" or "coxph"')

cutoff_arg <- grep('^--cutoffs=', args, value = TRUE)
if (length(cutoff_arg) > 0) {
  range_str <- sub('--cutoffs=', '', cutoff_arg[1])
  parts     <- as.integer(strsplit(range_str, ':')[[1]])
  cutoffs   <- seq(parts[1], parts[2])
} else {
  cutoffs <- seq(50, 95, 5)
}

horizon_arg <- grep('^--horizons=', args, value = TRUE)
horizons    <- if (length(horizon_arg) > 0)
  as.numeric(strsplit(sub('--horizons=', '', horizon_arg[1]), ',')[[1]]) else 5

nboot_arg <- grep('^--n-boot=', args, value = TRUE)
n_boot    <- if (length(nboot_arg) > 0) as.integer(sub('--n-boot=', '', nboot_arg[1])) else 1000

nboot_sel_arg <- grep('^--n-boot-select=', args, value = TRUE)
n_boot_select <- if (length(nboot_sel_arg) > 0) as.integer(sub('--n-boot-select=', '', nboot_sel_arg[1])) else 1000

cores_arg <- grep('^--cores=', args, value = TRUE)
n_cores   <- if (length(cores_arg) > 0) as.integer(sub('--cores=', '', cores_arg[1])) else 8

suffix_arg <- grep('^--suffix=', args, value = TRUE)
if (length(suffix_arg) > 0) {
  suffix <- sub('--suffix=', '', suffix_arg[1])
  if (!startsWith(suffix, '_')) suffix <- paste0('_', suffix)
  if (apoe4)
    cat('NOTE: --suffix= overrides --apoe4; using', suffix, '\n')
} else {
  suffix <- if (apoe4) '_withAPOE4' else '_noAPOE4'
}

in_path  <- paste0('dataset_all_percentile_BAI', suffix, '.csv')
out_csv  <- paste0('cutoff_results_lococv', suffix, '_', model_type, '.csv')
out_test_csv  <- paste0('cutoff_results_lococv_test', suffix, '_', model_type, '.csv')
out_folds_csv <- paste0('cutoff_results_lococv_folds', suffix, '_', model_type, '.csv')
out_boot_csv  <- paste0('cutoff_bootstrap_selection', suffix, '_', model_type, '.csv')
out_boot_hr_csv <- paste0('cutoff_bootstrap_HR', suffix, '_', model_type, '.csv')

# ---- Load data ----------------------------------------------------------
cat('Reading', in_path, '\n')
if (!file.exists(in_path))
  stop('Input not found: ', in_path, ' -- run step2_compute_percentile.py first')
df <- read.csv(in_path)
cat('  Total rows:', nrow(df), '\n')
df$event <- factor(df$event, levels = c('censor', 'dementia', 'death'))
df$id    <- seq_len(nrow(df))   # one row = one subject; used for the Fine-Gray cluster

cohorts <- sort(unique(df$dataset))
cat('  Cohorts    :', paste(cohorts, collapse = ', '), '\n')
cat('  Model type :', model_type, '\n')
cat('  Cutoffs    :', paste(cutoffs, collapse = ' '), '\n')
cat('  Horizons   :', paste(horizons, collapse = ' '), 'years\n\n')

# ---- Helper: fit chosen model, return HR/CI/pvalue ----------------------
# Fine-Gray notes:
#   strata(dataset) makes finegray estimate the censoring distribution (and hence
#   the IPCW weights) separately within each cohort, since follow-up and censoring
#   patterns differ substantially across ARIC/FHS/MESA/MrOS/SOF. It does not
#   stratify the baseline hazard of the Cox fit, so the estimand is unchanged.
#   cluster = id is required because finegray expands each subject into many rows;
#   without it the sandwich variance treats those rows as independent and the
#   standard error is too small (see ?finegray, citing Geskus 2011).
fit_model <- function(dat, mtype) {
  if (mtype == 'finegray') {
    pdata <- finegray(Surv(time2event, event) ~ BAIPercentileBinary + id + strata(dataset),
                      data = dat, etype = 'dementia')
    mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~ BAIPercentileBinary,
                 weight = fgwt, data = pdata, cluster = id)
  } else {
    dat$event_cs <- as.integer(dat$event == 'dementia')
    mod <- coxph(Surv(time2event, event_cs) ~ BAIPercentileBinary, data = dat)
  }
  s  <- summary(mod)$coefficients
  ci <- exp(confint(mod))
  list(HR     = exp(s[1, 'coef']),
       HR_lo  = ci[1, 1],
       HR_hi  = ci[1, 2],
       pvalue = s[1, 'Pr(>|z|)'])
}

# ---- Helper: fit the model at every candidate cutoff -----------------------
sweep_cutoffs <- function(dat, mtype, verbose = TRUE) {
  sweep <- data.frame(cutoff = cutoffs, HR = NA_real_, HR_lo = NA_real_,
                      HR_hi = NA_real_, pvalue = NA_real_)
  for (i in seq_along(cutoffs)) {
    cut <- cutoffs[i]
    dat$BAIPercentileBinary <- as.integer(dat$BAIPercentile >= cut)
    if (length(unique(dat$BAIPercentileBinary)) < 2) next
    tryCatch({
      res <- fit_model(dat, mtype)
      sweep[i, c('HR', 'HR_lo', 'HR_hi', 'pvalue')] <- c(res$HR, res$HR_lo, res$HR_hi, res$pvalue)
      if (verbose)
        cat(sprintf('  cutoff=%d  HR=%.3f [%.3f, %.3f]  p=%.4g\n',
                    cut, res$HR, res$HR_lo, res$HR_hi, res$pvalue))
    }, error = function(e) {
      if (verbose) cat(sprintf('  cutoff=%d: ERROR - %s\n', cut, conditionMessage(e)))
    })
  }
  sweep
}

# ---- Helper: time-dependent sensitivity / PPV / AUC of a binary cutoff ----
# Competing-risk definitions at horizon t (Blanche et al. 2013, as in timeROC):
#   cases     = dementia by t
#   non-cases = everyone else, i.e. still at risk at t or died (competing event) before t
#   Se  = P(BAIPercentile >= cut | case)
#   Sp  = P(BAIPercentile <  cut | non-case)
#   PPV = P(case | BAIPercentile >= cut)
#   AUC = (Se + Sp) / 2 for a binary marker
# Subjects censored before t are handled by inverse probability of censoring
# weighting (IPCW). The censoring distribution is a Kaplan-Meier estimate within
# each cohort, consistent with strata(dataset) in the Fine-Gray model.
# (timeROC::SeSpPPVNPV gives the same estimates with weighting = 'marginal', but
# cannot stratify the censoring model by cohort.)
ipcw_weights <- function(dat, t) {
  w <- numeric(nrow(dat))
  for (coh in unique(dat$dataset)) {
    idx  <- which(dat$dataset == coh)
    km   <- survfit(Surv(dat$time2event[idx], dat$event[idx] == 'censor') ~ 1)
    G    <- stepfun(km$time, c(1, km$surv))
    Gt   <- G(t)
    if (Gt <= 0)
      stop(sprintf('Horizon %g y exceeds follow-up in cohort %s', t, coh))
    Ti   <- dat$time2event[idx]
    ev   <- dat$event[idx]
    # G(T-) for events before t; G(t) for those still at risk at t; 0 if censored before t
    Gmin <- G(Ti - 1e-8)
    w[idx] <- ifelse(Ti > t, 1 / Gt, ifelse(ev != 'censor', 1 / Gmin, 0))
  }
  w
}

# Returns a matrix [cutoff x (Se, PPV, AUC)] for one horizon
clinical_metrics <- function(dat, t, cuts) {
  w    <- ipcw_weights(dat, t)
  case <- dat$time2event <= t & dat$event == 'dementia'
  t(sapply(cuts, function(cut) {
    pos <- dat$BAIPercentile >= cut
    se  <- sum(w[case & pos])  / sum(w[case])
    sp  <- sum(w[!case & !pos]) / sum(w[!case])
    ppv <- sum(w[case & pos])  / sum(w[pos])
    c(Se = se, PPV = ppv, AUC = (se + sp) / 2)
  }))
}

# Point estimates on the full data plus percentile bootstrap CIs (resampling
# subjects within cohort). Returns a data.frame with one row per cutoff x horizon.
clinical_metrics_boot <- function(dat, horizons, cuts, B) {
  est  <- lapply(horizons, function(t) clinical_metrics(dat, t, cuts))
  boot <- lapply(horizons, function(t) array(NA_real_, c(B, length(cuts), 3)))
  coh_idx <- split(seq_len(nrow(dat)), dat$dataset)
  for (b in seq_len(B)) {
    bi <- unlist(lapply(coh_idx, function(ix) ix[sample.int(length(ix), replace = TRUE)]))
    bd <- dat[bi, ]
    for (h in seq_along(horizons)) boot[[h]][b, , ] <- clinical_metrics(bd, horizons[h], cuts)
  }
  do.call(rbind, lapply(seq_along(horizons), function(h) {
    out <- data.frame(cutoff = cuts, horizon = horizons[h])
    for (m in seq_along(c('Se', 'PPV', 'AUC'))) {
      nm <- c('Se', 'PPV', 'AUC')[m]
      ci <- apply(boot[[h]][, , m, drop = FALSE], 2, quantile, c(0.025, 0.975), na.rm = TRUE)
      out[[nm]]            <- est[[h]][, m]
      out[[paste0(nm, '_lo')]] <- ci[1, ]
      out[[paste0(nm, '_hi')]] <- ci[2, ]
    }
    out
  }))
}

# ---- LOCO-CV ------------------------------------------------------------
fold_summary  <- list()   # per-fold best-cutoff info
held_out_list <- list()   # held-out rows annotated with BAIPercentileBinary
sweep_list    <- list()   # full training-set sweep per fold (all cutoffs)
test_sweep_list <- list() # held-out cohort sweep per fold (all cutoffs)

cat('=== LOCO-CV ===\n')
for (held_cohort in cohorts) {
  cat(sprintf('\n--- Fold: held-out = %s ---\n', held_cohort))

  train_df <- df[df$dataset != held_cohort, ]
  test_df  <- df[df$dataset == held_cohort, ]

  sweep <- sweep_cutoffs(train_df, model_type)

  valid <- sweep[!is.na(sweep$pvalue), ]
  if (nrow(valid) == 0) {
    cat('  No valid cutoffs found – skipping fold.\n')
    next
  }

  best <- valid[which.min(valid$pvalue), ]
  cat(sprintf('  --> Selected cutoff: %d  (train HR=%.3f [%.3f, %.3f], p=%.4g)\n',
              best$cutoff, best$HR, best$HR_lo, best$HR_hi, best$pvalue))

  # Store full sweep for this fold (all cutoffs, flagging the selected one)
  sweep$held_cohort <- held_cohort
  sweep$selected    <- sweep$cutoff == best$cutoff
  sweep_list[[held_cohort]] <- sweep[, c('held_cohort', 'cutoff', 'HR', 'HR_lo', 'HR_hi', 'pvalue', 'selected')]

  # Evaluate every cutoff in the held-out cohort (evaluation only; the cutoff is
  # still selected on the training cohorts above)
  cat('  Held-out cohort, all cutoffs:\n')
  test_sweep <- sweep_cutoffs(test_df, model_type)
  test_sweep$held_cohort <- held_cohort
  test_sweep$selected    <- test_sweep$cutoff == best$cutoff
  test_sweep_list[[held_cohort]] <- test_sweep[, c('held_cohort', 'cutoff', 'HR', 'HR_lo', 'HR_hi', 'pvalue', 'selected')]

  # Apply selected cutoff to held-out cohort
  test_df$BAIPercentileBinary <- as.integer(test_df$BAIPercentile >= best$cutoff)

  # Estimate HR in held-out cohort (informational; may be noisy in small cohorts)
  test_hr_info <- tryCatch({
    if (length(unique(test_df$BAIPercentileBinary)) < 2) stop('no variation')
    fit_model(test_df, model_type)
  }, error = function(e) {
    cat(sprintf('  (test-set HR unavailable: %s)\n', conditionMessage(e)))
    list(HR = NA, HR_lo = NA, HR_hi = NA, pvalue = NA)
  })

  fold_summary[[held_cohort]] <- data.frame(
    held_cohort  = held_cohort,
    best_cutoff  = best$cutoff,
    train_n      = nrow(train_df),
    train_HR     = best$HR,
    train_HR_lo  = best$HR_lo,
    train_HR_hi  = best$HR_hi,
    train_pvalue = best$pvalue,
    test_n       = nrow(test_df),
    test_HR      = test_hr_info$HR,
    test_HR_lo   = test_hr_info$HR_lo,
    test_HR_hi   = test_hr_info$HR_hi,
    test_pvalue  = test_hr_info$pvalue,
    test_n_dementia   = sum(test_df$event == 'dementia'),
    test_n_death      = sum(test_df$event == 'death'),
    test_person_years = sum(test_df$time2event),
    # median potential follow-up (reverse Kaplan-Meier)
    test_median_followup = unname(summary(survfit(Surv(time2event, event == 'censor') ~ 1,
                                                   data = test_df))$table['median']),
    stringsAsFactors = FALSE
  )

  held_out_list[[held_cohort]] <- test_df
}

# ---- Cross-validated HR (pooled held-out cohorts) -----------------------
cat('\n=== Cross-validated HR (pooled held-out cohorts) ===\n')
pooled <- do.call(rbind, held_out_list)
cat('  Pooled rows:', nrow(pooled),
    '  n_high:', sum(pooled$BAIPercentileBinary == 1),
    '  n_low:',  sum(pooled$BAIPercentileBinary == 0), '\n')

cv_res <- tryCatch(
  fit_model(pooled, model_type),
  error = function(e) {
    cat('ERROR fitting pooled model:', conditionMessage(e), '\n')
    NULL
  }
)

if (!is.null(cv_res)) {
  cat(sprintf('  CV HR = %.3f [%.3f, %.3f]  p = %.4g\n',
              cv_res$HR, cv_res$HR_lo, cv_res$HR_hi, cv_res$pvalue))
}

# ---- Save LOCO-CV results -----------------------------------------------
# Full training-set sweep across all folds (one row per fold × cutoff)
sweep_df <- do.call(rbind, sweep_list)

# Append CV-pooled summary row
cv_row <- NULL
if (!is.null(cv_res)) {
  cv_row <- data.frame(held_cohort = 'CV_pooled',
                       cutoff      = NA_integer_,
                       HR          = cv_res$HR,
                       HR_lo       = cv_res$HR_lo,
                       HR_hi       = cv_res$HR_hi,
                       pvalue      = cv_res$pvalue,
                       selected    = NA,
                       stringsAsFactors = FALSE)
}

write.csv(rbind(sweep_df, cv_row), out_csv, row.names = FALSE)
cat('LOCO-CV results saved to:', out_csv, '\n')

# Held-out cohorts of all folds pooled together, at every cutoff
cat('\n=== Pooled held-out cohorts, all cutoffs ===\n')
pooled_sweep <- sweep_cutoffs(pooled, model_type)
pooled_sweep$held_cohort <- 'Pooled'
pooled_sweep$selected    <- NA
test_sweep_list[['Pooled']] <- pooled_sweep[, c('held_cohort', 'cutoff', 'HR', 'HR_lo', 'HR_hi', 'pvalue', 'selected')]

write.csv(do.call(rbind, test_sweep_list), out_test_csv, row.names = FALSE)
cat('LOCO-CV held-out cohort sweep saved to:', out_test_csv, '\n')

fold_df <- do.call(rbind, fold_summary)
write.csv(fold_df, out_folds_csv, row.names = FALSE)
cat('LOCO-CV per-fold summary saved to:', out_folds_csv, '\n')
print(fold_df, row.names = FALSE)

# ---- Full-data sweep (all cohorts, for reference) -----------------------
cat('\n=== Full-data sweep (all cohorts, reference) ===\n')
results_full <- data.frame(cutoff = cutoffs,
                           HR     = NA_real_, HR_lo = NA_real_, HR_hi = NA_real_,
                           pvalue = NA_real_, n_high = NA_integer_, n_low = NA_integer_)

cat(sprintf('%-10s %-10s %-20s %-12s %-8s %-8s\n',
            'Cutoff', 'HR', '95% CI', 'p-value', 'n_high', 'n_low'))

for (i in seq_along(cutoffs)) {
  cut <- cutoffs[i]
  df$BAIPercentileBinary <- as.integer(df$BAIPercentile >= cut)
  tryCatch({
    res  <- fit_model(df, model_type)
    n_hi <- sum(df$BAIPercentileBinary == 1)
    n_lo <- sum(df$BAIPercentileBinary == 0)
    results_full[i, ] <- list(cut,
                              res$HR, res$HR_lo, res$HR_hi, res$pvalue,
                              n_hi, n_lo)
    cat(sprintf('%-10d %-10.3f [%-6.3f, %-6.3f]    %-12.4g %-8d %-8d\n',
                cut, res$HR, res$HR_lo, res$HR_hi, res$pvalue, n_hi, n_lo))
  }, error = function(e) {
    cat(sprintf('Cutoff %d: ERROR - %s\n', cut, conditionMessage(e)))
  })
}

best_full <- results_full[which.min(results_full$pvalue), ]
cat(sprintf('\nBest cutoff (full data): %d  (HR=%.3f, p=%.4g)\n',
            best_full$cutoff, best_full$HR, best_full$pvalue))

# ---- Bootstrap of the cutoff selection (full data) ------------------------
# Resample subjects within cohort, repeat the sweep, and record which cutoff has
# the smallest p-value. Shows how stable the selected cutoff is to sampling.
# The HRs from the same replicates give a bootstrap 95% CI of the HR at each cutoff.
if (n_boot_select > 0) {
  cat(sprintf('\n=== Bootstrap of cutoff selection (%d replicates, %d cores) ===\n',
              n_boot_select, n_cores))
  coh_idx <- split(seq_len(nrow(df)), df$dataset)
  boot_select_one <- function(b) {
    bi <- unlist(lapply(coh_idx, function(ix) ix[sample.int(length(ix), replace = TRUE)]))
    bd <- df[bi, ]
    bd$id <- seq_len(nrow(bd))   # each resampled row is its own subject
    sweep_cutoffs(bd, model_type, verbose = FALSE)
  }
  RNGkind("L'Ecuyer-CMRG")
  set.seed(2024)
  boot_sweeps <- parallel::mclapply(seq_len(n_boot_select), boot_select_one,
                                    mc.cores = n_cores, mc.set.seed = TRUE)
  boot_sel <- sapply(boot_sweeps, function(sw)
    if (all(is.na(sw$pvalue))) NA_integer_ else sw$cutoff[which.min(sw$pvalue)])
  boot_hr  <- sapply(boot_sweeps, function(sw) sw$HR)   # cutoff x replicate
  boot_tab <- data.frame(cutoff     = cutoffs,
                         n_selected = sapply(cutoffs, function(cc) sum(boot_sel == cc, na.rm = TRUE)))
  boot_tab$pct_selected <- round(100 * boot_tab$n_selected / sum(!is.na(boot_sel)), 1)
  write.csv(boot_tab, out_boot_csv, row.names = FALSE)
  cat(sprintf('  failed replicates: %d\n', sum(is.na(boot_sel))))
  print(boot_tab, row.names = FALSE)
  cat('Bootstrap selection frequencies saved to:', out_boot_csv, '\n')

  boot_hr_ci <- t(apply(boot_hr, 1, quantile, c(0.025, 0.975), na.rm = TRUE))
  boot_hr_df <- data.frame(cutoff = cutoffs, HR = results_full$HR,
                           HR_boot_lo = boot_hr_ci[, 1], HR_boot_hi = boot_hr_ci[, 2])
  write.csv(boot_hr_df, out_boot_hr_csv, row.names = FALSE)
  print(boot_hr_df, row.names = FALSE)
  cat('Bootstrap HR 95% CIs saved to:', out_boot_hr_csv, '\n')
}

# ---- Sensitivity / PPV / AUC for each cutoff (full data) ----------------
cat(sprintf('\n=== Sensitivity / PPV / AUC (full data, %d stratified bootstrap replicates) ===\n', n_boot))
set.seed(2024)
clin <- clinical_metrics_boot(df, horizons, cutoffs, n_boot)
fmt_pct_ci <- function(x, lo, hi) sprintf('%.1f%% (%.1f–%.1f)', 100 * x, 100 * lo, 100 * hi)
fmt_auc_ci <- function(x, lo, hi) sprintf('%.3f (%.3f–%.3f)', x, lo, hi)

out_clin <- paste0('cutoff_clinical_metrics_all', suffix, '.csv')
write.csv(clin, out_clin, row.names = FALSE)
cat('Raw Se/PPV/AUC saved to:', out_clin, '\n')

# Formatted Table 1 (full-data sweep)
fmt_pval <- function(p) if (is.na(p)) 'NA' else if (p < 0.001) sprintf('%.2e', p) else sprintf('%.3f', p)
model_label <- if (model_type == 'finegray') 'SHR' else 'HR'
table1 <- data.frame(
  Cutoff  = paste0(results_full$cutoff, 'th'),
  HR_col  = sprintf('%.2f', results_full$HR),
  CI_95   = sprintf('%.2f–%.2f', results_full$HR_lo, results_full$HR_hi),
  p_value = sapply(results_full$pvalue, fmt_pval),
  optimal = ifelse(!is.na(results_full$cutoff) & results_full$cutoff == best_full$cutoff, '*', ''),
  stringsAsFactors = FALSE
)
colnames(table1)[2] <- model_label
for (t in horizons) {
  ct <- clin[clin$horizon == t, ]
  ct <- ct[match(results_full$cutoff, ct$cutoff), ]
  table1[[sprintf('Sensitivity_%gy', t)]] <- fmt_pct_ci(ct$Se,  ct$Se_lo,  ct$Se_hi)
  table1[[sprintf('PPV_%gy', t)]]         <- fmt_pct_ci(ct$PPV, ct$PPV_lo, ct$PPV_hi)
  table1[[sprintf('AUC_%gy', t)]]         <- fmt_auc_ci(ct$AUC, ct$AUC_lo, ct$AUC_hi)
}

out_table1 <- paste0('cutoff_results_all', suffix, '_', model_type, '.csv')
write.csv(table1, out_table1, row.names = FALSE)
cat('Table 1 (full-data) saved to:', out_table1, '\n')

# Console print
cat(sprintf('\nTable 1. %s model results (full data)\n',
            if (model_type == 'finegray') 'Fine-Gray' else 'Cause-specific Cox PH'))
print(table1, row.names = FALSE)
