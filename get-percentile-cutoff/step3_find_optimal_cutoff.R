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
cat('  Cutoffs    :', paste(cutoffs, collapse = ' '), '\n\n')

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

# ---- LOCO-CV ------------------------------------------------------------
fold_summary  <- list()   # per-fold best-cutoff info
held_out_list <- list()   # held-out rows annotated with BAIPercentileBinary
sweep_list    <- list()   # full training-set sweep per fold (all cutoffs)

cat('=== LOCO-CV ===\n')
for (held_cohort in cohorts) {
  cat(sprintf('\n--- Fold: held-out = %s ---\n', held_cohort))

  train_df <- df[df$dataset != held_cohort, ]
  test_df  <- df[df$dataset == held_cohort, ]

  sweep <- data.frame(cutoff  = cutoffs,
                      HR      = NA_real_,
                      HR_lo   = NA_real_,
                      HR_hi   = NA_real_,
                      pvalue  = NA_real_)

  for (i in seq_along(cutoffs)) {
    cut <- cutoffs[i]
    train_df$BAIPercentileBinary <- as.integer(train_df$BAIPercentile >= cut)
    if (length(unique(train_df$BAIPercentileBinary)) < 2) next
    tryCatch({
      res              <- fit_model(train_df, model_type)
      sweep$HR[i]     <- res$HR
      sweep$HR_lo[i]  <- res$HR_lo
      sweep$HR_hi[i]  <- res$HR_hi
      sweep$pvalue[i] <- res$pvalue
      cat(sprintf('  cutoff=%d  HR=%.3f [%.3f, %.3f]  p=%.4g\n',
                  cut, res$HR, res$HR_lo, res$HR_hi, res$pvalue))
    }, error = function(e) {
      cat(sprintf('  cutoff=%d: ERROR - %s\n', cut, conditionMessage(e)))
    })
  }

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

out_table1 <- paste0('cutoff_results_all', suffix, '_', model_type, '.csv')
write.csv(table1, out_table1, row.names = FALSE)
cat('Table 1 (full-data) saved to:', out_table1, '\n')

# Console print
cat(sprintf('\nTable 1. %s model results (full data)\n',
            if (model_type == 'finegray') 'Fine-Gray' else 'Cause-specific Cox PH'))
cat(sprintf('%-10s %-8s %-18s %-10s\n', 'Cutoff', model_label, '95% CI', 'p-value'))
cat(strrep('-', 50), '\n')
for (i in seq_len(nrow(table1))) {
  marker <- if (table1$optimal[i] == '*') ' <-- optimal' else ''
  cat(sprintf('%-10s %-8s %-18s %-10s%s\n',
              table1$Cutoff[i], table1[i, 2], table1$CI_95[i],
              table1$p_value[i], marker))
}
