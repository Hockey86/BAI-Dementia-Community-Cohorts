#!/usr/bin/env Rscript
# Step 3: Find optimal BAI percentile cutoff that maximizes association
# with incident dementia (minimizes Fine-Gray competing-risk HR p-value).
#
# Usage:
#   Rscript step3_find_optimal_cutoff.R                  # noAPOE4 version
#   Rscript step3_find_optimal_cutoff.R --apoe4          # withAPOE4 version
#   Rscript step3_find_optimal_cutoff.R --cutoffs 50:95  # custom range (default 50:95)

library(survival)
library(cmprsk)

# ---- Parse CLI arguments -----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

apoe4    <- '--apoe4' %in% args
cutoff_arg <- grep('^--cutoffs=', args, value = TRUE)
if (length(cutoff_arg) > 0) {
  range_str <- sub('--cutoffs=', '', cutoff_arg)
  parts <- as.integer(strsplit(range_str, ':')[[1]])
  cutoffs <- seq(parts[1], parts[2])
} else {
  cutoffs <- seq(50,90,5)
}

suffix   <- if (apoe4) '_withAPOE4' else '_noAPOE4'
in_path  <- paste0('dataset_all_percentile_BAI', suffix, '.csv')
out_csv  <- paste0('optimal_cutoff_results', suffix, '.csv')
out_plot <- paste0('optimal_cutoff_pvalue', suffix, '.pdf')

# ---- Load data ---------------------------------------------------------
cat('Reading', in_path, '\n')
df <- read.csv(in_path)
cat('  Total rows:', nrow(df), '\n')

df <- df[df$n_similar >= 50, ]
cat('  Rows with n_similar >= 50:', nrow(df), '\n')

df$event <- factor(df$event, levels = c('censor', 'dementia', 'death'))

# ---- Sweep cutoffs -----------------------------------------------------
results <- data.frame(
  cutoff  = integer(0),
  HR      = numeric(0),
  HR_lo   = numeric(0),
  HR_hi   = numeric(0),
  pvalue  = numeric(0),
  n_high  = integer(0),
  n_low   = integer(0)
)

cat(sprintf('%-10s %-10s %-20s %-12s %-8s %-8s\n',
            'Cutoff', 'HR', '95% CI', 'p-value', 'n_high', 'n_low'))

for (cutoff in cutoffs) {
  df$BAIPercentileBinary <- as.integer(df$BAIPercentile >= cutoff)

  tryCatch({
    formula.str <- 'Surv(time2event, event) ~ BAIPercentileBinary'
    pdata <- finegray(as.formula(formula.str), data = df, etype = 'dementia')

    formula.str2 <- 'Surv(fgstart, fgstop, fgstatus) ~ BAIPercentileBinary'
    model <- coxph(as.formula(formula.str2), weight = fgwt, data = pdata)

    s    <- summary(model)$coefficients
    ci   <- exp(confint(model))
    hr   <- exp(s[1, 'coef'])
    pval <- s[1, 'Pr(>|z|)']
    n_hi <- sum(df$BAIPercentileBinary == 1)
    n_lo <- sum(df$BAIPercentileBinary == 0)

    results <- rbind(results, data.frame(
      cutoff  = cutoff,
      HR      = hr,
      HR_lo   = ci[1, 1],
      HR_hi   = ci[1, 2],
      pvalue  = pval,
      n_high  = n_hi,
      n_low   = n_lo
    ))

    cat(sprintf('%-10d %-10.3f [%-6.3f, %-6.3f]    %-12.4g %-8d %-8d\n',
                cutoff, hr, ci[1,1], ci[1,2], pval, n_hi, n_lo))
  }, error = function(e) {
    cat(sprintf('Cutoff %d: ERROR - %s\n', cutoff, conditionMessage(e)))
  })
}

# ---- Best cutoff -------------------------------------------------------
best_idx    <- which.min(results$pvalue)
best_cutoff <- results$cutoff[best_idx]
best_hr     <- results$HR[best_idx]
best_pval   <- results$pvalue[best_idx]

cat(sprintf('\n=== Best cutoff: %d  (HR = %.3f, p = %.4g) ===\n',
            best_cutoff, best_hr, best_pval))

# ---- Save raw results --------------------------------------------------
write.csv(results, out_csv, row.names = FALSE)
cat('Results saved to:', out_csv, '\n')

# ---- Save formatted Table 1 --------------------------------------------
fmt_pval <- function(p) {
  if (p < 0.001) sprintf('%.2e', p)
  else           sprintf('%.3f', p)
}

table1 <- data.frame(
  Cutoff   = paste0(results$cutoff, 'th'),
  SHR      = sprintf('%.2f', results$HR),
  CI_95    = sprintf('%.2f\u2013%.2f', results$HR_lo, results$HR_hi),
  p_value  = sapply(results$pvalue, fmt_pval),
  stringsAsFactors = FALSE
)
# Mark the optimal row
table1$optimal <- ifelse(results$cutoff == best_cutoff, '*', '')

out_table1 <- paste0('cutoff_results', suffix, '.csv')
write.csv(table1, out_table1, row.names = FALSE)
cat('Table 1 saved to:', out_table1, '\n')

# Print Table 1 to console
cat('\nTable 1. Fine-Gray model results across candidate cutoffs\n')
cat(sprintf('%-10s %-8s %-18s %-10s\n', 'Cutoff', 'SHR', '95% CI', 'p-value'))
cat(strrep('-', 50), '\n')
for (i in seq_len(nrow(table1))) {
  marker <- if (table1$optimal[i] == '*') ' <-- optimal' else ''
  cat(sprintf('%-10s %-8s %-18s %-10s%s\n',
              table1$Cutoff[i], table1$SHR[i], table1$CI_95[i],
              table1$p_value[i], marker))
}

# ---- Plot p-value vs cutoff --------------------------------------------
pdf(out_plot, width = 6, height = 4)
plot(results$cutoff, -log10(results$pvalue),
     type  = 'b', pch = 16, col = 'steelblue',
     xlab  = 'BAI Percentile Cutoff',
     ylab  = expression(-log[10](p-value)),
     main  = paste0('Association with incident dementia\n(Fine-Gray, ', suffix, ')'))
abline(v   = best_cutoff, col = 'red', lty = 2)
abline(h   = -log10(0.05), col = 'gray', lty = 3)
legend('topright',
       legend = c(paste0('Best cutoff = ', best_cutoff),
                  'p = 0.05'),
       col    = c('red', 'gray'),
       lty    = c(2, 3), bty = 'n')
dev.off()
cat('Plot saved to:', out_plot, '\n')
