library(survival)
library(tdROC)

# Help message
show_help <- function() {
  cat("Usage: Rscript step3_analysis-survival.R <exposure> [data.type] [stratification]\n\n")
  cat("Arguments:\n")
  cat("  exposure       Exposure variable name (required)\n")
  cat("                 Examples: BAI, alpha_bandpower_mean_C_N1, delta_theta_mean_C_N3\n\n")
  cat("  data.type      Data filtering type (optional, default: withoutAPOE)\n")
  cat("                 Options: withAPOE, withoutAPOE, withoutAPOE-matched\n\n")
  cat("  stratification Subgroup analysis (optional, default: none)\n")
  cat("                 Options: young, old, male, female\n\n")
  cat("Examples:\n")
  cat("  Rscript step3_analysis-survival.R BAI\n")
  cat("  Rscript step3_analysis-survival.R BAI withAPOE\n")
  cat("  Rscript step3_analysis-survival.R BAI withoutAPOE male\n")
}

dataset <- 'ARIC'
args <- commandArgs(trailingOnly = TRUE)

# Check for help or missing arguments
if (length(args) == 0 || args[1] %in% c('-h', '--help', 'help')) {
  show_help()
  quit(status = 0)
}

exposure <- args[[1]]
data.type <- if (length(args) >= 2) args[[2]] else 'withoutAPOE'
stratification <- if (length(args) >= 3) args[[3]] else ''

# Validate arguments
valid_data_types <- c('withAPOE', 'withoutAPOE', 'withoutAPOE-matched')
valid_stratifications <- c('', 'young', 'old', 'male', 'female')

if (!data.type %in% valid_data_types) {
  cat("Error: Invalid data.type. Must be one of:", paste(valid_data_types, collapse = ', '), "\n")
  quit(status = 1)
}

if (!stratification %in% valid_stratifications) {
  cat("Error: Invalid stratification. Must be one of:", paste(valid_stratifications, collapse = ', '), "\n")
  quit(status = 1)
}

cat(sprintf("Running analysis with exposure=%s, data.type=%s, stratification=%s\n", 
            exposure, data.type, stratification))

# Get current working directory and construct paths dynamically
current_dir <- getwd()
df <- read.csv(file.path(current_dir, 'dataset_ARIC.csv'))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))

df <- df[df$prevalent_dementia==0,]

if (exposure=='BAI') {
   df[,exposure] <- df[,exposure]/10
} else {
    df[,exposure] <- df[,exposure]/sd(df[,exposure], na.rm=T)
}

if (stratification=='young') {
  df <- df[df$age<70,]
} else if (stratification=='old') {
  df <- df[df$age>=70,]
} else if (stratification=='male') {
  df <- df[df$sexM==1,]
} else if (stratification=='female') {
  df <- df[df$sexM==0,]
}
print(dim(df))

formula.str <- sprintf('Surv(time2event, event) ~ %s + age + sexM', exposure)
pdata <- finegray(as.formula(formula.str), data=df, etype='dementia')
formula.str <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + sexM', exposure)
model1 <- coxph(as.formula(formula.str), weight=fgwt, data=pdata)
sm <- summary(model1)
print(sm)
print(AIC(model1))

res1 <- c(Model='AgeOnly',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)



if (data.type=='withAPOE') {
  cols <- c(exposure, 'age', 'sexM', 'educcollege', 'BMI', 'sleepmed', 'exercise', 'APOE4Count')
} else if (data.type=='withoutAPOE-matched') {
  cols <- c(exposure, 'age', 'sexM', 'educcollege', 'BMI', 'sleepmed', 'exercise', 'APOE4Count')
  df <- df[complete.cases(df[,c('time2event','event',cols)]),]
  cols <- c(exposure, 'age', 'sexM', 'educcollege', 'BMI', 'sleepmed', 'exercise')
} else {
  cols <- c(exposure, 'age', 'sexM', 'educcollege', 'BMI', 'sleepmed', 'exercise')
}
print(colSums(is.na(df[,c('time2event','event',cols)])))
df <- df[complete.cases(df[,c('time2event','event',cols)]),]
formula.str1 <- sprintf('Surv(time2event, event) ~ %s', paste0(cols, collapse='+'))
formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s', paste0(cols, collapse='+'))
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model2 <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model2)
print(sm)
print(AIC(model2))

res2 <- c(Model='Intermediate',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)#TODO N=
print(sprintf('%s: Intermediate model N = %d', dataset, nrow(df)))
#write.csv(df, file.path(current_dir, sprintf('dataset_%s_table1.csv', dataset)))

res <- survfit(Surv(time2event, event)~1, data=df)
png(file.path(current_dir, sprintf('AJ_%s.png', dataset)), width = 3*300, height = 1.5*300)
plot(res, xlab='Years since sleep study', ylab='Cumulative incidence', conf.int=F, cumhaz=1)
dev.off()

if (data.type=='withAPOE') {
  formula.str1 <- sprintf('Surv(time2event, event) ~ %s + age + sexM + educcollege + BMI + race_Black + race_Other + sleepmed + exercise + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI + APOE4Count', exposure)
  formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + sexM + educcollege + BMI + race_Black + race_Other + sleepmed + exercise + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI + APOE4Count', exposure)
} else {
  formula.str1 <- sprintf('Surv(time2event, event) ~ %s + age + sexM + educcollege + BMI + race_Black + race_Other + sleepmed + exercise + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI', exposure)
  formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + sexM + educcollege + BMI + race_Black + race_Other + sleepmed + exercise + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI', exposure)
}
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model3 <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model3)
print(sm)
print(AIC(model3))

res3 <- c(Model='Full',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)

res <- rbind(res1, res2, res3)
print(res)
write.csv(res, file.path(current_dir, sprintf('BAI_results_ARIC-survival-%s-%s%s.csv', exposure, data.type, stratification)), row.names=F)

