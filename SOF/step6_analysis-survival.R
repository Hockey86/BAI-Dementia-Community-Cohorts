library(survival)
library(tdROC)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || any(c('-h', '--help') %in% args)) {
  cat("Usage: Rscript step6_analysis-survival.R <exposure> [data.type]\n")
  cat("\nArguments:\n")
  cat("  exposure       Exposure variable name (required)\n")
  cat("  data.type      Data type: 'withAPOE', 'withoutAPOE', 'withoutAPOE-matched' (default: 'withoutAPOE')\n")
  cat("\nNote: SOF is female-only and elderly population, so no stratification options\n")
  cat("\nExample:\n")
  cat("  Rscript step6_analysis-survival.R BAI\n")
  cat("  Rscript step6_analysis-survival.R BAI withAPOE\n")
  quit(save = "no", status = 0)
}

if (length(args) < 1) {
  stop("Error: exposure argument is required. Use -h for help.")
}

exposure <- args[1]
data.type <- ifelse(length(args) >= 2, args[2], 'withoutAPOE')

valid_data_types <- c('withAPOE', 'withoutAPOE', 'withoutAPOE-matched')

if (!data.type %in% valid_data_types) {
  stop(sprintf("Error: Invalid data.type '%s'. Valid options: %s", data.type, paste(valid_data_types, collapse=", ")))
}

dataset <- 'SOF'
folder <- file.path(getwd())
df <- read.csv(file.path(folder, sprintf('dataset_%s.csv', dataset)))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))
df <- df[df$prevalent_dementia==0,]

if (exposure=='BAI') {
   df[,exposure] <- df[,exposure]/10
} else {
    df[,exposure] <- df[,exposure]/sd(df[,exposure], na.rm=T)
}

print(dim(df))

formula.str <- sprintf('Surv(time2event, event) ~ %s + age', exposure)
pdata <- finegray(as.formula(formula.str), data=df, etype='dementia')
formula.str <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age', exposure)
model1 <- coxph(as.formula(formula.str), weight=fgwt, data=pdata)
sm <- summary(model1)
print(sm)
print(c(AIC(model1), BIC(model1)))

res1 <- c(Model='AgeOnly',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)


if (data.type=='withAPOE') {
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'sleepmed', 'walking', 'smoke_current', 'APOE4')
} else if (data.type=='withoutAPOE-matched') {
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'sleepmed', 'walking', 'smoke_current', 'APOE4')
    df <- df[complete.cases(df[,c('time2event','event',cols)]),]
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'sleepmed', 'walking', 'smoke_current')
    
} else {
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'sleepmed', 'walking', 'smoke_current')
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
print(c(AIC(model2), BIC(model2)))

#formula.str1 <- sprintf('Surv(time2event, event) ~ %s+I(BAI^2)+I(BAI^3)', paste0(cols, collapse='+'))
#formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s+I(BAI^2)+I(BAI^3)', paste0(cols, collapse='+'))
#formula1 <- as.formula(formula.str1)
#formula2 <- as.formula(formula.str2)
#pdata <- finegray(formula1, data=df, etype='dementia')
#model2p <- coxph(formula2, weight=fgwt, data=pdata)
#sm <- summary(model2p)
#print(sm)
#print(c(AIC(model2p), BIC(model2p)))

res2 <- c(Model='Intermediate',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)#TODO N=
print(sprintf('%s: Intermediate model N = %d', dataset, nrow(df)))
write.csv(df, file.path(folder, sprintf('dataset_%s_table1.csv', dataset)))

res <- survfit(Surv(time2event, event)~1, data=df)
png(file.path(folder, sprintf("AJ_%s.png", dataset)), width = 3*300, height = 1.5*300)
plot(res, xlab='Years since sleep study', ylab='Cumulative incidence', conf.int=F, cumhaz=1)
dev.off()


if (data.type=='withAPOE') {
    formula.str1 <- sprintf('Surv(time2event, event) ~ %s + age + educcollege + BMI + sleepmed + walking + smoke_current +diabetes + hypertension + heartattack + stroke + depression + mmse + AHI + APOE4', exposure)
    formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + educcollege + BMI + sleepmed + walking + smoke_current + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI + APOE4', exposure)
} else {
    formula.str1 <- sprintf('Surv(time2event, event) ~ %s + age + educcollege + BMI + sleepmed + walking + smoke_current + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI', exposure)
    formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + educcollege + BMI + sleepmed + walking + smoke_current + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI', exposure)
}
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model3 <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model3)
print(sm)
print(c(AIC(model3), BIC(model3)))

res3 <- c(Model='Full',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)

res <- rbind(res1, res2, res3)
print(res)
write.csv(res, file.path(folder, sprintf('BAI_results_%s-survival-%s-%s.csv', dataset, exposure, data.type)), row.names=F)

