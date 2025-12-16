library(survival)


dataset <- 'SOF'
folder <- file.path(getwd())
df <- read.csv(file.path(folder, sprintf('dataset_%s.csv', dataset)))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))
df <- df[df$prevalent_dementia==0,]

cols <- c('age', 'educcollege', 'BMI', 'sleepmed', 'walking', 'smoke_current')

formula.str1 <- sprintf('Surv(time2event, event) ~ BAI*APOE4 + %s', paste0(cols, collapse='+'))
formula.str2 <- gsub('Surv(time2event, event)', 'Surv(fgstart, fgstop, fgstatus)', formula.str1, fixed=T)
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model)
print(sm)
print(AIC(model))

res <- c(sm$coefficients['BAI:APOE4', c('exp(coef)', 'Pr(>|z|)')],
         sm$conf.int['BAI:APOE4', c('lower .95', 'upper .95')],
         Np=sm$nevent, N=nrow(df))
print(res)
write.csv(t(res), file.path(folder, sprintf('interaction_results_%s_APOE4-survival.csv', dataset)), row.names=F)

