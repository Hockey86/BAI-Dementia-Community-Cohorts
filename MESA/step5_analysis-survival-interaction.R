library(survival)


dataset <- 'MESA'
folder <- file.path(getwd())
df <- read.csv(file.path(folder, sprintf('dataset_%s.csv', dataset)))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))
df <- df[df$prevalent_dementia==0,]

#data.type <- 'withoutAPOE'
cols <- c('sexM', 'race_Chinese', 'race_Black', 'race_Hispanic', 'educcollege', 'BMI', 'sleepmed', 'walk_min_per_wk', 'smoke_current')
df <- df[complete.cases(df[,c('time2event','event','BAI','age',cols)]),]
df[,'age_group'] <- as.integer(df$age>=70)
print(dim(df))

formula.str1 <- sprintf('Surv(time2event, event) ~ BAI*age_group + %s', paste0(cols, collapse='+'))
formula.str2 <- gsub('Surv(time2event, event)', 'Surv(fgstart, fgstop, fgstatus)', formula.str1, fixed=T)
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model)
print(sm)
print(AIC(model))

res <- c(sm$coefficients['BAI:age_group', c('exp(coef)', 'Pr(>|z|)')],
         sm$conf.int['BAI:age_group', c('lower .95', 'upper .95')],
         Np=sm$nevent, N=nrow(df))
print(res)
write.csv(t(res), file.path(folder, sprintf('interaction_results_%s_age-survival.csv', dataset)), row.names=F)


cols <- c('age', 'race_Chinese', 'race_Black', 'race_Hispanic', 'educcollege', 'BMI', 'sleepmed', 'walk_min_per_wk', 'smoke_current')

formula.str1 <- sprintf('Surv(time2event, event) ~ BAI*sexM + %s', paste0(cols, collapse='+'))
formula.str2 <- gsub('Surv(time2event, event)', 'Surv(fgstart, fgstop, fgstatus)', formula.str1, fixed=T)
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model)
print(sm)
print(AIC(model))

res <- c(sm$coefficients['BAI:sexM', c('exp(coef)', 'Pr(>|z|)')],
         sm$conf.int['BAI:sexM', c('lower .95', 'upper .95')],
         Np=sm$nevent, N=nrow(df))
print(res)
write.csv(t(res), file.path(folder, sprintf('interaction_results_%s_sex-survival.csv', dataset)), row.names=F)


cols <- c('age', 'sexM', 'race_Chinese', 'race_Black', 'race_Hispanic', 'educcollege', 'BMI', 'sleepmed', 'walk_min_per_wk', 'smoke_current')

formula.str1 <- sprintf('Surv(time2event, event) ~ BAI*APOE4Count + %s', paste0(cols, collapse='+'))
formula.str2 <- gsub('Surv(time2event, event)', 'Surv(fgstart, fgstop, fgstatus)', formula.str1, fixed=T)
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model)
print(sm)
print(AIC(model))

res <- c(sm$coefficients['BAI:APOE4Count', c('exp(coef)', 'Pr(>|z|)')],
         sm$conf.int['BAI:APOE4Count', c('lower .95', 'upper .95')],
         Np=sm$nevent, N=nrow(df))
print(res)
write.csv(t(res), file.path(folder, sprintf('interaction_results_%s_APOE4-survival.csv', dataset)), row.names=F)

