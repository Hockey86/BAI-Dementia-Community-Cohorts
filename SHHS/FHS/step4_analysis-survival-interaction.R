library(survival)


dataset <- 'FHS'
current_dir <- getwd()
df <- read.csv(file.path(current_dir, sprintf('dataset_%s.csv', dataset)))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))
df <- df[df$prevalent_dementia==0,]

data.type <- 'withoutAPOE'
cols <- c('sexM', 'educcollege', 'BMI', 'race_Other', 'sleepmed', 'pascore')
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
write.csv(t(res), file.path(current_dir, sprintf('age_interaction_results_%s-survival-%s.csv', dataset, data.type)), row.names=F)


cols <- c('age', 'educcollege', 'BMI', 'race_Other', 'sleepmed', 'pascore')

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
write.csv(t(res), file.path(current_dir, sprintf('sex_interaction_results_%s-survival-%s.csv', dataset, data.type)), row.names=F)
