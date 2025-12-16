library(survival)


exposure <- 'age'
dataset <- 'ARIC'
data.type <- 'withoutAPOE'
stratification <- ''
folder <- file.path(getwd())
df <- read.csv(file.path(folder, sprintf('dataset_%s.csv', dataset)))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))
df <- df[df$prevalent_dementia==0,]

df[,exposure] <- df[,exposure]/10

print(dim(df))

formula.str <- sprintf('Surv(time2event, event) ~ %s', exposure)
pdata <- finegray(as.formula(formula.str), data=df, etype='dementia')
formula.str <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s', exposure)
model0 <- coxph(as.formula(formula.str), weight=fgwt, data=pdata)
sm <- summary(model0)
print(sm)
print(c(AIC(model0), BIC(model0)))

res <- c(Model='age',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)
res <- rbind(res)

print(res)
write.csv(res, file.path(folder, sprintf('BAI_results_%s-survival-%s-%s%s.csv', dataset, exposure, data.type, stratification)), row.names=F)

