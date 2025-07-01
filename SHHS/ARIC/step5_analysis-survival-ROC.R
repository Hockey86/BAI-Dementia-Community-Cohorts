library(survival)
library(cmprsk)
library(survminer)
library(tdROC)
library(xlsx)


nbt <- 1000
dataset <- 'ARIC'
current_dir <- getwd()
df <- read.csv(file.path(current_dir, sprintf('dataset_%s.csv', dataset)))
df$event <- factor(df$event, levels=c('censor', 'dementia', 'death'))
df <- df[df$prevalent_dementia==0,]
data.type <- 'withoutAPOE'

cols.basic <- c('age', 'sexM', 'educcollege', 'BMI', 'sleepmed', 'exercise')
df <- df[complete.cases(df[,c('time2event','event','BAI', cols.basic)]),]
print(dim(df))


print('===================')
print('BAI')
z <- df$BAI
roc.res0 <- tdROC.cr(as.numeric(z), df$time2event, as.integer(df$event)-1, nbt, method='integral', output='AUC', nboot=nbt)
youden.idx <- which.max(roc.res0$main_res$ROC$sens-(1-roc.res0$main_res$ROC$specB))
perf0 <- c(AUC=roc.res0$main_res$AUC.B.integral,
           AUC.lb=roc.res0$boot_res$bAUC.B.integral$CIlow,
           AUC.ub=roc.res0$boot_res$bAUC.B.integral$CIhigh,
           sens.youden=roc.res0$main_res$ROC$sens[youden.idx],
           spec.youden=roc.res0$main_res$ROC$specB[youden.idx],
           cutoff.youden=roc.res0$main_res$ROC$cut.off[youden.idx])
perf0 <- t(as.data.frame(perf0))
print(perf0)



print('===================')
print('cov')
cols <- cols.basic
formula.str1 <- sprintf('Surv(time2event, event) ~ %s', paste0(cols, collapse='+'))
formula.str2 <- gsub('Surv(time2event, event)', 'Surv(fgstart, fgstop, fgstatus)', formula.str1, fixed=T)
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model1 <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model1)

coef <- t(matrix(coefficients(model1)))
coef <- matrix(rep(t(coef),nrow(df)), ncol=ncol(coef), byrow=T)
z <- rowSums(coef*df[,cols])
roc.res1 <- tdROC.cr(as.numeric(z), df$time2event, as.integer(df$event)-1, nbt, method='integral', output='AUC', nboot=nbt)
youden.idx <- which.max(roc.res1$main_res$ROC$sens-(1-roc.res1$main_res$ROC$specB))
perf1 <- c(AUC=roc.res1$main_res$AUC.B.integral,
           AUC.lb=roc.res1$boot_res$bAUC.B.integral$CIlow,
           AUC.ub=roc.res1$boot_res$bAUC.B.integral$CIhigh,
           sens.youden=roc.res1$main_res$ROC$sens[youden.idx],
           spec.youden=roc.res1$main_res$ROC$specB[youden.idx],
           cutoff.youden=roc.res1$main_res$ROC$cut.off[youden.idx])
perf1 <- t(as.data.frame(perf1))
print(perf1)



print('===================')
print('BAI+cov')
cols <- c('BAI', cols.basic)
formula.str1 <- sprintf('Surv(time2event, event) ~ %s', paste0(cols, collapse='+'))
formula.str2 <- gsub('Surv(time2event, event)', 'Surv(fgstart, fgstop, fgstatus)', formula.str1, fixed=T)
formula1 <- as.formula(formula.str1)
formula2 <- as.formula(formula.str2)
pdata <- finegray(formula1, data=df, etype='dementia')
model2 <- coxph(formula2, weight=fgwt, data=pdata)
sm <- summary(model2)

coef <- t(matrix(coefficients(model2)))
coef <- matrix(rep(t(coef),nrow(df)), ncol=ncol(coef), byrow=T)
z <- rowSums(coef*df[,cols])
roc.res2 <- tdROC.cr(as.numeric(z), df$time2event, as.integer(df$event)-1, nbt, method='integral', output='AUC', nboot=nbt)
youden.idx <- which.max(roc.res2$main_res$ROC$sens-(1-roc.res2$main_res$ROC$specB))
perf2 <- c(AUC=roc.res2$main_res$AUC.B.integral,
           AUC.lb=roc.res2$boot_res$bAUC.B.integral$CIlow,
           AUC.ub=roc.res2$boot_res$bAUC.B.integral$CIhigh,
           sens.youden=roc.res2$main_res$ROC$sens[youden.idx],
           spec.youden=roc.res2$main_res$ROC$specB[youden.idx],
           cutoff.youden=roc.res2$main_res$ROC$cut.off[youden.idx])
perf2 <- t(as.data.frame(perf2))
print(perf2)

perf <- rbind(BAI_only=perf0, cov_only=perf1, BAI_cov=perf2)

save.path <- file.path(current_dir, sprintf("ROC_results_%s_nbt%d.xlsx", dataset, nbt))
write.xlsx(perf, file=save.path, sheetName="AUC", row.names=T)
write.xlsx(roc.res0$main_res$ROC, file=save.path, append=TRUE, sheetName="BAI-only-ROC", row.names=FALSE)
write.xlsx(roc.res0$boot_res$bROC, file=save.path, append=TRUE, sheetName=sprintf("BAI-only-ROC-bt%d",nbt), row.names=FALSE)
write.xlsx(roc.res1$main_res$ROC, file=save.path, append=TRUE, sheetName="cov-only-ROC", row.names=FALSE)
write.xlsx(roc.res1$boot_res$bROC, file=save.path, append=TRUE, sheetName=sprintf("cov-only-ROC-bt%d",nbt), row.names=FALSE)
write.xlsx(roc.res2$main_res$ROC, file=save.path, append=TRUE, sheetName="BAI_cov-ROC", row.names=FALSE)
write.xlsx(roc.res2$boot_res$bROC, file=save.path, append=TRUE, sheetName=sprintf("BAI_cov-ROC-bt%d",nbt), row.names=FALSE)
