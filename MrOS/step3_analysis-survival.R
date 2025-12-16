library(survival)
library(cmprsk)
library(survminer)
library(tdROC)
library(ggplot2)
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 || any(c('-h', '--help') %in% args)) {
  cat("Usage: Rscript step3_analysis-survival.R <exposure> [data.type] [stratification]\n")
  cat("\nArguments:\n")
  cat("  exposure       Exposure variable name (required)\n")
  cat("  data.type      Data type: 'withAPOE', 'withoutAPOE', 'withoutAPOE-matched' (default: 'withoutAPOE')\n")
  cat("  stratification Stratification: '', 'young', 'old' (default: '') - Note: MrOS is male-only\n")
  cat("\nExample:\n")
  cat("  Rscript step3_analysis-survival.R BAI\n")
  cat("  Rscript step3_analysis-survival.R BAI withAPOE young\n")
  quit(save = "no", status = 0)
}

if (length(args) < 1) {
  stop("Error: exposure argument is required. Use -h for help.")
}

exposure <- args[1]
data.type <- ifelse(length(args) >= 2, args[2], 'withoutAPOE')
stratification <- ifelse(length(args) >= 3, args[3], '')

valid_data_types <- c('withAPOE', 'withoutAPOE', 'withoutAPOE-matched')
valid_stratifications <- c('', 'young', 'old')

if (!data.type %in% valid_data_types) {
  stop(sprintf("Error: Invalid data.type '%s'. Valid options: %s", data.type, paste(valid_data_types, collapse=", ")))
}

if (!stratification %in% valid_stratifications) {
  stop(sprintf("Error: Invalid stratification '%s'. Valid options: %s (Note: MrOS is male-only)", stratification, paste(valid_stratifications, collapse=", ")))
}

dataset <- 'MrOS'
folder <- file.path(getwd())
df <- read.csv(file.path(folder, sprintf('dataset_%s.csv', dataset)))
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
}
print(dim(df))

#df$event2 <- 0
#df$event2[df$event=='censor'] <- 0
#df$event2[df$event=='dementia'] <- 1
#df$event2[df$event=='death'] <- 2
#df$group <- ''
#df$group[df$BAI>quantile(df$BAI,2/3)] <- 'Third quartile in BAI'
#df$group[df$BAI<quantile(df$BAI,1/3)] <- 'First quartile in BAI'

#df2 <- df[df$group!='',]

#cif <- cuminc(ftime = df2$time2event, fstatus = df2$event2, group = df2$group)
#cif.filter <- cif[grepl("1", names(cif))]
#names(cif.filter) <- gsub("1", "Cognitive Impairment", names(cif.filter))
#ggcompetingrisks(cif.filter, conf.int = TRUE,
#                 xlab = "Year since sleep test", ylab = "Cumulative Incidence (%)",
#                 ggtheme = theme_minimal())


#model1 <- crr(ftime = df$time2event, fstatus = df$event2, cov1 = df[,c('BAI', 'age')], failcode=1, cencode=0)

# New data frame for predictions
#new_data <- df[1:10,c('BAI','age')]
#new_data <- new_data %>% 
#  mutate(pred = predict(model1, cov1 = new_data))

# Convert to long format for plotting
#long_data <- new_data %>%
#  pivot_longer(cols = starts_with("pred"), names_to = "covariate", values_to = "cif")

# Plot CIF using ggplot2
#ggplot(long_data, aes(x = time, y = cif, color = covariate)) +
#  geom_line() +
#  labs(title = "Cumulative Incidence Function",
#       x = "Time",
#       y = "Cumulative Incidence",
#       color = "Covariate") +
#  theme_minimal()


formula.str <- sprintf('Surv(time2event, event) ~ %s', exposure)
pdata <- finegray(as.formula(formula.str), data=df, etype='dementia')
formula.str <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s', exposure)
model0 <- coxph(as.formula(formula.str), weight=fgwt, data=pdata)
sm <- summary(model0)
print(sm)
print(c(AIC(model0), BIC(model0)))

res0 <- c(Model='Unadjusted',
          sm$coefficients[exposure, c('exp(coef)', 'Pr(>|z|)')],
          sm$conf.int[exposure, c('lower .95', 'upper .95')],
          Np=sm$nevent)

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
          Np=sm$nevent)#TODO N=


if (data.type=='withAPOE') {
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'race_NonWhite', 'sleepmed', 'pascore', 'smoke_current', 'APOE4Count')
} else if (data.type=='withoutAPOE-matched') {
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'race_NonWhite', 'sleepmed', 'pascore', 'smoke_current', 'APOE4Count')
    df <- df[complete.cases(df[,c('time2event','event',cols)]),]
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'race_NonWhite', 'sleepmed', 'pascore', 'smoke_current')
    
} else {
    cols <- c(exposure, 'age', 'educcollege', 'BMI', 'race_NonWhite', 'sleepmed', 'pascore', 'smoke_current')
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

#res <- survfit(Surv(time2event, event)~1, data=df)
#png(file.path(folder, sprintf("AJ_%s.png", dataset)), width = 3*300, height = 1.5*300)
#plot(res, xlab='Years since sleep study', ylab='Cumulative incidence', conf.int=F, cumhaz=1)
#dev.off()

# Extract time points and cumulative events
#cum_events <- data.frame(
#  time = res$time,
#  n.event = cumsum(res$n.event)  # cumulative events
#)

# Make cumulative hazard plot
#p <- ggsurvplot(
#  res,
#  fun = "cumhaz",
#  conf.int = FALSE,
#  xlab = "Years since sleep study",
#  ylab = "Cumulative incidence",
#  risk.table = FALSE # we'll add our own "events table"
#)
#print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")

# Create a table of cumulative events
#tbl <- ggplot(cum_events, aes(x = time, y = 0, label = n.event)) +
#  geom_text(vjust = 1.5) +
#  theme_void() +
#  scale_x_continuous(limits = range(p$plot$data$time))
#print("bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb")
## Arrange plot and table
#final_plot <- p$plot / tbl + plot_layout(heights = c(3, 1))

#png(file.path(folder, sprintf("AJ_%s.png", dataset)), width = 3*300, height = 1.5*300)
#print(final_plot)
#dev.off()



if (data.type=='withAPOE') {
    formula.str1 <- sprintf('Surv(time2event, event) ~ %s + age + educcollege + BMI + race_NonWhite + sleepmed + pascore + smoke_current + diabetes + hypertension + heartattack + stroke + depression + tms + AHI + APOE4Count', exposure)
    formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + educcollege + BMI + race_NonWhite + sleepmed + pascore + smoke_current + diabetes + hypertension + heartattack + stroke + depression + tms + AHI + APOE4Count', exposure)
} else {
    formula.str1 <- sprintf('Surv(time2event, event) ~ %s + age + educcollege + BMI + race_NonWhite + sleepmed + pascore + smoke_current + diabetes + hypertension + heartattack + stroke + depression + tms + AHI', exposure)
    formula.str2 <- sprintf('Surv(fgstart, fgstop, fgstatus) ~ %s + age + educcollege + BMI + race_NonWhite + sleepmed + pascore + smoke_current + diabetes + hypertension + heartattack + stroke + depression + tms + AHI', exposure)
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

res <- rbind(res0, res1, res2, res3)
print(res)
write.csv(res, file.path(folder, sprintf('BAI_results_%s-survival-%s-%s%s.csv', dataset, exposure, data.type, stratification)), row.names=F)

