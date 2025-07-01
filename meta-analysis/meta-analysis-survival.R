library(meta)
library(metap)

# Prepare the data
data.dir <- "/data/haoqisun/BAI_dementia_community"
data.type <- "withoutAPOE"
stratification <- ''
res.dir <- file.path(data.dir, 'meta-analysis', data.type)

df.mesa <- read.csv(file.path(data.dir, sprintf('MESA/BAI_results_MESA-survival-%s%s.csv', data.type, stratification)))
df.aric <- read.csv(file.path(data.dir, sprintf('SHHS/ARIC/BAI_results_ARIC-survival-%s%s.csv', data.type, stratification)))
df.fhs <- read.csv(file.path(data.dir, sprintf('SHHS/FHS/BAI_results_FHS-survival-%s%s.csv', data.type, stratification)))
df.mros <- read.csv(file.path(data.dir, sprintf('MrOS/BAI_results_MrOS-survival-%s%s.csv', data.type, stratification)))
df.sof <- read.csv(file.path(data.dir, sprintf('SOF/BAI_results_SOF-survival-%s%s.csv', data.type, stratification)))

model.types <- c('AgeOnly', 'Intermediate')#, 'Full')
scale <- 1

studies <- c("MESA", "ARIC", "FHS", "MrOS", "SOF")
for (model.type in model.types) {
  hrs <- c()
  ci_lowers <- c()
  ci_uppers <- c()
  p_values <- c()
  for (study in studies) {
      if (study=='MrOS') {
          df_ <- df.mros
      } else if (study=='SOF') {
          df_ <- df.sof
      } else if (study=='MESA') {
          df_ <- df.mesa
      } else if (study=='FHS') {
          df_ <- df.fhs
      } else if (study=='ARIC') {
          df_ <- df.aric
      }
      hrs <- c(hrs, df_$exp.coef.[df_$Model==model.type])
      ci_lowers <- c(ci_lowers, df_$lower..95[df_$Model==model.type])
      ci_uppers <- c(ci_uppers, df_$upper..95[df_$Model==model.type])
      p_values <- c(p_values, df_$Pr...z..[df_$Model==model.type])
  }
  data <- data.frame(study=studies, hr=hrs, ci_lower=ci_lowers, ci_upper=ci_uppers, p_value=p_values)
  
  # Convert Odds Ratios to Log Odds Ratios
  data$log_hr <- log(data$hr)*scale
  data$log_hr_se <- (log(data$ci_upper) - log(data$ci_lower)) / (2 * 1.96)*scale
  
  print('#######')
  print(model.type)
  print(data)
  
  # Perform the meta-analysis
  meta_result <- metagen(TE = data$log_hr, 
                         seTE = data$log_hr_se, 
                         studlab = data$study,
                         sm = "HR")
  
  summary_meta <- summary(meta_result)
  summary_text <- capture.output(summary_meta)
  
  # Combine p-values using Fisher's method
  p_combined <- metap::sumlog(data$p_value)
  summary_text2 <- capture.output(p_combined)
  
  summary_text <- c(summary_text, summary_text2)
  writeLines(summary_text, file.path(res.dir, sprintf("meta_analysis_result_%s%s-survival.txt", model.type, stratification)))
  #print(summary_text)
  
  save.path <- file.path(res.dir, sprintf("forest_plot_%s%s-survival.png", model.type, stratification))
  png(save.path, width=2400, height=900, res=300)
  forest(meta_result, common=F,
         col.square = c("pink", "lightblue", "lightgreen", "orchid", "gold"),
         col.study = c("red", "blue", "green", "purple", "orange"),
         print.pval = TRUE,      # Print p-values
         leftcols = c("studlab", "effect", "ci", "pval"), 
         leftlabs = c("Study", "HR", "95% CI", "P"),
         rightcols= c("w.random"),
         rightlabs= c("Weight"),
         digits.pval = 3,
         spacing=1.3,
         fontsize=14, digits=3)#plotwidth="24cm")
  dev.off()
  #file=save.path,
  
}


