library(meta)
library(metafor)
library(metap)

# Prepare the data
data.dir <- "/data/haoqisun/BAI_dementia_community"

df.mesa <- read.csv(file.path(data.dir, 'MESA/BAI_results_MESA-binary.csv'))
df.aric <- read.csv(file.path(data.dir, 'SHHS/ARIC/BAI_results_ARIC-binary.csv'))
df.fhs <- read.csv(file.path(data.dir, 'SHHS/FHS/BAI_results_FHS-binary.csv'))
df.mros <- read.csv(file.path(data.dir, 'MrOS/BAI_results_MrOS-binary.csv'))
df.sof <- read.csv(file.path(data.dir, 'SOF/BAI_results_SOF-binary.csv'))

model.types <- c('AgeOnly', 'Mediate', 'Full')

for (model.type in model.types) {
  data <- data.frame(
    study = c("MESA", "ARIC", "FHS", "MrOS", "SOF"),
    or = c(df.mesa$OddsRatio[df.mesa$Model==model.type],
           df.aric$OddsRatio[df.aric$Model==model.type],
           df.fhs$OddsRatio[df.fhs$Model==model.type],
           df.mros$OddsRatio[df.mros$Model==model.type],
           df.sof$OddsRatio[df.sof$Model==model.type]),
    ci_lower = c(df.mesa$LB[df.mesa$Model==model.type],
           df.aric$LB[df.aric$Model==model.type],
           df.fhs$LB[df.fhs$Model==model.type],
           df.mros$LB[df.mros$Model==model.type],
           df.sof$LB[df.sof$Model==model.type]),
    ci_upper = c(df.mesa$UB[df.mesa$Model==model.type],
                 df.aric$UB[df.aric$Model==model.type],
                 df.fhs$UB[df.fhs$Model==model.type],
                 df.mros$UB[df.mros$Model==model.type],
                 df.sof$UB[df.sof$Model==model.type]),
    p_value = c(df.mesa$P[df.mesa$Model==model.type],
                 df.aric$P[df.aric$Model==model.type],
                 df.fhs$P[df.fhs$Model==model.type],
                 df.mros$P[df.mros$Model==model.type],
                 df.sof$P[df.sof$Model==model.type])
  )
  
  # Convert Odds Ratios to Log Odds Ratios
  data$log_or <- log(data$or)
  data$log_or_se <- (log(data$ci_upper) - log(data$ci_lower)) / (2 * 1.96)
  
  print('#######')
  print(model.type)
  print(data)
  
  # Perform the meta-analysis
  meta_result <- metagen(TE = data$log_or, 
                         seTE = data$log_or_se, 
                         studlab = data$study,
                         sm = "OR")
  
  summary_meta <- summary(meta_result)
  summary_text <- capture.output(summary_meta)
  
  # Combine p-values using Fisher's method
  p_combined <- metap::sumlog(data$p_value)
  summary_text2 <- capture.output(p_combined)
  
  summary_text <- c(summary_text, summary_text2)
  writeLines(summary_text, sprintf("meta_analysis_result_%s-binary.txt", model.type))
  #print(summary_text)
  
  #png(sprintf("forest_plot_%s.png", model.type), width = 700*2, height = 300*2)#, dpi=300)
  #forest(meta_result, plotwidth="24cm")
  #dev.off()
  forest(meta_result, file=sprintf("forest_plot_%s-binary.png", model.type), width=1000, fontsize=16, digits=3)#plotwidth="24cm")
  
}




