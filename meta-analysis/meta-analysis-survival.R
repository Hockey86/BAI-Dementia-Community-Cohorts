library(meta)
library(metap)

args <- commandArgs(trailingOnly = TRUE)
exposure <- args[1]
data.type <- ifelse(length(args) >= 2, args[2], 'withoutAPOE')
stratification <- ifelse(length(args) >= 3, args[3], '')
#save_fig_format <- ifelse(length(args) >= 4, args[4], 'png')  # 'png' or 'pdf'
save_fig_format <- 'png'


# Prepare the data
data.dir <- "/data/haoqisun/BAI_dementia_community"
res.dir <- file.path(data.dir, 'meta-analysis', data.type)

mesa_file <- file.path(data.dir, sprintf('MESA/BAI_results_MESA-survival-%s-%s%s.csv', exposure, data.type, stratification))
df.mesa <- if (file.exists(mesa_file)) read.csv(mesa_file) else NULL

aric_file <- file.path(data.dir, sprintf('SHHS/ARIC/BAI_results_ARIC-survival-%s-%s%s.csv', exposure, data.type, stratification))
df.aric <- if (file.exists(aric_file)) read.csv(aric_file) else NULL

fhs_file <- file.path(data.dir, sprintf('SHHS/FHS/BAI_results_FHS-survival-%s-%s%s.csv', exposure, data.type, stratification))
df.fhs <- if (file.exists(fhs_file)) read.csv(fhs_file) else NULL

if (stratification=='male') {
    mros_file <- file.path(data.dir, sprintf('MrOS/BAI_results_MrOS-survival-%s-%s%s.csv', exposure, data.type, ''))
} else {
    mros_file <- file.path(data.dir, sprintf('MrOS/BAI_results_MrOS-survival-%s-%s%s.csv', exposure, data.type, stratification))
}
df.mros <- if (file.exists(mros_file)) read.csv(mros_file) else NULL

if (stratification=='female' || stratification=='old') {
    sof_file <- file.path(data.dir, sprintf('SOF/BAI_results_SOF-survival-%s-%s%s.csv', exposure, data.type, ''))
} else {
    sof_file <- file.path(data.dir, sprintf('SOF/BAI_results_SOF-survival-%s-%s%s.csv', exposure, data.type, stratification))
}
df.sof <- if (file.exists(sof_file)) read.csv(sof_file) else NULL

if (exposure=='age') {
    model.types <- c('age')
} else {
    model.types <- c('AgeOnly', 'Intermediate', 'Full')
}
scale <- 1

all_studies <- c("MESA", "ARIC", "FHS-OS", "MrOS", "SOF")
for (model.type in model.types) {
  hrs <- c()
  ci_lowers <- c()
  ci_uppers <- c()
  p_values <- c()
  available.studies <- c()
  available.cols.square <- c()
  available.cols.study <- c()
  
  for (study in all_studies) {
      if (study=='MrOS') {
          df_ <- df.mros
          col.square <- "orchid"
          col.study <- "purple"
      } else if (study=='SOF') {
          df_ <- df.sof
          col.square <- "gold"
          col.study <- "orange"
      } else if (study=='MESA') {
          df_ <- df.mesa
          col.square <- "pink"
          col.study <- "red"
      } else if (study=='FHS-OS') {
          df_ <- df.fhs
          col.square <- "lightgreen"
          col.study <- "green"
      } else if (study=='ARIC') {
          df_ <- df.aric
          col.square <- "lightblue"
          col.study <- "blue"
      }
      
      if (!is.null(df_)) {
          hrs <- c(hrs, df_$exp.coef.[df_$Model==model.type])
          ci_lowers <- c(ci_lowers, df_$lower..95[df_$Model==model.type])
          ci_uppers <- c(ci_uppers, df_$upper..95[df_$Model==model.type])
          p_values <- c(p_values, df_$Pr...z..[df_$Model==model.type])
          available.studies <- c(available.studies, study)
          available.cols.square <- c(available.cols.square, col.square)
          available.cols.study <- c(available.cols.study, col.study)
      } else {
          cat(sprintf("Warning: File for study %s does not exist, skipping...\n", study))
      }
  }
  data <- data.frame(study=available.studies, hr=hrs, ci_lower=ci_lowers, ci_upper=ci_uppers, p_value=p_values)
  
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
  writeLines(summary_text, file.path(res.dir, sprintf("meta_analysis_result-survival-%s-%s%s.txt", exposure, model.type, stratification)))
  print(summary_text)
  
  save.path <- file.path(res.dir, sprintf("forest_plot_survival-%s-%s%s.%s", exposure, model.type, stratification, save_fig_format))

  if (save_fig_format == 'pdf') {
    pdf(save.path, width=8, height=3)
  } else {
    png(save.path, width=2400, height=900, res=300)
  }

  forest(meta_result, common=F,
         col.square = available.cols.square,
         col.study = available.cols.study,
         print.pval = TRUE,      # Print p-values
         leftcols = c("studlab", "effect", "ci", "pval"),
         leftlabs = c("Study", "HR", "95% CI", "P"),
         rightcols= c("w.random"),
         rightlabs= c("Weight"),
         digits.pval = 3,
         spacing=1.3,
         fontsize=14, digits=3)#plotwidth="24cm")
  dev.off()
  ##file=save.path,
  
}


