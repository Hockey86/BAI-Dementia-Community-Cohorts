library(meta)
library(metafor)
library(metap)

args <- commandArgs(trailingOnly = TRUE)
col.interaction <- ifelse(length(args) >= 1, args[1], 'APOE4')

# Prepare the data
data.dir <- "/data/haoqisun/BAI_dementia_community"
res.dir <- file.path(data.dir, 'meta-analysis', 'interaction')

mesa_file <- file.path(data.dir, sprintf('MESA/interaction_results_MESA_%s-survival.csv', col.interaction))
df.mesa <- if (file.exists(mesa_file)) read.csv(mesa_file) else NULL

aric_file <- file.path(data.dir, sprintf('SHHS/ARIC/interaction_results_ARIC_%s-survival.csv', col.interaction))
df.aric <- if (file.exists(aric_file)) read.csv(aric_file) else NULL

fhs_file <- file.path(data.dir, sprintf('SHHS/FHS/interaction_results_FHS_%s-survival.csv', col.interaction))
df.fhs <- if (file.exists(fhs_file)) read.csv(fhs_file) else NULL

mros_file <- file.path(data.dir, sprintf('MrOS/interaction_results_MrOS_%s-survival.csv', col.interaction))
df.mros <- if (file.exists(mros_file)) read.csv(mros_file) else NULL

sof_file <- file.path(data.dir, sprintf('SOF/interaction_results_SOF_%s-survival.csv', col.interaction))
df.sof <- if (file.exists(sof_file)) read.csv(sof_file) else NULL

scale <- 1

all_studies <- c("MESA", "ARIC", "FHS", "MrOS", "SOF")

hrs <- c()
ci_lowers <- c()
ci_uppers <- c()
p_values <- c()
available_studies <- c()

for (study in all_studies) {
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
  
  if (!is.null(df_)) {
    hrs <- c(hrs, df_$exp.coef.)
    ci_lowers <- c(ci_lowers, df_$lower..95)
    ci_uppers <- c(ci_uppers, df_$upper..95)
    p_values <- c(p_values, df_$Pr...z..)
    available_studies <- c(available_studies, study)
  } else {
    cat(sprintf("Warning: Interaction file for study %s with %s does not exist, skipping...\n", study, col.interaction))
  }
}
data <- data.frame(study=available_studies, hr=hrs, ci_lower=ci_lowers, ci_upper=ci_uppers, p_value=p_values)

# Convert Odds Ratios to Log Odds Ratios
data$log_hr <- log(data$hr)*scale
data$log_hr_se <- (log(data$ci_upper) - log(data$ci_lower)) / (2 * 1.96)*scale

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
writeLines(summary_text, file.path(res.dir, sprintf("meta_analysis_%s_interaction_result-survival.txt", col.interaction)))
print(summary_text)

#png(sprintf("forest_plot_%s.png", model.type), width = 700*2, height = 300*2)#, dpi=300)
#forest(meta_result, plotwidth="24cm")
#dev.off()
forest(meta_result, common=F,
       file=file.path(res.dir, sprintf("forest_plot_%s_interaction-survival.png", col.interaction)),
       width=1000, fontsize=16, digits=3)#plotwidth="24cm")
