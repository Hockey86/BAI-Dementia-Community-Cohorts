Rscript meta-analysis-survival.R BAI
Rscript meta-analysis-survival.R BAI withAPOE
Rscript meta-analysis-survival.R BAI withoutAPOE-matched
Rscript meta-analysis-survival.R BAI withoutAPOE young
Rscript meta-analysis-survival.R BAI withoutAPOE old
Rscript meta-analysis-survival.R BAI withoutAPOE male
Rscript meta-analysis-survival.R BAI withoutAPOE female

for x in COUPL_OVERLAP_C DENS_C alpha_bandpower_kurtosis_C_N2 alpha_bandpower_mean_C_N1 delta_alpha_mean_C_N3 delta_bandpower_kurtosis_C_N2 delta_bandpower_mean_C_N3 delta_theta_mean_C_N3 kurtosis_N2_C kurtosis_N3_C sigma_bandpower_kurtosis_C_N2 theta_bandpower_kurtosis_C_N2 theta_bandpower_kurtosis_C_N3
do
    Rscript meta-analysis-survival.R $x
done

Rscript meta-analysis-survival-interaction.R age
Rscript meta-analysis-survival-interaction.R sex
Rscript meta-analysis-survival-interaction.R APOE4
