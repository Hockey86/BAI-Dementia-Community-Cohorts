import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from step1_generate_healthy_ref_dataset import *
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


# Fine-Gray model (dementia as event, death as competing risk), same as in step3_find_optimal_cutoff.R:
# strata(dataset) estimates the censoring weights within each cohort; cluster = id gives robust SE.
ro.r('''
suppressMessages(library(survival))
fit_finegray <- function(dat) {
  dat$event <- factor(dat$event, levels = c('censor', 'dementia', 'death'))
  dat$id    <- seq_len(nrow(dat))
  pdata <- finegray(Surv(time2event, event) ~ exposure + id + strata(dataset),
                    data = dat, etype = 'dementia')
  mod <- coxph(Surv(fgstart, fgstop, fgstatus) ~ exposure,
               weight = fgwt, data = pdata, cluster = id)
  s  <- summary(mod)$coefficients
  ci <- exp(confint(mod))
  c(exp(s[1, 'coef']), ci[1, 1], ci[1, 2], s[1, 'Pr(>|z|)'])
}
''')

def fit_finegray(dat):
    with localconverter(ro.default_converter + pandas2ri.converter):
        rdat = ro.conversion.py2rpy(dat[['dataset','event','time2event','exposure']])
    return list(ro.globalenv['fit_finegray'](rdat))

def format_hr(hr, lo, hi, p):
    p_str = 'p<0.001' if p<0.001 else f'p={p:.3f}'
    return f'{hr:.2f} [{lo:.2f}, {hi:.2f}], {p_str}'


df = pd.read_csv('dataset_all_percentile_BAI_noAPOE4.csv')
df['BAI'] = df.BAI*10 # the BAI was divided by 10 in the original paper, recover it
BAI_median = 0.
cutoff = 90
BAI_perc = np.percentile(df.BAI, cutoff)
#BAI_median = df.BAI.median()
#print(BAI_median)

dfs = [
    load_mesa(healthy=False),
    load_mros(healthy=False),
    load_sof(healthy=False),
    load_fhs(healthy=False),
    load_aric(healthy=False),
]
df2 = pd.concat(dfs, ignore_index=True, axis=0)
df2['BAI'] = df2.BAI*10 # the BAI was divided by 10 in the original paper, recover it
df = df.merge(df2[['dataset', 'BAI', 'age', 'sexM', 'BMI','AHI', 'healthy_indicator']], on=['dataset','BAI'],how='inner',validate='1:1')

df['logAHI'] = np.log1p(df.AHI.values)
X = df.loc[df.healthy_indicator==1, ['age','sexM','BMI','logAHI']].values
y = df.BAI[df.healthy_indicator==1].values
Xmean=X.mean(axis=0)
Xstd=X.std(axis=0)
assert Xstd.min()>0
X2=(X-Xmean)/Xstd
model=LinearRegression().fit(X2,y)
df['BAI_res'] = df.BAI.values-model.predict((df[['age','sexM','BMI','logAHI']].values-Xmean)/Xstd)
BAI_res_perc = np.percentile(df.BAI_res, cutoff)

rows = []

groups = [
    (f'BAI>=0 & Percentile>={cutoff}', (df.BAI>=BAI_median) & (df.BAIPercentile>=cutoff)),
    (f'BAI>=0 & Percentile<{cutoff}',  (df.BAI>=BAI_median) & (df.BAIPercentile<cutoff)),
    (f'BAI<0 & Percentile>={cutoff}',  (df.BAI<BAI_median)  & (df.BAIPercentile>=cutoff)),
    (f'BAI<0 & Percentile<{cutoff}',   (df.BAI<BAI_median)  & (df.BAIPercentile<cutoff)),
    (f'Percentile>={cutoff}', (df.BAIPercentile>=cutoff)),
    (f'Percentile<{cutoff}', (df.BAIPercentile<cutoff)),
    (f'BAI>=0',   (df.BAI>=BAI_median) ),
    (f'BAI<0',   (df.BAI<BAI_median) ),
    (f'BAI>=BAI {cutoff} Percentile ({BAI_perc:.1f}y)',   (df.BAI>=BAI_perc) ),
    (f'BAI<BAI {cutoff} Percentile ({BAI_perc:.1f}y)',   (df.BAI<BAI_perc) ),
    (f'BAI_res>=BAI_res {cutoff} Percentile ({BAI_res_perc:.1f}y)',   (df.BAI_res>=BAI_res_perc) ),
    (f'BAI_res<BAI_res {cutoff} Percentile ({BAI_res_perc:.1f}y)',   (df.BAI_res<BAI_res_perc) ),
]

# compare each pair of rows (1st vs 2nd, 3rd vs 4th, ...): exposure = 1st row, reference = 2nd row
assert len(groups)%2==0
hr_strs = []
for (label1, mask1), (label2, mask2) in zip(groups[0::2], groups[1::2]):
    assert not (mask1 & mask2).any()
    dat = df[mask1 | mask2].copy()
    dat['exposure'] = mask1[mask1 | mask2].astype(int).values
    hr_strs.extend([format_hr(*fit_finegray(dat)), 'ref'])

def n_pct(m, denom):
    return f'{m.sum()} ({m.sum()/denom*100:.1f}%)'

for (label, mask), hr_str in zip(groups, hr_strs):
    n = mask.sum()
    rows.append({
        'group':            label,
        'N (%)':            n_pct(mask, len(df)),
        'Dementia, n (%)':  n_pct(df.event[mask] == "dementia", n),
        'Death, n (%)':     n_pct(df.event[mask] == "death", n),
        'HR [95% CI], p':   hr_str,
    })

output_file = f'reclassification_results_noAPOE4_cutoff{cutoff}.csv'
results = pd.DataFrame(rows)
print(results)
results.to_csv(output_file, index=False)
print(f'Saved to {output_file}')