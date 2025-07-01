import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})
import seaborn as sns
sns.set_style('ticks')


feat_names = [
'BAI',
'COUPL_OVERLAP_C', 'DENS_C', 'sigma_bandpower_kurtosis_C_N2',

'delta_bandpower_kurtosis_C_N2',
'theta_bandpower_kurtosis_C_N2',
'alpha_bandpower_kurtosis_C_N2',
'kurtosis_N2_C',

'delta_bandpower_mean_C_N3',
'delta_alpha_mean_C_N3',
'delta_theta_mean_C_N3',
'theta_bandpower_kurtosis_C_N3',
'kurtosis_N3_C',

#'alpha_bandpower_mean_C_N1',
]
feat_names2 = [
'BAI',
'spindle-SO overlap', 'spindle density', 'sigma power kurtosis during N2',

'delta power kurtosis during N2',
'theta power kurtosis during N2',
'alpha power kurtosis during N2',
'signal kurtosis during N2',

'delta power during N3',
'delta/alpha ratio during N3',
'delta/theta ratio during N3',
'theta power kurtosis during N3',
'signal kurtosis during N3'

#'alpha_bandpower_mean_C_N1',
]

df = {'Name':[], 'HR':[], 'LB':[], 'UB':[], 'P':[]}
for fn in feat_names:
    with open(os.path.join('withoutAPOE', f'meta_analysis_result_{fn}-Intermediate-survival.txt'), 'r') as ff:
        for l in ff:
            if l.startswith('Random effects model'):
                #try:
                xx = l[len('Random effects model '):].split()
                df['Name'].append(fn)
                df['HR'].append( float(xx[0]) )
                df['LB'].append( float(xx[1][1:-1]) )
                df['UB'].append( float(xx[2][:-1]) )
                df['P'].append(  ''.join(xx[4:]) )
                if df['P'][-1]=='<0.0001':
                    df['P'][-1] = 0.
                else:
                    df['P'][-1] = float(df['P'][-1])
                #except Exception as ee:
                #    print(fn, l, str(ee))
                break

df = pd.DataFrame(data=df)
print(df)
breakpoint()
df.to_csv('BAI_feature_coefs.csv', index=False)

