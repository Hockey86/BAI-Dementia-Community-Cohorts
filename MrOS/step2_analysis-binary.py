import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import f_oneway
import statsmodels.api as sm
import statsmodels.formula.api as smf


def main():
    df = pd.read_csv('dataset_MrOS.csv')
    df = df[df.prevalent_dementia==0].reset_index(drop=True)
    #import matplotlib.pyplot as plt
    #plt.scatter(df.age, df.BA,s=4,c='k')
    #plt.plot([30,120],[30,120],c='r')
    #plt.show()
    
    df_res = defaultdict(list)
    """
    q1, q2 = np.nanpercentile(df.BAI, (100/3, 200/3))
    
    ids = {}
    ids[0] = df.BAI<q1
    ids[1] = (df.BAI<q2)&(df.BAI>=q1)
    ids[2] = df.BAI>=q2
    
    df_res['Name'].append('N')
    for i in range(3):
        df_res[f'Q{i+1}'].append(ids[i].sum())
    df_res['P'].append(np.nan)
        
    df_res['Name'].append('Age')
    col = 'vsage1'
    for i in range(3):
        df_ = df[ids[i]].reset_index(drop=True)
        df_res[f'Q{i+1}'].append(f'{df_.vsage1.mean():.1f} ({df_.vsage1.std():.1f})')
    df_res['P'].append(f_oneway(df.loc[ids[0],col].dropna(), df.loc[ids[1],col].dropna(), df.loc[ids[2],col].dropna()).pvalue)
    """
    
    model2formula = {
    'AgeOnly':'inc_dementia ~ BAI + age',
    
    'Mediate':'inc_dementia ~ BAI + age + educcollege + BMI + race_NonWhite + sleepmed + pascore + APOE4Count',
    
    'Full':   'inc_dementia ~ BAI + age + educcollege + BMI + race_NonWhite + sleepmed + pascore + diabetes + hypertension + heartattack + stroke + depression + tms + AHI + APOE4Count',
    }
    for mn, formula in model2formula.items():
        print(mn)
        model = smf.logit(formula=formula, data=df).fit()
        print(model.summary())
        df_res['Model'].append(mn)
        df_res['OddsRatio'].append(np.exp(model.params.BAI))
        ci = model.conf_int()
        df_res['LB'].append(np.exp(ci.loc['BAI',0]))
        df_res['UB'].append(np.exp(ci.loc['BAI',1]))
        df_res['P'].append(model.pvalues.BAI)
        df_res['N'].append(model.nobs)
        cols = [formula.split('~')[0]]+formula.split('~')[-1].split('+')
        cols = [x.strip() for x in cols]
        df2 = df[cols].dropna()
        assert len(df2)==model.nobs
        df_res['N+'].append((df2.inc_dementia==1).sum())
    df_res = pd.DataFrame(data=df_res)
    
    print(df_res)
    df_res.to_csv('BAI_results_MrOS-binary.csv', index=False)
    
        

if __name__=='__main__':
    main()

