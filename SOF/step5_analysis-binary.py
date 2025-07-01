import numpy as np
import pandas as pd
from collections import defaultdict
from scipy.stats import f_oneway
import statsmodels.api as sm
import statsmodels.formula.api as smf


def main():
    df = pd.read_csv('dataset_SOF.csv')
    print(df.shape)
    df = df[df.prevalent_dementia==0].reset_index(drop=True)
    print(df.shape)
    import matplotlib.pyplot as plt
    plt.scatter(df.age, df.BA,s=4,c='k')
    plt.plot([30,120],[30,120],c='r')
    plt.show()
    
    df_res = defaultdict(list)
    
    model2formula = {
    'AgeOnly':'inc_dementia ~ BAI + age',
    
    'Mediate':'inc_dementia ~ BAI + age + educcollege + BMI + sleepmed + walking + APOE4',
    
    'Full':   'inc_dementia ~ BAI + age + educcollege + BMI + sleepmed + walking + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI + APOE4',
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
    df_res.to_csv('BAI_results_SOF-binary.csv', index=False)
    
        

if __name__=='__main__':
    main()

