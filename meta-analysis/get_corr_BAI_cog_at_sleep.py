import os
import numpy as np
from scipy.stats import *
import pandas as pd
from pingouin import partial_corr


base_dir = '/data/haoqisun/BAI_dementia_community'
datasets = ['MESA', 'ARIC', 'FHS', 'MrOS', 'SOF']

for dataset in datasets:
    if dataset == 'MESA':
        col = 'CASIscore'
        data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
        covar = ['age', 'sexM']
    elif dataset == 'MrOS':
        col = 'tms'
        data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
        covar = ['age']
    elif dataset == 'SOF':
        col = 'mmse'
        data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
        covar = ['age']
    elif dataset in ['FHS','ARIC']:
        col = 'mmse'
        data_path = os.path.join(base_dir, 'SHHS', dataset, f'dataset_{dataset}_table1.csv')
        covar = ['age', 'sexM']
    df = pd.read_csv(data_path)
    #val = df[['BAI',col]].dropna()
    #val = val.values.astype(float)
    #bai = val[:,0]
    #cog = val[:,1]
    #corr1 = pearsonr(bai, cog)
    #corr2 = spearmanr(bai, cog)
    #print(dataset, len(val), corr1, corr2)
    res1 = partial_corr(data=df, x='BAI', y=col, covar=covar, method='pearson')
    res2 = partial_corr(data=df, x='BAI', y=col, covar=covar, method='spearman')
    print(dataset)
    print(res1)
    print(res2)
    print()
