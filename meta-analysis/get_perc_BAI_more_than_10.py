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
    df = df[['BAI']].dropna()
    df *= 10  # it was /10 when saved
    print(f'{dataset}: {(df.BAI>=10).mean()*100:.1f}%')
    #print(df)
    print()
