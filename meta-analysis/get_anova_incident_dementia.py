import os
import numpy as np
from scipy.stats import *
import pandas as pd


base_dir = '/data/haoqisun/BAI_dementia_community'
datasets = ['SOF', 'MESA', 'ARIC', 'FHS', 'MrOS']

vals = []
for dataset in datasets:
    if dataset == 'MESA':
        data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
    elif dataset == 'MrOS':
        data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
    elif dataset == 'SOF':
        data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
    elif dataset in ['FHS','ARIC']:
        data_path = os.path.join(base_dir, 'SHHS', dataset, f'dataset_{dataset}_table1.csv')
    df = pd.read_csv(data_path)
    breakpoint()
    vals.append([(df.event=='dementia').sum(), (df.event!='dementia').sum()])
vals = np.array(vals)
res = chi2_contingency(vals)
print(res)
