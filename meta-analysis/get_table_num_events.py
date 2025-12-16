import numpy as np
import pandas as pd


datasets = ['MESA', 'ARIC', 'FHS-OS', 'MrOS', 'SOF']
paths = ['../MESA/dataset_MESA_table1.csv',
         '../SHHS/ARIC/dataset_ARIC_table1.csv',
         '../SHHS/FHS/dataset_FHS_table1.csv',
        '../MrOS/dataset_MrOS_table1.csv',
        '../SOF/dataset_SOF_table1.csv',]

times = np.arange(1,25+1,1)
res = {'time':times}
for dataset, path in zip(datasets, paths):
    df = pd.read_csv(path)
    res[dataset+'_dementia'] = []
    res[dataset+'_death'] = []
    for time in times:
        e = df.event[(df.time2event<=time)&(df.time2event>0)]
        count1 = (e=='dementia').sum()
        count2 = (e=='death').sum()
        res[dataset+'_dementia'].append(count1)
        res[dataset+'_death'].append(count2)
df_res = pd.DataFrame(data=res)
print(df_res)
df_res['total_dementia'] = df_res[[f'{d}_dementia' for d in datasets]].values.sum(axis=1)
df_res['total_death'] = df_res[[f'{d}_death' for d in datasets]].values.sum(axis=1)

for d in datasets+['total']:
    df_res[d] = df_res[f'{d}_dementia'].astype(str)+', '+df_res[f'{d}_death'].astype(str)
    df_res = df_res.drop(columns=[f'{d}_dementia', f'{d}_death'])
print(df_res)
df_res.to_excel('table_num_events.xlsx', index=False)


