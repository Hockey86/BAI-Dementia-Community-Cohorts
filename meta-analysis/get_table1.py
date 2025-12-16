import os
import numpy as np
import pandas as pd


def rindex(alist, value):
    return len(alist) - alist[-1::-1].index(value) -1


base_dir = '/data/haoqisun/BAI_dementia_community'
datasets = ['MrOS', 'SOF', 'MESA', 'FHS', 'ARIC']
rows = [
'N',
'Age, mean (SD), y',
'Sex',
'    Female, n(%)',
'    Male, n(%)',
'Race & Ethnicity',
'    Asian, n(%)',
'    Black, n(%)',
'    Hispanic, n(%)',
'    White, n(%)',
'    Other, n(%)',
'Education',
'    <= High School, n(%)',
'    >= College, n(%)',
'APOE',
'    Genotype available, n(%)',
'    e4 carrier among available, n(%)',
'Incident dementia',
'    Incidence, n(%)',
'    Time from sleep recording, median (IQR), y',
'All-cause mortality',
'    Incidence, n(%)',
'    Time from sleep recording, median (IQR), y',
'Follow-up time, median (min-max), y',
'Body mass index (BMI), median (IQR), kg/m2',
'Baseline cognitive score, median (% normal)',
'Current smoker, n(%)',
'Hypertension, n(%)',
'Diabetes, n(%)',
'Myocardial infarction, n(%)',
'Stroke, n(%)',
'Depression, n(%)',
'Apnea-hypopnea index (AHI), median (IQR), /hour',
'Sleep medication use , n(%)',
'Sleep EEG-based brain age index (BAI), mean (SD), y',
]

df_res = pd.DataFrame(data={'Name':rows})
for d in datasets:
    df_res[d] = np.nan
    

# MrOS
dataset = 'MrOS'
di = datasets.index(dataset)+1

data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
df = pd.read_csv(data_path)
#df = df[df.prevalent_dementia==0].reset_index(drop=True)
#df = df[df.APOE4Count.notna()&df.BAI.notna()].reset_index(drop=True)
df_res.iloc[rows.index('N'), di] = len(df)
df_res.iloc[rows.index('Age, mean (SD), y'), di] = f'{df.age.mean():.1f} ({df.age.std():.1f})'
df_res.iloc[rows.index('    Female, n(%)'), di] = '0 (0%)'
df_res.iloc[rows.index('    Male, n(%)'), di] = f'{len(df)} (100%)'
print(dataset)
print(f'Female N(event) = {0}')
print(f'Male N(event) = {(df.event=="dementia").sum()}')
print(f'Age<70: N = {(df.age<70).sum()}, N(event) = {((df.age<70)&(df.event=="dementia")).sum()}')
print(f'Age>=70: N = {(df.age>=70).sum()}, N(event) = {((df.age>=70)&(df.event=="dementia")).sum()}')
df_res.iloc[rows.index('    Asian, n(%)'), di] = f'{df.race_Asian.sum()} ({df.race_Asian.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Black, n(%)'), di] = f'{df.race_Black.sum()} ({df.race_Black.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Hispanic, n(%)'), di] = '--'
df_res.iloc[rows.index('    White, n(%)'), di] = f'{len(df)-df.race_NonWhite.sum()} ({100-df.race_NonWhite.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Other, n(%)'), di] = f'{df.race_Other.sum()} ({df.race_Other.mean()*100:.1f}%)'
df_res.iloc[rows.index('    <= High School, n(%)'), di] = f'{(df.educcollege==0).sum()} ({(df.educcollege==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    >= College, n(%)'), di] = f'{df.educcollege.sum()} ({df.educcollege.mean()*100:.1f}%)'
n=(df.APOE4Count>0).sum()
d=df.APOE4Count.notna().sum()
df_res.iloc[rows.index('    Genotype available, n(%)'),di] = f'{d} ({df.APOE4Count.notna().mean()*100:.1f}%)'
df_res.iloc[rows.index('    e4 carrier among available, n(%)'),di] = f'{n} ({n/d*100:.1f}%)'
df_res.iloc[rows.index('    Incidence, n(%)'), di] = f'{(df.event=="dementia").sum()} ({(df.event=="dementia").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="dementia"], (25,50,75))
df_res.iloc[rows.index('    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rindex(rows, '    Incidence, n(%)'), di] = f'{(df.event=="death").sum()} ({(df.event=="death").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="death"], (25,50,75))
df_res.iloc[rindex(rows, '    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.time2event, (0,50,100))
df_res.iloc[rindex(rows, 'Follow-up time, median (min-max), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.BMI, (25,50,75))
df_res.iloc[rows.index('Body mass index (BMI), median (IQR), kg/m2'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df2 = df.dropna(subset='tms').reset_index(drop=True)
df_res.iloc[rows.index('Baseline cognitive score, median (% normal)'), di] = f'3MS {df2.tms.median():.0f} ({(df2.tms>=81).mean()*100:.1f}%)'
df_res.iloc[rows.index('Current smoker, n(%)'), di] = f'{df.smoke_current.sum()}({df.smoke_current.mean()*100:.1f}%)'
df_res.iloc[rows.index('Hypertension, n(%)'), di] = f'{df.hypertension.sum():.0f} ({df.hypertension.mean()*100:.1f}%)'
df_res.iloc[rows.index('Diabetes, n(%)'), di] = f'{df.diabetes.sum():.0f} ({df.diabetes.mean()*100:.1f}%)'
df_res.iloc[rows.index('Myocardial infarction, n(%)'), di] = f'{df.heartattack.sum():.0f} ({df.heartattack.mean()*100:.1f}%)'
df_res.iloc[rows.index('Stroke, n(%)'), di] = f'{df.stroke.sum():.0f} ({df.stroke.mean()*100:.1f}%)'
df_res.iloc[rows.index('Depression, n(%)'), di] = f'{df.depression.sum():.0f} ({df.depression.mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.AHI, (25,50,75))
df_res.iloc[rows.index('Apnea-hypopnea index (AHI), median (IQR), /hour'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rows.index('Sleep medication use , n(%)'), di] = f'{df.sleepmed.sum()}({df.sleepmed.mean()*100:.1f}%)'
df_res.iloc[rows.index('Sleep EEG-based brain age index (BAI), mean (SD), y'), di] = f'{df.BAI.mean():.1f} ({df.BAI.std():.1f})'




# SOF
dataset = 'SOF'
di = datasets.index(dataset)+1

data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
df = pd.read_csv(data_path)
#df = df[df.prevalent_dementia==0].reset_index(drop=True)
#df = df[df.APOE4.notna()&df.BAI.notna()].reset_index(drop=True)
df_res.iloc[rows.index('N'), di] = len(df)
df_res.iloc[rows.index('Age, mean (SD), y'), di] = f'{df.age.mean():.1f} ({df.age.std():.1f})'
df_res.iloc[rows.index('    Female, n(%)'), di] = f'{len(df)}(100%)'
df_res.iloc[rows.index('    Male, n(%)'), di] = '0(0%)'
print(dataset)
print(f'Female N(event) = {(df.event=="dementia").sum()}')
print(f'Male N(event) = {0}')
print(f'Age<70: N = {(df.age<70).sum()}, N(event) = {((df.age<70)&(df.event=="dementia")).sum()}')
print(f'Age>=70: N = {(df.age>=70).sum()}, N(event) = {((df.age>=70)&(df.event=="dementia")).sum()}')
df_res.iloc[rows.index('    Asian, n(%)'), di] = '0(0%)'
df_res.iloc[rows.index('    Black, n(%)'), di] = '0(0%)'
df_res.iloc[rows.index('    Hispanic, n(%)'), di] = '--'
df_res.iloc[rows.index('    White, n(%)'), di] = f'{(df.Race==1).sum()}({(df.Race==1).mean()*100:.1f}%)'
df_res.iloc[rows.index('    Other, n(%)'), di] = f'{(df.Race==2).sum()}({(df.Race==2).mean()*100:.1f}%)'
df_res.iloc[rows.index('    <= High School, n(%)'), di] = f'{(df.educcollege==0).sum()}({(df.educcollege==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    >= College, n(%)'), di] = f'{df.educcollege.sum()}({df.educcollege.mean()*100:.1f}%)'
n=(df.APOE4>0).sum()
d=df.APOE4.notna().sum()
df_res.iloc[rows.index('    Genotype available, n(%)'),di] = f'{d} ({df.APOE4.notna().mean()*100:.1f}%)'
df_res.iloc[rows.index('    e4 carrier among available, n(%)'),di] = f'{n} ({n/d*100:.1f}%)'
df_res.iloc[rows.index('    Incidence, n(%)'), di] = f'{(df.event=="dementia").sum()}({(df.event=="dementia").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="dementia"], (25,50,75))
df_res.iloc[rows.index('    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rindex(rows, '    Incidence, n(%)'), di] = f'{(df.event=="death").sum()}({(df.event=="death").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="death"], (25,50,75))
df_res.iloc[rindex(rows, '    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.time2event, (0,50,100))
df_res.iloc[rindex(rows, 'Follow-up time, median (min-max), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.BMI, (25,50,75))
df_res.iloc[rows.index('Body mass index (BMI), median (IQR), kg/m2'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df2 = df.dropna(subset='mmse').reset_index(drop=True)
df_res.iloc[rows.index('Baseline cognitive score, median (% normal)'), di] = f'MMSE {df2.mmse.median():.0f} ({(df2.mmse>=24).mean()*100:.1f}%)'
df_res.iloc[rows.index('Current smoker, n(%)'), di] = f'{df.smoke_current.sum()}({df.smoke_current.mean()*100:.1f}%)'
df_res.iloc[rows.index('Hypertension, n(%)'), di] = f'{df.hypertension.sum():.0f}({df.hypertension.mean()*100:.1f}%)'
df_res.iloc[rows.index('Diabetes, n(%)'), di] = f'{df.diabetes.sum():.0f}({df.diabetes.mean()*100:.1f}%)'
df_res.iloc[rows.index('Myocardial infarction, n(%)'), di] = f'{df.heartattack.sum():.0f}({df.heartattack.mean()*100:.1f}%)'
df_res.iloc[rows.index('Stroke, n(%)'), di] = f'{df.stroke.sum():.0f}({df.stroke.mean()*100:.1f}%)'
df_res.iloc[rows.index('Depression, n(%)'), di] = f'{df.depression.sum():.0f}({df.depression.mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.AHI, (25,50,75))
df_res.iloc[rows.index('Apnea-hypopnea index (AHI), median (IQR), /hour'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rows.index('Sleep medication use , n(%)'), di] = f'{df.sleepmed.sum()}({df.sleepmed.mean()*100:.1f}%)'
df_res.iloc[rows.index('Sleep EEG-based brain age index (BAI), mean (SD), y'), di] = f'{df.BAI.mean():.1f} ({df.BAI.std():.1f})'




# MESA
dataset = 'MESA'
di = datasets.index(dataset)+1

data_path = os.path.join(base_dir, dataset, f'dataset_{dataset}_table1.csv')
df = pd.read_csv(data_path)
#df = df[df.prevalent_dementia==0].reset_index(drop=True)
#df = df[df.APOE4Count.notna()&df.BAI.notna()].reset_index(drop=True)
df_res.iloc[rows.index('N'), di] = len(df)
df_res.iloc[rows.index('Age, mean (SD), y'), di] = f'{df.age.mean():.1f} ({df.age.std():.1f})'
df_res.iloc[rows.index('    Female, n(%)'), di] = f'{(df.sexM==0).sum()}({(df.sexM==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    Male, n(%)'), di] = f'{(df.sexM==1).sum()}({(df.sexM==1).mean()*100:.1f}%)'
print(dataset)
print(f'Female N(event) = {((df.event=="dementia")&(df.sexM==0)).sum()}')
print(f'Male N(event) = {((df.event=="dementia")&(df.sexM==1)).sum()}')
print(f'Age<70: N = {(df.age<70).sum()}, N(event) = {((df.age<70)&(df.event=="dementia")).sum()}')
print(f'Age>=70: N = {(df.age>=70).sum()}, N(event) = {((df.age>=70)&(df.event=="dementia")).sum()}')
df_res.iloc[rows.index('    Asian, n(%)'), di] = f'{df.race_Chinese.sum()}({df.race_Chinese.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Black, n(%)'), di] = f'{df.race_Black.sum()}({df.race_Black.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Hispanic, n(%)'), di] = f'{df.race_Hispanic.sum()}({df.race_Hispanic.mean()*100:.1f}%)'
df_res.iloc[rows.index('    White, n(%)'), di] = f'{df.race_White.sum()}({df.race_White.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Other, n(%)'), di] = '0(0%)'
df_res.iloc[rows.index('    <= High School, n(%)'), di] = f'{(df.educcollege==0).sum()}({(df.educcollege==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    >= College, n(%)'), di] = f'{df.educcollege.sum()}({df.educcollege.mean()*100:.1f}%)'
n=(df.APOE4Count>0).sum()
d=df.APOE4Count.notna().sum()
df_res.iloc[rows.index('    Genotype available, n(%)'),di] = f'{d} ({df.APOE4Count.notna().mean()*100:.1f}%)'
df_res.iloc[rows.index('    e4 carrier among available, n(%)'),di] = f'{n} ({n/d*100:.1f}%)'
df_res.iloc[rows.index('    Incidence, n(%)'), di] = f'{(df.event=="dementia").sum()}({(df.event=="dementia").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="dementia"], (25,50,75))
df_res.iloc[rows.index('    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rindex(rows, '    Incidence, n(%)'), di] = f'{(df.event=="death").sum()}({(df.event=="death").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="death"], (25,50,75))
df_res.iloc[rindex(rows, '    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.time2event, (0,50,100))
df_res.iloc[rindex(rows, 'Follow-up time, median (min-max), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.BMI, (25,50,75))
df_res.iloc[rows.index('Body mass index (BMI), median (IQR), kg/m2'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df2 = df.dropna(subset='CASIscore').reset_index(drop=True)
df_res.iloc[rows.index('Baseline cognitive score, median (% normal)'), di] = f'CASI {df2.CASIscore.median():.0f} ({(df2.CASIscore>=77).mean()*100:.1f}%)'
df_res.iloc[rows.index('Current smoker, n(%)'), di] = f'{df.smoke_current.sum()}({df.smoke_current.mean()*100:.1f}%)'
df_res.iloc[rows.index('Hypertension, n(%)'), di] = f'{df.hypertension.sum():.0f}({df.hypertension.mean()*100:.1f}%)'
df_res.iloc[rows.index('Diabetes, n(%)'), di] = f'{df.diabetes.sum():.0f}({df.diabetes.mean()*100:.1f}%)'
df_res.iloc[rows.index('Myocardial infarction, n(%)'), di] = f'{df.heartattack.sum():.0f}({df.heartattack.mean()*100:.1f}%)'
df_res.iloc[rows.index('Stroke, n(%)'), di] = f'{df.stroke.sum():.0f}({df.stroke.mean()*100:.1f}%)'
df_res.iloc[rows.index('Depression, n(%)'), di] = f'{df.depression.sum():.0f}({df.depression.mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.AHI, (25,50,75))
df_res.iloc[rows.index('Apnea-hypopnea index (AHI), median (IQR), /hour'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rows.index('Sleep medication use , n(%)'), di] = f'{df.sleepmed.sum()}({df.sleepmed.mean()*100:.1f}%)'
df_res.iloc[rows.index('Sleep EEG-based brain age index (BAI), mean (SD), y'), di] = f'{df.BAI.mean():.1f} ({df.BAI.std():.1f})'





# FHS
dataset = 'FHS'
di = datasets.index(dataset)+1
data_path = os.path.join(base_dir, 'SHHS', dataset, f'dataset_{dataset}_table1.csv')
df = pd.read_csv(data_path)
#df = df[df.prevalent_dementia==0].reset_index(drop=True)
#df = df[df.APOE4Count.notna()&df.BAI.notna()].reset_index(drop=True)
df_res.iloc[rows.index('N'), di] = len(df)
df_res.iloc[rows.index('Age, mean (SD), y'), di] = f'{df.age.mean():.1f} ({df.age.std():.1f})'
df_res.iloc[rows.index('    Female, n(%)'), di] = f'{(df.sexM==0).sum()}({(df.sexM==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    Male, n(%)'), di] = f'{(df.sexM==1).sum()}({(df.sexM==1).mean()*100:.1f}%)'
print(dataset)
print(f'Female N(event) = {((df.event=="dementia")&(df.sexM==0)).sum()}')
print(f'Male N(event) = {((df.event=="dementia")&(df.sexM==1)).sum()}')
print(f'Age<70: N = {(df.age<70).sum()}, N(event) = {((df.age<70)&(df.event=="dementia")).sum()}')
print(f'Age>=70: N = {(df.age>=70).sum()}, N(event) = {((df.age>=70)&(df.event=="dementia")).sum()}')
df_res.iloc[rows.index('    Asian, n(%)'), di] = f'{df.race_Asian.sum()}({df.race_Asian.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Black, n(%)'), di] = f'{df.race_Black.sum()}({df.race_Black.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Hispanic, n(%)'), di] = f'{df.ethnicity_Hispanic.sum()}({df.ethnicity_Hispanic.mean()*100:.1f}%)'
df_res.iloc[rows.index('    White, n(%)'), di] = f'{df.race_White.sum()}({df.race_White.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Other, n(%)'), di] = f'{df.race_Other.sum()}({df.race_Other.mean()*100:.1f}%)'
df_res.iloc[rows.index('    <= High School, n(%)'), di] = f'{(df.educcollege==0).sum()}({(df.educcollege==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    >= College, n(%)'), di] = f'{df.educcollege.sum()}({df.educcollege.mean()*100:.1f}%)'
n=(df.APOE4Count>0).sum()
d=df.APOE4Count.notna().sum()
df_res.iloc[rows.index('    Genotype available, n(%)'),di] = f'{d} ({df.APOE4Count.notna().mean()*100:.1f}%)'
df_res.iloc[rows.index('    e4 carrier among available, n(%)'),di] = f'{n} ({n/d*100:.1f}%)'
df_res.iloc[rows.index('    Incidence, n(%)'), di] = f'{(df.event=="dementia").sum()}({(df.event=="dementia").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="dementia"], (25,50,75))
df_res.iloc[rows.index('    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rindex(rows, '    Incidence, n(%)'), di] = f'{(df.event=="death").sum()}({(df.event=="death").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="death"], (25,50,75))
df_res.iloc[rindex(rows, '    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.time2event, (0,50,100))
df_res.iloc[rindex(rows, 'Follow-up time, median (min-max), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.BMI, (25,50,75))
df_res.iloc[rows.index('Body mass index (BMI), median (IQR), kg/m2'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df2 = df.dropna(subset='mmse').reset_index(drop=True)
df_res.iloc[rows.index('Baseline cognitive score, median (% normal)'), di] = f'MMSE {df2.mmse.median():.0f} ({(df2.mmse>=24).mean()*100:.1f}%)'
df_res.iloc[rows.index('Current smoker, n(%)'), di] = f'{df.smoke_current.sum()}({df.smoke_current.mean()*100:.1f}%)'
df_res.iloc[rows.index('Hypertension, n(%)'), di] = f'{df.hypertension.sum():.0f}({df.hypertension.mean()*100:.1f}%)'
df_res.iloc[rows.index('Diabetes, n(%)'), di] = f'{df.diabetes.sum():.0f}({df.diabetes.mean()*100:.1f}%)'
df_res.iloc[rows.index('Myocardial infarction, n(%)'), di] = f'{df.heartattack.sum():.0f}({df.heartattack.mean()*100:.1f}%)'
df_res.iloc[rows.index('Stroke, n(%)'), di] = f'{df.stroke.sum():.0f}({df.stroke.mean()*100:.1f}%)'
df_res.iloc[rows.index('Depression, n(%)'), di] = f'{df.depression.sum():.0f}({df.depression.mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.AHI, (25,50,75))
df_res.iloc[rows.index('Apnea-hypopnea index (AHI), median (IQR), /hour'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rows.index('Sleep medication use , n(%)'), di] = f'{df.sleepmed.sum()}({df.sleepmed.mean()*100:.1f}%)'
df_res.iloc[rows.index('Sleep EEG-based brain age index (BAI), mean (SD), y'), di] = f'{df.BAI.mean():.1f} ({df.BAI.std():.1f})'




# ARIC
dataset = 'ARIC'
di = datasets.index(dataset)+1
data_path = os.path.join(base_dir, 'SHHS', dataset, f'dataset_{dataset}_table1.csv')
df = pd.read_csv(data_path)
#df = df[df.prevalent_dementia==0].reset_index(drop=True)
#df = df[df.APOE4Count.notna()&df.BAI.notna()].reset_index(drop=True)
df = df[df.BAI.notna()].reset_index(drop=True)
df_res.iloc[rows.index('N'), di] = len(df)
df_res.iloc[rows.index('Age, mean (SD), y'), di] = f'{df.age.mean():.1f} ({df.age.std():.1f})'
df_res.iloc[rows.index('    Female, n(%)'), di] = f'{(df.sexM==0).sum()}({(df.sexM==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    Male, n(%)'), di] = f'{(df.sexM==1).sum()}({(df.sexM==1).mean()*100:.1f}%)'
print(dataset)
print(f'Female N(event) = {((df.event=="dementia")&(df.sexM==0)).sum()}')
print(f'Male N(event) = {((df.event=="dementia")&(df.sexM==1)).sum()}')
print(f'Age<70: N = {(df.age<70).sum()}, N(event) = {((df.age<70)&(df.event=="dementia")).sum()}')
print(f'Age>=70: N = {(df.age>=70).sum()}, N(event) = {((df.age>=70)&(df.event=="dementia")).sum()}')
df_res.iloc[rows.index('    Asian, n(%)'), di] = '0(0%)'
df_res.iloc[rows.index('    Black, n(%)'), di] = f'{df.race_Black.sum()}({df.race_Black.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Hispanic, n(%)'), di] = '--'
df_res.iloc[rows.index('    White, n(%)'), di] = f'{df.race_White.sum()}({df.race_White.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Other, n(%)'), di] = f'{df.race_Other.sum()}({df.race_Other.mean()*100:.1f}%)'
df_res.iloc[rows.index('    <= High School, n(%)'), di] = f'{(df.educcollege==0).sum()}({(df.educcollege==0).mean()*100:.1f}%)'
df_res.iloc[rows.index('    >= College, n(%)'), di] = f'{df.educcollege.sum()}({df.educcollege.mean()*100:.1f}%)'
df_res.iloc[rows.index('    Genotype available, n(%)'),di] = f'0 (0%)'
df_res.iloc[rows.index('    e4 carrier among available, n(%)'),di] = '--'
df_res.iloc[rows.index('    Incidence, n(%)'), di] = f'{(df.event=="dementia").sum()}({(df.event=="dementia").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="dementia"], (25,50,75))
df_res.iloc[rows.index('    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rindex(rows, '    Incidence, n(%)'), di] = f'{(df.event=="death").sum()}({(df.event=="death").mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.time2event[df.event=="death"], (25,50,75))
df_res.iloc[rindex(rows, '    Time from sleep recording, median (IQR), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.time2event, (0,50,100))
df_res.iloc[rindex(rows, 'Follow-up time, median (min-max), y'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
q1, q2, q3 = np.nanpercentile(df.BMI, (25,50,75))
df_res.iloc[rows.index('Body mass index (BMI), median (IQR), kg/m2'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df2 = df.dropna(subset='mmse').reset_index(drop=True)
df_res.iloc[rows.index('Baseline cognitive score, median (% normal)'), di] = f'MMSE {df2.mmse.median():.0f} ({(df2.mmse>=24).mean()*100:.1f}%)'
df_res.iloc[rows.index('Current smoker, n(%)'), di] = f'{df.smoke_current.sum()}({df.smoke_current.mean()*100:.1f}%)'
df_res.iloc[rows.index('Hypertension, n(%)'), di] = f'{df.hypertension.sum():.0f}({df.hypertension.mean()*100:.1f}%)'
df_res.iloc[rows.index('Diabetes, n(%)'), di] = f'{df.diabetes.sum():.0f}({df.diabetes.mean()*100:.1f}%)'
df_res.iloc[rows.index('Myocardial infarction, n(%)'), di] = f'{df.heartattack.sum():.0f}({df.heartattack.mean()*100:.1f}%)'
df_res.iloc[rows.index('Stroke, n(%)'), di] = f'{df.stroke.sum():.0f}({df.stroke.mean()*100:.1f}%)'
df_res.iloc[rows.index('Depression, n(%)'), di] = f'{df.depression.sum():.0f}({df.depression.mean()*100:.1f}%)'
q1, q2, q3 = np.nanpercentile(df.AHI, (25,50,75))
df_res.iloc[rows.index('Apnea-hypopnea index (AHI), median (IQR), /hour'), di] = f'{q2:.1f} ({q1:.1f}-{q3:.1f})'
df_res.iloc[rows.index('Sleep medication use , n(%)'), di] = f'{df.sleepmed.sum()}({df.sleepmed.mean()*100:.1f}%)'
df_res.iloc[rows.index('Sleep EEG-based brain age index (BAI), mean (SD), y'), di] = f'{df.BAI.mean():.1f} ({df.BAI.std():.1f})'


print(df_res)
breakpoint()
df_res.to_excel('table1.xlsx', index=False)


