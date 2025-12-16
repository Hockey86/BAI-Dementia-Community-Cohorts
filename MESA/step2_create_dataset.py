import os
import numpy as np
import pandas as pd


def main():
    dataset_dir = '/data/haoqisun/dataset_MESA_2024'
    
    df_ex1 = pd.read_stata(os.path.join(dataset_dir, 'MESAE1FinalLabel20240307.dta'), convert_categoricals=False)
    df_ex2 = pd.read_stata(os.path.join(dataset_dir, 'MESAe2FinalLabel20240209.dta'))
    df_ex3 = pd.read_stata(os.path.join(dataset_dir, 'MESAe3FinalLabel20240212.dta'))
    df_ex4 = pd.read_stata(os.path.join(dataset_dir, 'MESAe4FinalLabel20240307.dta'))
    df_ex5 = pd.read_stata(os.path.join(dataset_dir, 'MESAe5_FinalLabel_20240213.dta'))
    df_ex6 = pd.read_sas(os.path.join(dataset_dir, 'MESAe6_FinalLabel_20220125.sas7bdat'))
    
    # get CASI
    df_ex5_casi = pd.read_sas(os.path.join(dataset_dir, 'MESAe5_CASI_20190823.sas7bdat'))
    df_ex5_casi.loc[(df_ex5_casi.valid5!=1)|(df_ex5_casi.flagcasi5c==1), 'casisum5c'] = np.nan
    df_ex5 = df_ex5.drop(columns='casisum5c').merge(df_ex5_casi[['idno', 'casisum5c']], on='idno', how='left', validate='1:1')
    
    df_ex6_casi = pd.read_stata(os.path.join(dataset_dir, 'MESAe6_CASI_20200721.dta'))
    df_ex6_casi.loc[(df_ex6_casi.valid6!='Valid')|(df_ex6_casi.flagcasi6c=='Yes'), 'casisum6c'] = np.nan
    df_ex6 = df_ex6.merge(df_ex6_casi[['idno', 'casisum6c']], on='idno', how='left', validate='1:1')
    
    # BAI
    df = pd.read_csv('BAI_MESA-with-features.csv')
    df = df.dropna(subset='BAI', ignore_index=True)
    dfm = pd.read_csv(os.path.join(dataset_dir, 'mesa_nsrr_bridge_ids.csv'))
    df = df.merge(dfm, on='mesaid', how='inner', validate='1:1')

    # outcome
    
    # use meds alz1c--alz6c per guideline
    df_alz = df_ex1[['idno','alzh1c']].merge(df_ex2[['idno','alzh2c']], on='idno', how='left', validate='1:1')
    df_alz = df_alz.merge(df_ex3[['idno','alzh3c']], on='idno', how='left', validate='1:1')
    df_alz = df_alz.merge(df_ex4[['idno','alzh4c']], on='idno', how='left', validate='1:1')
    df_alz = df_alz.merge(df_ex5[['idno','alzh5c']], on='idno', how='left', validate='1:1')
    df_alz = df_alz.merge(df_ex6[['idno','alzh6c']], on='idno', how='left', validate='1:1')
    df_alz['alzh2c'] = df_alz.alzh2c.map({'No':0, 'Yes':1}).astype(float)
    df_alz['alzh3c'] = df_alz.alzh3c.map({'0: NO':0, '1: YES':1}).astype(float)
    df_alz['alzh4c'] = df_alz.alzh4c.map({'0: NO':0, '1: YES':1}).astype(float)
    df_alz['alzh5c'] = df_alz.alzh5c.map({'0: NO':0, '1: YES':1}).astype(float)
    df = df.merge(df_alz, on='idno', how='left', validate='1:1')
    
    df2 = pd.read_stata(os.path.join(dataset_dir, 'EventsDementiaThru2018_20210826.dta'))
    df2 = df2.merge(df_ex6[['idno', 'e16dyc']], on='idno', how='left', validate='1:1')
    df2['days_bl_to_dementia'] = np.nan
    df2['days_bl_to_death'] = np.nan
    df2['days_bl_to_ltfu'] = np.nan
    ids = df2.icddementia_dh=='1: Yes'
    df2.loc[ids, 'days_bl_to_dementia'] = df2.icddementia_dhtt[ids].values
    ids = (df2.icddementia_dh=='0: No')&(df2.tot_d=='1: Yes')
    df2.loc[ids, 'days_bl_to_death'] = df2.dthcode_dtt[ids].values
    ids = (df2.icddementia_dh=='0: No')&(df2.tot_d=='0: No')
    df2.loc[ids, 'days_bl_to_ltfu'] = df2.icddementia_dhtt[ids].values
    
    # add dementia based on change in CASI from ex5 to ex6 >=1.5std
    df_ex56 = df_ex5.merge(df_ex6, on='idno', how='inner', validate='1:1')
    df_ex56['casi56change'] = df_ex56.casisum6c-df_ex56.casisum5c
    cutoff = df_ex56.casi56change.std()*1.5  # 9.81
    sids_dementia_cog = df_ex56.idno[df_ex56.casi56change <= -cutoff]
    ids = df2.days_bl_to_dementia.isna()&np.in1d(df2.idno, sids_dementia_cog)
    df2.loc[ids, 'days_bl_to_dementia'] = df2.e16dyc[ids]
    df2.loc[ids, 'days_bl_to_death'] = np.nan
    df2.loc[ids, 'days_bl_to_ltfu'] = np.nan
    
    # days of sleep visit since BL
    df3 = pd.read_stata(os.path.join(dataset_dir, 'MESAe5_SleepPolysomn_20160922.dta'))  # sleep visit is around ex5
    df3 = df3[(df3.havepsg5=='1: YES')&(df3.status_psg5=='1: PASSED')&(np.in1d(df3.quo2m15, ['5: SIGNAL GOOD >= 95% OF SLEEP TIME', '4: SIGNAL GOOD FOR 75-94% OF SLEEP TIME', '3:  SIGNAL GOOD FOR 50-74% OF SLEEP TIME']))&(df3.slewake5=='0: NO')].reset_index(drop=True)
    df3 = df3.merge(df_ex5[['idno', 'e15dyc']], on='idno', how='inner', validate='1:1')
    df2 = df2.merge(df3, on='idno', how='inner', validate='1:1')
    df2['days_bl_to_sleepstudy'] = df2.e15dyc+df2.stdypdy5c
    df2['days_sl_to_dementia'] = df2.days_bl_to_dementia-df2.days_bl_to_sleepstudy
    df2['days_sl_to_death'] = df2.days_bl_to_death-df2.days_bl_to_sleepstudy
    df2['days_sl_to_ltfu'] = df2.days_bl_to_ltfu-df2.days_bl_to_sleepstudy

    df = df.merge(df2[['idno', 'days_bl_to_sleepstudy', 'days_sl_to_dementia', 'days_sl_to_death', 'days_sl_to_ltfu']], on='idno', how='inner', validate='1:1')
    """
# compare age and sex of missing BAI vs not
# first comment line 27: #df = df.dropna(subset='BAI', ignore_index=True)
from scipy.stats import ttest_ind

df3=df2.merge(df,on='idno',how='left',validate='1:1')

aa = df3.age[df3.BA.notna()].dropna().values
bb = df3.age[df3.BA.isna()].dropna().values
print(aa.mean(), aa.std(), bb.mean(), bb.std())
# 69.45425188374597 9.100403101178479 70.5813953488372 10.234989140402856
print(ttest_ind(aa, bb))
# TtestResult(statistic=-0.8001267697892754, pvalue=0.4237374321000379, df=1899.0)

aa=df3.sexM[df3.BA.notna()].dropna().values
bb=df3.sexM[df3.BA.isna()].dropna().values
print(aa.mean()*100, bb.mean()*100)
# 47.0% 41.9%
print(proportions_ztest([aa.sum(),bb.sum()],[len(aa),len(bb)])[1])
# 0.5054952810669111
    """
    df['inc_dementia'] = ((df.days_sl_to_dementia>0)|(df.alzh6c==1)).astype(int)
    df['prevalent_dementia'] = ((df.days_sl_to_dementia<=0)|(df[['alzh1c', 'alzh2c', 'alzh3c', 'alzh4c', 'alzh5c']]==1).any(axis=1)).astype(int)
    
    # convert days_sl_to* to time2event and event
    df['time2event'] = df[['days_sl_to_dementia', 'days_sl_to_death', 'days_sl_to_ltfu']].sum(axis=1)/365.25
    df['event'] = ''
    df.loc[df.days_sl_to_dementia.notna(), 'event'] = 'dementia'
    df.loc[df.days_sl_to_death.notna(), 'event'] = 'death'
    df.loc[df.days_sl_to_ltfu.notna(), 'event'] = 'censor'
    df = df.drop(columns=['days_sl_to_dementia', 'days_sl_to_death', 'days_sl_to_ltfu'])
    
    # education and race
    
    df_ex1['educcollege'] = (df_ex1.educ1>=7).astype(float)
    df_ex1.loc[df_ex1.educ1.isna(), 'educcollege'] = np.nan
    df_ex1['race_White'] = (df_ex1.race1c==1).astype(int)
    df_ex1['race_Chinese'] = (df_ex1.race1c==2).astype(int)
    df_ex1['race_Black'] = (df_ex1.race1c==3).astype(int)
    df_ex1['race_Hispanic'] = (df_ex1.race1c==4).astype(int)
    df = df.merge(df_ex1[['idno', 'educcollege', 'race_White', 'race_Chinese', 'race_Black', 'race_Hispanic']], on='idno', how='left', validate='1:1')
    
    # BMI and physical activity, diabetes, hypertension, depression, cognitive score
    
    df_ex5 = df_ex5.rename(columns={'bmi5c':'BMI', 'walkmn5c':'walk_min_per_wk', 'casisum5c':'CASIscore'})
    df_ex5['diabetes'] = df_ex5.dm035c.map({
        'NORMAL':0, 'IMPAIRED FASTING GLUCOSE':1,
        'UNTREATED DIABETES':2, 'TREATED DIABETES':3}).astype(float)
    df_ex5['diabetes'] = (df_ex5.diabetes>=2).astype(float)
    df_ex5.loc[df_ex5.dm035c.isna(), 'diabetes'] = np.nan
    df_ex5['hypertension'] = df_ex5.htn5c.map({
        '0: NO':0, '1: YES':1}).astype(float)
    df_ex5.loc[df_ex5.htn5c.isna(), 'hypertension'] = np.nan
    # cut-off=16: https://www.apa.org/pi/about/publications/caregivers/practice-settings/assessment/tools/depression-scale
    df_ex5['depression'] = (df_ex5.cesd5c>=16).astype(float)
    df_ex5.loc[df_ex5.cesd5c.isna(), 'depression'] = np.nan
    df_ex5.loc[df_ex5.CASIscore<20, 'CASIscore'] = np.nan
    df_ex5['smoke_current'] = (df_ex5.smkstat5=='Current smoker').astype(int)
    df = df.merge(df_ex5[['idno', 'BMI', 'walk_min_per_wk', 'diabetes', 'hypertension', 'depression', 'smoke_current', 'CASIscore']], on='idno', how='left', validate='1:1')
    
    #df2 = pd.read_sas(os.path.join(dataset_dir, 'MESAe5_CASI_20190823.sas7bdat'))
    
    # heart attack, stroke
    
    df2 = pd.read_stata(os.path.join(dataset_dir, 'MESAEvThru2019_20220727.dta'))
    df2['mi'] = df2.mi.map({'No':0, 'Yes':1}).astype(float)
    df2['strk'] = df2.strk.map({'No':0, 'Yes':1}).astype(float)
    df = df.merge(df2[['idno', 'mi', 'mitt', 'strk', 'strktt', 'prebase', 'fuptt']], on='idno', how='left', validate='1:1')
    df['heartattack'] = ((df.mi==1)&(df.mitt<=df.days_bl_to_sleepstudy)).astype(float)
    df.loc[df.mi.isna()|df.mitt.isna(), 'heartattack'] = np.nan
    df.loc[df.prebase=='mi', 'heartattack'] = 1
    df['stroke'] = ((df.strk==1)&(df.strktt<=df.days_bl_to_sleepstudy)).astype(float)
    df.loc[df.strk.isna()|df.strktt.isna(), 'stroke'] = np.nan
    
    # meds
    
    df2 = pd.read_stata(os.path.join(dataset_dir, 'MESAe5_SleepQ_20140617.dta'))
    df2['sleepmed'] = df2.slpngpills5.map({
        '1: NO, NOT IN THE PAST 4 WEEKS':1,
        '2: YES, LESS THAN ONCE A WEEK':2,
        '3: YES, 1 OR 2 TIMES A WEEK':3,
        '4: YES, 3 OR 4 TIMES A WEEK':4,
        '5: YES, 5 OR MORE TIMES A WEEK':5 }).astype(float)
    df2['sleepmed'] = (df2.sleepmed>=3).astype(float)
    df2.loc[df2.slpngpills5.isna(), 'sleepmed'] = np.nan
    df = df.merge(df2[['idno', 'sleepmed']], on='idno', how='left', validate='1:1')
    
    # AHI
    df_nsrr = pd.read_csv('mesa-sleep-dataset-haoqi.csv')
    df_nsrr = df_nsrr.rename(columns={'ahi_a0h4':'AHI'})
    df = df.merge(df_nsrr[['idno', 'AHI']], on='idno', how='left', validate='1:1')
    
    # APOE
    df2 = pd.read_sas(os.path.join(dataset_dir, 'MESA_ApoE_03102014.sas7bdat'))
    df2['APOE4Count'] = df2.ApoE.map({33:0, 34:1, 23:0, 44:2, 24:1, 22:0})
    df = df.merge(df2[['idno', 'APOE4Count']], on='idno', how='left', validate='1:1')
    breakpoint()
    
    df.to_csv('dataset_MESA.csv', index=False)
    
    

if __name__=='__main__':
    main()

