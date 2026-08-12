"""
v2: use 3MS decline>=1.5std (as in the paper) OR alz medication, use variables from the .dta file if possible
v3: available from .dta (as in the paper)
v4: available from .dta (as in the paper)
vs2: create it, as in the paper (same as v2)

* as in the paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6699896/): 
report of physician-diagnosed dementia: MHALZH

use of dementia medication (verified by clinic staff based on examination of pill bottles): M1ALZH

having a change in 3MS scores ≥ 1.5 standard deviation worse than the mean change from baseline to any follow-up visit (equal to a decline of 7.32, 9.43 and 13.62 points on 3MS from baseline to visit 2, visit 3 and visit 4, respectively)
"""
import numpy as np
import pandas as pd


def main():
    #"""
    df = pd.read_stata('final_Dataset_psg_AD_230616.dta')
    df = df.rename(columns={x:x.replace('3ms','tms') for x in df.columns})  # to use statsmodels.formula
    df['id'] = df.id.str.lower()

    df['race_White'] = 0; df.loc[df.girace==1, 'race_White'] = 1
    df['race_Black'] = 0; df.loc[df.girace==2, 'race_Black'] = 1
    df['race_Asian'] = 0; df.loc[df.girace==3, 'race_Asian'] = 1
    df['race_Other'] = 0; df.loc[df.girace>=4, 'race_Other'] = 1
    df.loc[df.pqpslmed<=1, 'pqpslmed'] = 0
    df.loc[df.pqpslmed>=2, 'pqpslmed'] = 1
    df.loc[df.gieduc<=4, 'gieduc'] = 0
    df.loc[df.gieduc>=5, 'gieduc'] = 1
    
    #df2 = pd.read_csv('mastersheet_MrOS.csv')
    #df2 = df2.rename(columns={'SID':'id'})
    #df = df.merge(df2, on='id', how='left', validate='1:1')

    #df2 = pd.read_csv('/data/haoqisun/BAI_dementia_community/extract_brain_age/brain_age_using_luna/BAI_results.csv')
    df2 = pd.read_csv('BAI_MrOS_with_features.csv')
    df2 = df2.rename(columns={'SID':'id'})
    df = df.merge(df2.drop(columns='Age'), on='id', how='left', validate='1:1')
    
    sites = ['sd', 'pi', 'pa', 'bi', 'po', 'mn']
    for site in sites:
        df[f'site_{site}'] = df.id.str.startswith(site).astype(int)

    cols = ['id', 'prevalent_dementia', 'inc_dementia', 'vsage1',
       'actnap2p', 'aclsepmp', 'acsminmp', 'acslatmp', 'acwasomp',
       'acsefnmp', 'acnap5mp', 'gieduc', 'race_White', 'race_Black', 'race_Asian', 'race_Other',
       'v1park', 'posllatp',
       'popcsa90', 'powaso', 'potmst1p', 'potmst2p', 'potms34p', 'potmremp',
       'poai_all', 'poordi3', 'pocai4p', 'poslpeff', 'poordi4', 'poavgplm',
       'poavplma', 'vscpap', 'pqpslmed', 'slsa', 'mhdiab', 'vspark', 'mhparkt',
       'mhmi', 'mhstrk', 'mhbp', 'tusmknow', 'tusmkcgn', 'slnap', 'slnaphr',
       'vstbs', 'dpgdsyn', 'slnaphwk', 'epepwort', 'epeds',
       'qlfxst51', 'vstms', 'hwbmi', 'pascore', 'dpgds15', 'tursmoke',
       'vs2park', 'v2tbs', 'v2park', 'v2age1', 'v2tms', 'v2vsfytm', 'v3park',
       'v3dement', 'v3dementt', 'v3tbs', 'v3age1', 'v3tms', 'v3vsfytm',
       'v4park', 'v4dement', 'v4dementt', 'v4tbs', 'v4age1', 'v4tms',
       'tmms4sc', 'v4sfytm', 'efstatus', 'dadead', 'fucytime', 'efvsstat',
       'efv2stat', 'efv3stat', 'efvs2sta', 'efv4stat', 'fuv2yt', 'fuvsyt',
       'fuv3yt', 'fuv4yt', 'v1alzh', 'v2alzh', 'v3alzh', 'v4alzh', 'vsbenzo',
       'vsalzh', 'vsslpmed', 'acnap10mp', 'dement3',
       'tmms3sc', 'tmms2sc', 'dement4', 'BA', 'BAI',
       'COUPL_OVERLAP_C', 'DENS_C',
       'alpha_bandpower_kurtosis_C_N2', 'alpha_bandpower_mean_C_N1',
       'delta_alpha_mean_C_N3', 'delta_bandpower_kurtosis_C_N2',
       'delta_bandpower_mean_C_N3', 'delta_theta_mean_C_N3', 'kurtosis_N2_C',
       'kurtosis_N3_C', 'sigma_bandpower_kurtosis_C_N2',
       'theta_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N3',]
    
    df = df[cols]
    breakpoint()
    df = df.rename(columns={'vsage1':'age', 'gieduc':'educcollege', 'hwbmi':'BMI', 'pqpslmed':'sleepmed', 'mhdiab':'diabetes', 'mhbp':'hypertension', 'mhmi':'heartattack', 'mhstrk':'stroke', 'dpgdsyn':'depression', 'vstms':'tms', 'poordi4':'AHI'})
    #"""
    df = pd.read_csv('dataset_MrOS_.csv')
    df['race_NonWhite'] = ((df.race_Black+df.race_Asian+df.race_Other)>0).astype(int)
    df['smoke_current'] = (df.tursmoke==2).astype(int)
    cols = ['id', 'prevalent_dementia', 'inc_dementia', 'dement3', 'dement4', 'BAI',
       'COUPL_OVERLAP_C', 'DENS_C',
       'alpha_bandpower_kurtosis_C_N2', 'alpha_bandpower_mean_C_N1',
       'delta_alpha_mean_C_N3', 'delta_bandpower_kurtosis_C_N2',
       'delta_bandpower_mean_C_N3', 'delta_theta_mean_C_N3', 'kurtosis_N2_C',
       'kurtosis_N3_C', 'sigma_bandpower_kurtosis_C_N2',
       'theta_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N3',
    'age', 'educcollege', 'BMI', 'race_Black', 'race_Asian', 'race_Other', 'race_NonWhite',
    'sleepmed', 'pascore', 'diabetes', 'hypertension', 'heartattack', 'stroke', 'depression', 'smoke_current',
    'tms', 'AHI']#, 'v2alzh']
    df = df[cols]
    df = df[df.BAI.notna()].reset_index(drop=True)
    
    # get APOE e4
    df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/GNAUG15.SAS7BDAT')
    df2['APOE4Count'] = df2.rs429358.astype(str).str.count('C').astype(float)
    df2.loc[df2.rs429358==b'.M', 'APOE4Count'] = np.nan
    df2['id'] = df2.ID.astype(str).str.lower()
    df = df.merge(df2[['id', 'APOE4Count']], on='id', how='inner', validate='1:1')
    
    """
    # add time to death
    #df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/effeb24.sas7bdat')
    df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/efaug18.sas7bdat')
    df2['id'] = df2.ID.astype(str).str.lower()
    df2 = df2[df2.FUCYTIME>0].reset_index(drop=True)  # remove subjects without follow up
    df = df.merge(df2[['id', 'EFSTATUS', 'FUCYTIME']], on='id', how='inner', validate='1:1')
    
    #qual to a decline of 7.32, 9.43 and 13.62 points on 3MS from baseline to visit 2, visit 3 and visit 4, respect
    
    # add time to dementia
    df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/v2feb24.sas7bdat')
    df2m = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/m2feb21.sas7bdat') # M1ALZH here is the same as v2alzh in .dta
    df2['id'] = df2.ID.astype(str).str.lower()
    df2m['id'] = df2m.ID.astype(str).str.lower()
    df2 = df2.merge(df2m, on='id', how='inner', validate='1:1')
    df2['TMM12SCR'] = (df2.TMM12SCR<=-7.32).astype(float)
    df2.loc[df2.TMM12SCR.isna(), 'TMM12SCR'] = np.nan
    df2['dement2'] = (df2[['TMM12SCR', 'M1ALZH']].sum(axis=1)>0).astype(float)
    df2.loc[df2.TMM12SCR.isna()&df2.M1ALZH.isna(), 'dement2'] = np.nan
    
    df3 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/v3feb24.sas7bdat')
    df3['id'] = df3.ID.astype(str).str.lower()
    #df3['dement3'] = ((df3[['MHALZH', 'MHALZHT']].sum(axis=1)>0)|(df3.TMM13SCR<-9.43)).astype(float)
    #df3.loc[df3[['MHALZH', 'MHALZHT']].isna().all(axis=1), 'dement3'] = np.nan
    
    df4 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/v4feb24.sas7bdat')
    df4['id'] = df4.ID.astype(str).str.lower()
    #df4['dement4'] = ((df4[['MHALZH', 'MHALZHT']].sum(axis=1)>0)|(df4.TMM14SC<-13.62)).astype(float)
    #df4.loc[df4[['MHALZH', 'MHALZHT']].isna().all(axis=1), 'dement4'] = np.nan

    #df.merge(dfvs2,on='id',how='left',validate='1:1').TMM1S2SC.std()*1.5
    #8.62
    dfvs2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/vs2feb24.sas7bdat')
    dfvs2m = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/MS2AUG16.SAS7BDAT')
    dfvs2['id'] = dfvs2.ID.astype(str).str.lower()
    dfvs2m['id'] = dfvs2m.ID.astype(str).str.lower()
    dfvs2 = dfvs2.merge(dfvs2m, on='id', how='inner', validate='1:1')
    dfvs2['TMM1S2SC'] = (dfvs2.TMM1S2SC<=-8.62).astype(float)
    dfvs2.loc[dfvs2.TMM1S2SC.isna(), 'TMM1S2SC'] = np.nan
    dfvs2['dementvs2'] = (dfvs2[['TMM1S2SC', 'M1ALZH', 'MHALZH']].sum(axis=1)>0).astype(float)
    dfvs2.loc[dfvs2.TMM1S2SC.isna()&dfvs2.M1ALZH.isna()&dfvs2.MHALZH.isna(), 'dementvs2'] = np.nan
    
    df = df.merge(df2[['id', 'dement2', 'V2V1FYTM']], on='id', how='left', validate='1:1')
    df = df.merge(df3[['id', 'V3V1FYTM']], on='id', how='left', validate='1:1')
    df = df.merge(dfvs2[['id', 'dementvs2', 'VS21FYTM']], on='id', how='left', validate='1:1')
    df = df.merge(df4[['id', 'V41FYTM']], on='id', how='left', validate='1:1')
    
    cols1 = ['dement2', 'dement3', 'dementvs2', 'dement4']
    cols2 = ['V2V1FYTM', 'V3V1FYTM', 'VS21FYTM', 'V41FYTM']
    df['event'] = ''
    df['time2event'] = np.nan
    for i in range(len(df)):
        if df.loc[i, cols1].sum()>0:
            df.loc[i,'event'] = 'dementia'
            df.loc[i,'time2event'] = df.loc[i, cols2[np.where(df.loc[i, cols1]==1)[0][0]]]
        elif df.EFSTATUS.iloc[i]==1:
            df.loc[i,'event'] = 'death'
            df.loc[i,'time2event'] = df.FUCYTIME.iloc[i]
        else:
            df.loc[i,'event'] = 'censor'
            df.loc[i,'time2event'] = df.FUCYTIME.iloc[i]
    
    # add time to PSG
    df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/POSFEB23.SAS7BDAT')
    df2['id'] = df2.ID.astype(str).str.lower()
    df = df.merge(df2[['id', 'PODAYSVS']], on='id', how='inner', validate='1:1')
    df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/vsfeb24.sas7bdat')
    df2['id'] = df2.ID.astype(str).str.lower()
    df = df.merge(df2[['id', 'VSV1FYTM']], on='id', how='inner', validate='1:1')
    df['time2event'] = df.time2event - (df.VSV1FYTM+df.PODAYSVS/365.25)
    """
    df2 = pd.read_csv('mros_general_data_w_pd_and_dem.csv')
    df2['id'] = df2['id'].astype(str).str.lower()
    df2 = df2.rename(columns={'inc_dem_or_death':'event', 'follow_up_time_dem_with_death':'time2event'})
    df2['time2event'] = df2['time2event']/365
    df2.loc[df2.event==0, 'event'] = 'censor'
    df2.loc[df2.event==1, 'event'] = 'dementia'
    df2.loc[df2.event==2, 'event'] = 'death'
    df = df.merge(df2[['id', 'event', 'time2event']], on='id', how='inner', validate='1:1')
    
    df['inc_dementia'] = ((df.event=='dementia')&(df.time2event>0)&(df.prevalent_dementia==0)).astype(int)
    df.loc[df.prevalent_dementia==1, 'inc_dementia'] = 0
    
    cols = ['id', 'prevalent_dementia', 'inc_dementia', 'BAI',
       'COUPL_OVERLAP_C', 'DENS_C',
       'alpha_bandpower_kurtosis_C_N2', 'alpha_bandpower_mean_C_N1',
       'delta_alpha_mean_C_N3', 'delta_bandpower_kurtosis_C_N2',
       'delta_bandpower_mean_C_N3', 'delta_theta_mean_C_N3', 'kurtosis_N2_C',
       'kurtosis_N3_C', 'sigma_bandpower_kurtosis_C_N2',
       'theta_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N3',
    'age', 'educcollege', 'BMI', 'race_Black', 'race_Asian', 'race_Other', 'race_NonWhite',
    'sleepmed', 'pascore', 'diabetes', 'hypertension', 'heartattack', 'stroke', 'depression', 'smoke_current',
    'tms', 'AHI', 'APOE4Count', 'event', 'time2event']
    df = df[cols]
    
    print(df)
    breakpoint()
    df.to_csv('dataset_MrOS.csv', index=False)


if __name__=='__main__':
    main()

