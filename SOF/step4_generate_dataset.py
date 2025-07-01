import os
import numpy as np
import pandas as pd


def main():
    """
     ~ BAI + age + educcollege + BMI + sleepmed + walking + diabetes + hypertension + heartattack + stroke + depression + mmse + AHI + APOE4
    """
    
    # age and BAI
    df = pd.read_excel('BAI_SOF.xlsx')
    df = df[df.BAI.notna()].reset_index(drop=True)
    df = df[df.SID<10000].reset_index(drop=True)  # since AA do not have mortality info
    df = df.rename(columns={'SID':'ID', 'Age':'age'})
    
    data_dir = '/data/haoqisun/dataset_MrOS_SOF_PSG/SOF-datasets'
    
    # sleep, v8, and v9 days
    #dfv1 = pd.read_sas(os.path.join(data_dir, 'v1demogr.sas7bdat'))  # (V1DAYS==0).all() == True
    dfv8 = pd.read_sas(os.path.join(data_dir, 'v8demogr.sas7bdat'))
    dfv9 = pd.read_sas(os.path.join(data_dir, 'v9demogr.sas7bdat'))
    dfsl = pd.read_sas(os.path.join(data_dir, 'v8sleep.sas7bdat'))
    dfsl = dfsl[dfsl.V8PDAYS.notna()].reset_index(drop=True)
    #df = df.merge(dfv1[['ID', 'V1DAYS']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv8[['ID', 'V8DAYS']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv9[['ID', 'V9DAYS']], on='ID', how='left', validate='1:1')
    df = df.merge(dfsl, on='ID', how='left', validate='1:1')
    df = df[(df.V8QUEEG1>=2)&(df.V8QUEEG2>=2)].reset_index(drop=True)
    df['SLDAYS'] = df.V8DAYS+df.V8PDAYS
    df['SL_V9_DAYS'] = df.V9DAYS-df.SLDAYS
    df['SL_V8_DAYS'] = df.V8DAYS-df.SLDAYS
    
    # outcome -- death
    dfv1 = pd.read_sas(os.path.join(data_dir, 'v1endpt.sas7bdat'))
    dfv1 = dfv1.rename(columns={'V1DEATH':'DEATH', 'V1FOLALL':'FOLALL'})
    df = df.merge(dfv1[['ID', 'DEATH', 'FOLALL']], on='ID', how='left', validate='1:1')
    df['FOLALL'] = df.FOLALL-df.SLDAYS
    
    # outcome -- dementia
    # diagnosis
    dfv4 = pd.read_sas(os.path.join(data_dir, 'v4medhx.sas7bdat'))
    dfv5 = pd.read_sas(os.path.join(data_dir, 'v5medhx.sas7bdat'))
    dfv6 = pd.read_sas(os.path.join(data_dir, 'v6medhx.sas7bdat'))
    dfv8 = pd.read_sas(os.path.join(data_dir, 'v8medhx.sas7bdat'))
    dfv9 = pd.read_sas(os.path.join(data_dir, 'v9medhx.sas7bdat'))
    df = df.merge(dfv4[['ID', 'V4EALZH']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv5[['ID', 'V5SALZH']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv6[['ID', 'V6SALZH']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv8[['ID', 'V8EALZH']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv9[['ID', 'V9EALZH']], on='ID', how='left', validate='1:1')
    assert df.V4EALZH.sum()==0
    assert df.V5SALZH.sum()==0
    assert df.V6SALZH.sum()==0
    df = df.drop(columns=['V4EALZH', 'V5SALZH', 'V6SALZH'])
    df.loc[df.V8EALZH==1, 'V9EALZH'] = 1
    
    # medication
    dfv9 = pd.read_sas(os.path.join(data_dir, 'v9meds.sas7bdat'))
    df = df.merge(dfv9[['ID', 'V9ALZHM']], on='ID', how='left', validate='1:1')
    
    # cog test
    dfv1 = pd.read_sas(os.path.join(data_dir, 'v1cogfxn.sas7bdat'))
    dfv8 = pd.read_sas(os.path.join(data_dir, 'v8cogfxn.sas7bdat'))
    dfv9 = pd.read_sas(os.path.join(data_dir, 'v9cogfxn.sas7bdat'))
    dfcog = dfv1.merge(dfv8, on='ID', how='inner', validate='1:1').merge(dfv9, on='ID', how='inner', validate='1:1')
    df = df.merge(dfv1[['ID', 'V1SHT3MS']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv8[['ID', 'V8SHT3MS']], on='ID', how='left', validate='1:1')
    df = df.merge(dfv9[['ID', 'V9SHT3MS']], on='ID', how='left', validate='1:1')
    
    # diagnosis OR medication OR MMSE<-1.5std
    v8_cog_cutoff = (dfcog.V8SHT3MS-dfcog.V1SHT3MS).std()*1.5  # 2.95
    v9_cog_cutoff = (dfcog.V9SHT3MS-dfcog.V1SHT3MS).std()*1.5  # 5.64
    df['dement8'] = ((df.V8EALZH==1)|(df.V8SHT3MS-df.V1SHT3MS<-v8_cog_cutoff)).astype(int)
    #df['dement9'] = ((df.V9EALZH==1)|(df.V9ALZHM==1)|(df.V9SHT3MS-df.V1SHT3MS<-v9_cog_cutoff)).astype(int)
    
    #df2 = pd.read_stata(os.path.join(data_dir, 'SOF_HB_Cog.dta'))
    df2 = pd.read_csv('sof_data_bai_fromSasha.csv')
    df = df.merge(df2[['ID', 'APOE4', 'V9PRMCOG']], on='ID', how='inner', validate='1:1')
    df['dement9'] = df.V9PRMCOG.map({1:0, 2:1, 3:1})
    df.loc[df.dement8==1, 'dement9'] = 1
    
    for i in range(len(df)):
        if df.dement8.iloc[i]==0 and df.dement9.iloc[i]==0:
            df.loc[i, 'prevalent_dementia'] = 0
            df.loc[i, 'inc_dementia'] = 0
            if df.DEATH.iloc[i] == 0:
                df.loc[i, 'event'] = 'censor'
                df.loc[i, 'time2event'] = df.FOLALL.iloc[i]/365.25
            elif df.DEATH.iloc[i] == 1:
                df.loc[i, 'event'] = 'death'
                df.loc[i, 'time2event'] = df.FOLALL.iloc[i]/365.25
        #elif df.dement8.iloc[i]==1 and df.dement9.iloc[i]==0:
        #    raise NotImplementedError
        elif df.dement8.iloc[i]==0 and df.dement9.iloc[i]==1:
            df.loc[i, 'prevalent_dementia'] = 0
            df.loc[i, 'inc_dementia'] = 1
            df.loc[i, 'event'] = 'dementia'
            df.loc[i, 'time2event'] = df.SL_V9_DAYS.iloc[i]/365.25
        elif df.dement8.iloc[i]==1 and df.dement9.iloc[i]==1:
            df.loc[i, 'prevalent_dementia'] = 1
            df.loc[i, 'inc_dementia'] = 0
            df.loc[i, 'event'] = 'dementia'
            df.loc[i, 'time2event'] = df.SL_V8_DAYS.iloc[i]/365.25
    
    # education (and race), race not used, since AA has no V9 mortality
    df2 = pd.read_sas(os.path.join(data_dir, 'v1demogr.sas7bdat'))
    df2 = df2.rename(columns={'V1RACE':'Race'})
    df2['educcollege'] = (df2.V1EDUC>=12).astype(float)
    df2.loc[df2.V1EDUC.isna(), 'educcollege'] = np.nan
    df = df.merge(df2[['ID', 'educcollege', 'Race']], on='ID', how='left', validate='1:1')

    # BMI
    df2 = pd.read_sas(os.path.join(data_dir, 'v8anthro.sas7bdat'))
    df2 = df2.rename(columns={'V8BMI':'BMI'})
    df = df.merge(df2[['ID', 'BMI']], on='ID', how='left', validate='1:1')
    
    # sleep medication
    df2 = pd.read_sas(os.path.join(data_dir, 'v8meds.sas7bdat'))
    df2 = df2.rename(columns={'V8SLPMD':'sleepmed'})
    df = df.merge(df2[['ID', 'sleepmed']], on='ID', how='left', validate='1:1')
    
    # walking
    df2 = pd.read_sas(os.path.join(data_dir, 'v8lifestyle.sas7bdat'))
    df2 = df2.rename(columns={'V8EXER':'walking'})
    df = df.merge(df2[['ID', 'walking']], on='ID', how='left', validate='1:1')
    
    # diabetes + hypertension + heartattack + stroke
    df2 = pd.read_sas(os.path.join(data_dir, 'v8medhx.sas7bdat'))
    df2 = df2.rename(columns={'V8EDIAB':'diabetes', 'V8EHYPER':'hypertension', 'V8ESTRK':'stroke', 'V8EHEART':'heartattack'})
    df = df.merge(df2[['ID', 'diabetes', 'hypertension', 'stroke', 'heartattack']], on='ID', how='left', validate='1:1')
    
    # depression
    df2 = pd.read_sas(os.path.join(data_dir, 'v8qol.sas7bdat'))
    df2['depression'] = (df2.V8GDS15>=6).astype(float)
    df2.loc[df2.V8GDS15.isna(), 'depression'] = np.nan
    df = df.merge(df2[['ID', 'depression']], on='ID', how='left', validate='1:1')
    
    # mmse
    df2 = pd.read_sas(os.path.join(data_dir, 'v8cogfxn.sas7bdat'))
    df2 = df2.rename(columns={'V8MMSE':'mmse'})
    df = df.merge(df2[['ID', 'mmse']], on='ID', how='left', validate='1:1')
    
    # AHI
    df2 = pd.read_sas(os.path.join(data_dir, 'v8sleep.sas7bdat'))
    df2 = df2.rename(columns={'V8RDI3P':'AHI'})
    df = df.merge(df2[['ID', 'AHI']], on='ID', how='left', validate='1:1')
    
    cols = ['ID', 'prevalent_dementia', 'inc_dementia', 'event', 'time2event',
    'BAI', 'BA',
    'COUPL_OVERLAP_C', 'DENS_C',
       'alpha_bandpower_kurtosis_C_N2', 'alpha_bandpower_mean_C_N1',
       'delta_alpha_mean_C_N3', 'delta_bandpower_kurtosis_C_N2',
       'delta_bandpower_mean_C_N3', 'delta_theta_mean_C_N3', 'kurtosis_N2_C',
       'kurtosis_N3_C', 'sigma_bandpower_kurtosis_C_N2',
       'theta_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N3',
    'age', 'educcollege', 'BMI', 'Race',
    'sleepmed', 'walking', 'diabetes', 'hypertension', 'heartattack', 'stroke', 'depression',
    'mmse', 'AHI', 'APOE4']
    df = df[cols]
    
    print(df)
    breakpoint()
    df.to_csv('dataset_SOF.csv', index=False)
    


if __name__=='__main__':
    main()

