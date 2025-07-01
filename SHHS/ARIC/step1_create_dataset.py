import os
import numpy as np
import pandas as pd


def main():
    df = pd.read_excel('../BAI_SHHS.xlsx')
    df = df[df.visitnumber=='shhs1'].reset_index(drop=True)
    
    df_shhs = pd.read_csv('/data/haoqisun/dataset_SHHS/shhs1-dataset-0.20.0.csv')
    df_shhs = df_shhs[df_shhs.overall_shhs1>=3].reset_index(drop=True)
    df_shhs['pptidr'] = df_shhs.pptidr.astype(str).str.strip()
    m = {0:0, 1:1, 2:0}
    df_shhs['smokstat_s1'] = df_shhs.smokstat_s1.apply(lambda x:m.get(x,x))
    #df_shhs.loc[df_shhs.alcoh>50, 'alcoh'] = 50  # https://sleepdata.org/datasets/shhs/variables/alcoh/known-issues
    df_shhs = df_shhs.rename(columns={'bmi_s1':'BMI', 'parrptdiab':'diabetes', 'htnderv_s1':'hypertension', 'mi15':'heartattack', 'stroke15':'stroke', 'ahi_a0h4':'AHI'})
    df_shhs.loc[df_shhs.heartattack==8, 'heartattack'] = np.nan
    df_shhs.loc[df_shhs.stroke==8, 'stroke'] = np.nan
    df_shhs['depression'] = ((df_shhs.blue25<=3)|(df_shhs.down25<=3)).astype(float)
    df_shhs.loc[df_shhs.blue25.isna()&df_shhs.down25.isna(), 'depression'] = np.nan
    df_shhs['sleepmed'] = ((df_shhs.benzod1==1)|(df_shhs.tca1==1)|(df_shhs.ntca1==1)).astype(float)
    df_shhs.loc[df_shhs.benzod1.isna()&df_shhs.tca1.isna()&df_shhs.ntca1.isna(), 'sleepmed'] = np.nan
    df_shhs['race_White'] = (df_shhs.race==1).astype(int)
    df_shhs['race_Black'] = (df_shhs.race==2).astype(int)
    df_shhs['race_Other'] = (df_shhs.race==3).astype(int)
    df_shhs['sexM'] = (df_shhs.gender==1).astype(int)
    df_shhs['educcollege'] = (df_shhs.educat>=3).astype(int)
    cols = ['sexM', 'race_White', 'race_Black', 'race_Other', 'educcollege', 'BMI', 'diabetes', 'hypertension', 'AHI', 'heartattack', 'stroke', 'depression', 'sleepmed', 'smokstat_s1']
    df = df.merge(df_shhs[['nsrrid', 'pptidr']+cols], on='nsrrid', how='inner', validate='1:1')

    df_link = pd.read_sas('/data/haoqisun/dataset_SHHS/parent_shhs_public_2016.sas7bdat')
    df_link = df_link[df_link.parent==b'ARIC'].reset_index(drop=True)
    df_link['pptidr'] = df_link.pptidr.astype(str).str.strip()
    df = df.merge(df_link[['pptidr', 'days_studyv1']], on='pptidr', how='inner', validate='1:1')
    
    # get demographic info
    data_dir = '/data/haoqisun/dataset_ARIC_2024b/Main_Study'

    # get dementia survival info by merging with dementia survival variables
    # get death info
    df_dem = pd.read_csv(os.path.join(data_dir, 'V7/CSV/status71.csv'))
    breakpoint()
    df_dem = df_dem.rename(columns={'ID_C':'pptidr', 'DEMDXL3CENS_71':'DEM_STATUS', 'COXDATE_DEMDXL3_71_FOLLOWUPDAYS':'DEM_SURVDATE', 'DATEOFDEATH_FOLLOWUPDAYS':'DATEDTH', 'STATUSDATE71_FOLLOWUPDAYS':'FOLLOWUPDAYS'})
    df = df.merge(df_dem[['pptidr', 'DEM_STATUS', 'DEM_SURVDATE', 'DATEDTH', 'FOLLOWUPDAYS']], on='pptidr', how='inner', validate='1:1')
    print(f'after merging with df_dem, df.shape = {df.shape}')
    
    for i in range(len(df)):
        if df.DEM_STATUS.iloc[i]==1:
            df.loc[i, 'event'] = 'dementia'
            df.loc[i, 'time2event'] = df.DEM_SURVDATE.iloc[i]/365.25
        elif pd.notna(df.DATEDTH.iloc[i]):
            df.loc[i, 'event'] = 'death'
            df.loc[i, 'time2event'] = df.DATEDTH.iloc[i]/365.25
        else:
            df.loc[i, 'event'] = 'censor'
            df.loc[i, 'time2event'] = max(df.DEM_SURVDATE.iloc[i], df.FOLLOWUPDAYS.iloc[i])/365.25
    df['time2event'] = df.time2event - df.days_studyv1/365.25
    df['prevalent_dementia'] = ((df.event=='dementia')&(df.time2event<=0)).astype(int)
    df['inc_dementia'] = ((df.event=='dementia')&(df.time2event>0)).astype(int)
    
    """
    df2 = pd.read_csv(os.path.join(data_dir, 'Longitudinal/CSV/v1_v5_analytes.csv'))
    df2 = df2.rename(columns={'ID_C':'pptidr'})
    res = df.merge(df2,on='pptidr', how='inner',validate='1:1')[['age','AGE_V1', 'AGE_V2', 'AGE_V3', 'AGE_V4', 'AGE_V5']]
    print(res)
    """
    # sleep happens at V4
    
    """
    df_pa = pd.read_csv(os.path.join(data_dir, 'v5/CSV/pac.csv'))
    df_pa = df_pa.rename(columns={'ID_C':'pptidr', 'PAC1':'exercise'})
    df_pa.loc[df_pa.exercise=='Y', 'exercise'] = 1
    df_pa.loc[df_pa.exercise=='N', 'exercise'] = 0
    df = df.merge(df_pa[['pptidr', 'exercise']], on='pptidr', how='inner', validate='1:1')
    """
    df_pa = pd.read_csv(os.path.join(data_dir, 'v3/CSV/derive37.csv'))
    df_pa = df_pa.rename(columns={'ID_C':'pptidr', 'SPRT_I31':'exercise'})
    df = df.merge(df_pa[['pptidr', 'exercise']], on='pptidr', how='left', validate='1:1')

    df_cog = pd.read_csv(os.path.join(data_dir, 'v5/CSV/derive51.csv'))
    df_cog = df_cog.rename(columns={'ID_C':'pptidr', 'PRORATEDMMS51':'mmse'})
    df = df.merge(df_cog[['pptidr', 'mmse']], on='pptidr', how='left', validate='1:1')
    
    cols = ['visitnumber', 'nsrrid', 'pptidr', 'days_studyv1',
        'prevalent_dementia', 'inc_dementia', 'time2event', 'event',
        'age', 'sexM', 'race_Black', 'race_White', 'race_Other',
        'educcollege', 'BMI', 'diabetes', 'hypertension', 'AHI', 'heartattack',
        'stroke', 'depression', 'sleepmed', 'smokstat_s1',  'exercise', 'mmse',
        'BA', 'BAI','COUPL_OVERLAP_C',
        'DENS_C', 'alpha_bandpower_kurtosis_C_N2', 'alpha_bandpower_mean_C_N1',
        'delta_alpha_mean_C_N3', 'delta_bandpower_kurtosis_C_N2',
        'delta_bandpower_mean_C_N3', 'delta_theta_mean_C_N3', 'kurtosis_N2_C',
        'kurtosis_N3_C', 'sigma_bandpower_kurtosis_C_N2',
        'theta_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N3']
    df = df[cols]
    print(df)
    breakpoint()
    df.to_csv('dataset_ARIC.csv', index=False)


if __name__=='__main__':
    main()

