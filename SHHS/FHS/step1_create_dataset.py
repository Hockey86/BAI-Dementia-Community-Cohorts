import os
import numpy as np
import pandas as pd


def main():
    df = pd.read_excel('../BAI_SHHS.xlsx')
    df = df[df.visitnumber=='shhs1'].reset_index(drop=True)
    
    df_shhs = pd.read_csv('../shhs1-dataset-0.20.0.csv')
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
    cols = ['BMI', 'diabetes', 'hypertension', 'AHI', 'heartattack', 'stroke', 'depression', 'sleepmed', 'smokstat_s1']
    df = df.merge(df_shhs[['nsrrid', 'pptidr']+cols], on='nsrrid', how='inner', validate='1:1')

    df_link = pd.read_sas('../parent_shhs_public_2016.sas7bdat')
    df_link = df_link[df_link.parent==b'FOFF'].reset_index(drop=True)#|(df_link.parent==b'OM1')
    df_link = df_link.rename(columns={'pid':'PID'})
    df_link['PID'] = df_link.PID.astype(int)
    df_link['parent'] = df_link.parent.astype(str)
    df_link['pptidr'] = df_link.pptidr.astype(str).str.strip()
    df = df.merge(df_link[['PID', 'pptidr', 'days_studyv1', 'parent']], on='pptidr', how='inner', validate='1:1')
    
    # get demographic info from FHS
    data_dir = '/data/haoqisun/dataset_Framingham_Offspring_2023b'
    df_cov = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/vr_dates_2019_a_1175d.csv'))
    df_cov['sexM'] = (df_cov.sex==1).astype(int)
    df_cov = df_cov[['PID', 'sexM']]

    df_race = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/race_1.csv'))
    df_race['race_Asian'] = (df_race.race=='A').astype(int)
    df_race['race_Black'] = (df_race.race=='B').astype(int)
    df_race['race_White'] = (df_race.race=='W').astype(int)
    df_race['race_Other'] = ((df_race.race=='O')|df_race.race.isna()).astype(int)
    #df_race['race_Unknown'] = df_race.race.isna().astype(int)
    df_race['ethnicity_Hispanic'] = (df_race.ethnicity=='Hisp').astype(float)
    df_race.loc[df_race.ethnicity.isna(), 'ethnicityHispanic'] = np.nan
    cols = ['race_Asian', 'race_Black', 'race_White', 'race_Other', 'ethnicity_Hispanic']#, 'race_Unknown'
    df_cov = df_cov.merge(df_race[['PID']+cols], on='PID', how='left', validate='1:1')

    df_edu = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/vr_educ_2018_a_1307d.csv'))
    df_edu = df_edu[pd.notna(df_edu.education)].reset_index(drop=True)
    educ_mapping = { 'a:noHSdeg':0, 'b:HSdeg':0, 'c:somecoll':1, 'd:collgrad':1, }
    df_edu['educcollege'] = df_edu.education.apply(lambda x:educ_mapping.get(x,x)).astype(int)
    cols = ['educcollege']
    df_cov = df_cov.merge(df_edu[['PID']+cols], on='PID', how='left', validate='1:1')

    df_pase = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/t_bmd_ex07_1_0104d_v1.csv'))
    df_pase = df_pase.rename(columns={'p6_7pase':'pascore'})
    cols = ['pascore']
    df_cov = df_cov.merge(df_pase[['PID']+cols], on='PID', how='left', validate='1:1')
    
    df_pase = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/e_exam_ex09_1b_0844d.csv'))
    df_pase = df_pase.rename(columns={'j638':'walkforexercise'})
    cols = ['walkforexercise']
    df_cov = df_cov.merge(df_pase[['PID']+cols], on='PID', how='left', validate='1:1')

    df_cog = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/vr_mmse_ex09_1b_0943d.csv'))
    df_cog = df_cog[df_cog.exam==6].reset_index(drop=True)
    df_cog = df_cog.rename(columns={'cogscr':'mmse'})
    cols = ['mmse']
    df_cov = df_cov.merge(df_cog[['PID']+cols], on='PID', how='left', validate='1:1')

    df = df.merge(df_cov, on='PID', how='inner', validate='1:1')

    df_dem = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/vr_demsurv_2018_a_1281d.csv'))
    cols_outcome = ['DEM_STATUS', 'DEM_SURVDATE']
    df = df.merge(df_dem[['PID']+cols_outcome], on='PID', how='inner', validate='1:1')  # outcome can be empty, so use how=inner
    df_dth = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/vr_survdth_2019_a_1337d.csv'))
    cols = ['DATEDTH']
    df = df.merge(df_dth[['PID']+cols], on='PID', how='inner', validate='1:1')
    for i in range(len(df)):
        if df.DEM_STATUS.iloc[i]==1:
            df.loc[i, 'event'] = 'dementia'
            df.loc[i, 'time2event'] = df.DEM_SURVDATE.iloc[i]/365.25
        elif pd.notna(df.DATEDTH.iloc[i]):
            df.loc[i, 'event'] = 'death'
            df.loc[i, 'time2event'] = df.DATEDTH.iloc[i]/365.25
        else:
            df.loc[i, 'event'] = 'censor'
            df.loc[i, 'time2event'] = df.DEM_SURVDATE.iloc[i]/365.25
    df['time2event'] = df.time2event - df.days_studyv1/365.25
    df['prevalent_dementia'] = ((df.event=='dementia')&(df.time2event<=0)).astype(int)
    df['inc_dementia'] = ((df.event=='dementia')&(df.time2event>0)).astype(int)
    
    # APOE
    
    df2 = pd.read_csv(os.path.join(data_dir, 'Datasets/CSV/eisoform1_4d.csv'))
    mapping = {0:0, 1:0, 2:1, 3:0, 4:1, 5:2, 6:0, 7:0, 8:np.nan, 9:np.nan}
    df2['APOE4Count'] = df2.E_TYPE.apply(lambda x:mapping[x])
    df = df.merge(df2[['PID', 'APOE4Count']], on='PID', how='left', validate='1:1')
    
    cols = ['visitnumber', 'nsrrid', 'PID','pptidr', 'parent', 'days_studyv1',
        'prevalent_dementia', 'inc_dementia', 'time2event', 'event',
        'age', 'sexM', 'race_Asian', 'race_Black',
        'race_White', 'race_Other', 'ethnicity_Hispanic',# 'race_Unknown',
        'educcollege', 'BMI', 'diabetes', 'hypertension', 'AHI', 'heartattack',
        'stroke', 'depression', 'sleepmed', 'smokstat_s1',  'pascore', 'walkforexercise', 'mmse',
        'BA', 'BAI', 'APOE4Count','COUPL_OVERLAP_C',
        'DENS_C', 'alpha_bandpower_kurtosis_C_N2', 'alpha_bandpower_mean_C_N1',
        'delta_alpha_mean_C_N3', 'delta_bandpower_kurtosis_C_N2',
        'delta_bandpower_mean_C_N3', 'delta_theta_mean_C_N3', 'kurtosis_N2_C',
        'kurtosis_N3_C', 'sigma_bandpower_kurtosis_C_N2',
        'theta_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N3']
    df = df[cols]
    print(df)
    breakpoint()
    df.to_csv('dataset_FHS.csv', index=False)


if __name__=='__main__':
    main()

