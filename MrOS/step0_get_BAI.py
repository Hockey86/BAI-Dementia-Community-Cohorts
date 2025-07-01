"""
originally the BAIs come from
/data/haoqisun/BAI_dementia_community/extract_brain_age/brain_age_using_luna/BAI_results.csv

But we need the features as well. So re-compute them here.
"""
import os, subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm


def compute_BAI(model_path, pid, edf_path, xml_path, age, ch_names):
    eeg_ch_names = ['C3-A2', 'C4-A1']
    ecg_ch_name = 'ECG L-ECG R'
    eeg_ch_names2 = ','.join(eeg_ch_names).replace('-','_').replace(' ','_')
    ecg_ch_name2 = ecg_ch_name.replace('-','_').replace(' ','_')
    cwd = os.getcwd()
    covar_path = os.path.join(cwd, f'covar_{pid}.txt')
    with open(covar_path, 'w') as f:
        f.write('ID\tage\tcen\tecg\tlinenoise\n')
        f.write(f'{pid}\t{age}\t{eeg_ch_names2}\t{ecg_ch_name2}\t60')
    
    lst_path = os.path.join(cwd, f's_{pid}.lst')
    with open(lst_path, 'w') as f:
        f.write(f'{pid}\t{edf_path}\t{xml_path}')
     
    model_path1 = os.path.join(model_path, 'm1-adult-age-luna1.txt')
    model_path2 = os.path.join(model_path, 'm1-adult-age-luna2.txt')
    ch_names1 = ['C3', 'C4', 'A1', 'A2', 'ECG L', 'ECG R']
    ch_names2 = ['C3-A2', 'C4-A1', 'ECG L-ECG R']
    if all([x in ch_names for x in ch_names1]):
        model_path_ = model_path1
    elif all([x in ch_names for x in ch_names2]):
        model_path_ = model_path2
    else:
        raise SystemExit

    output_db_path = os.path.join(cwd, f'output_{pid}.db')
    with open(model_path_ , 'r') as f:
        subprocess.run(['luna', lst_path, 'vars='+covar_path, 'th=3', 'mpath='+model_path, '-o', output_db_path], stdin=f)
    
    output_csv_path = os.path.join(cwd, f'output_{pid}.csv')
    output_csv_path2 = os.path.join(cwd, f'output2_{pid}.csv')
    rcode = f'''library(luna)
k <- ldb("{output_db_path}")
res <- lx(k, "PREDICT", "BL")
write.csv(res, "{output_csv_path}", row.names = F)
res <- lx(k,"PREDICT","FTR")[,c("FTR","X")]
write.csv(res, "{output_csv_path2}", row.names = F)
'''
    rcode_path = os.path.join(cwd, f'code_{pid}.R')
    with open(rcode_path,'w') as f:
        f.write(rcode)
    subprocess.run(['Rscript', rcode_path])
    
    df_res = pd.read_csv(output_csv_path)
    if 'Y1' not in df_res.columns:# or 'Y1' not in df_res.columns:
        ba = np.nan
    else:
        ba = float(df_res.Y1.iloc[0])
    #ba_uncorrected = float(df_res.Y.iloc[0])
    
    feat = pd.read_csv(output_csv_path2).T
    feat.columns = feat.iloc[0]
    feat = feat.iloc[1:].reset_index(drop=True)
    
    if os.path.exists(output_csv_path): os.remove(output_csv_path)
    if os.path.exists(output_csv_path2): os.remove(output_csv_path2)
    if os.path.exists(covar_path): os.remove(covar_path)
    if os.path.exists(lst_path): os.remove(lst_path)
    if os.path.exists(output_db_path): os.remove(output_db_path)
    if os.path.exists(rcode_path): os.remove(rcode_path)
    return ba, feat


def main():
    df = pd.read_csv('mastersheet_MrOS.csv')
    model_path = '/data/haoqisun/BAI_dementia_community/MrOS/BAI_model_luna'
    output_path = 'BAI_MrOS_with_features.csv'
    
    df2 = pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/MrOS-datasets/vsfeb24.sas7bdat')
    df2['SID'] = df2['ID'].astype(str).str.lower()
    df = df.merge(df2[['SID','VSAGE1']].rename(columns={'VSAGE1':'Age'}), on='SID', how='inner', validate='1:1')
    
    df['BA'] = np.nan
    df['BAI'] = np.nan
    for i in tqdm(range(len(df))):
        age = df.Age.iloc[i]
        channels = eval(df.Channel.iloc[i])
        #clean_path = clean_edf()
        
        ba, feat = compute_BAI(model_path,
            df.SID.iloc[i], df.EDFPath.iloc[i], df.AnnotPath.iloc[i],
            age, channels)
            
        df.loc[i,'BA'] = ba
        df.loc[i,'BAI'] = ba-age
        for col in feat.columns:
            df.loc[i,col] = feat.loc[0,col]
            
        if i%10==0:
            cols = ['SID', 'Age', 'BA', 'BAI']+list(feat.columns)
            df[cols].to_csv(output_path, index=False)
    
    df = df[cols]
    print(df)
    df.to_csv(output_path, index=False)


if __name__=='__main__':
    main()
    
