import os, subprocess
import numpy as np
import pandas as pd
from tqdm import tqdm


def main():
    df = pd.read_excel('SOF_mastersheet-after-step2.xlsx')
    print(df.shape)
    df = df[df.Age.notna()].reset_index(drop=True)
    print(df.shape)
    
    linenoise = 60
    model_path = 'BAI_model_luna'
    cwd = os.getcwd()
    for i in tqdm(range(len(df))):
        pid = df.SID.iloc[i]
        age = df.Age.iloc[i]
        #try:
        edf_path = df.CleanEDFPath.iloc[i]
        eeg_ch_names = 'C3A2,C4A1'
        ecg_ch_name = 'ECG'
        covar_path = os.path.join(cwd, f'covar_{pid}.txt')
        with open(covar_path, 'w') as f:
            f.write('ID\tage\tcen\tecg\tlinenoise\n')
            f.write(f'{pid}\t{age}\t{eeg_ch_names}\t{ecg_ch_name}\t{linenoise}')
        
        lst_path = os.path.join(cwd, f's_{pid}.lst')
        with open(lst_path, 'w') as f:
            f.write(f'{pid}\t{edf_path}')#\t{xml_path}
        
        output_db_path = os.path.join(cwd, f'output_{pid}.db')
        with open(os.path.join(model_path, 'm1-adult-age-luna.txt'), 'r') as f:
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
            continue
        ba = float(df_res.Y1.iloc[0])
        #ba_uncorrected = float(df_res.Y.iloc[0])
        
        df_features = pd.read_csv(output_csv_path2).T
        df_features.columns = df_features.iloc[0]
        df_features = df_features.iloc[1:].reset_index(drop=True)
        
        df.loc[i, 'BA'] = ba
        df.loc[i, 'BAI'] = ba-age
        for col in df_features.columns:
            df.loc[i,col] = df_features.loc[0,col]

        cols = ['SID', 'Age', 'BA', 'BAI']+list(df_features.columns)
        if i%10==9:
            df[cols].to_csv('BAI_SOF-tmp.csv', index=False)
        
        os.remove(covar_path)
        os.remove(lst_path)
        os.remove(output_db_path)
        os.remove(output_csv_path)
        os.remove(rcode_path)
        os.remove(output_csv_path2)

    print(df[cols])
    df[cols].to_excel('BAI_SOF.xlsx', index=False)
    

if __name__ == '__main__':
    main()

