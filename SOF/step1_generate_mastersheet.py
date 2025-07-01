import os, datetime
import numpy as np
import pandas as pd
import mne
#import pyedflib
from tqdm import tqdm
import warnings
    

def main():
    root_dir = '/data/haoqisun/dataset_MrOS_SOF_PSG/SOF'
    edf_paths = [x for x in os.listdir(root_dir) if x.lower().endswith('.edf')]
    ann_paths = [x for x in os.listdir(root_dir) if x.lower().endswith('_nsrr.xml')]
    sids = sorted(set([x[:-4] for x in edf_paths]) & set([x[:-9] for x in ann_paths]), key=int)
    
    df = pd.DataFrame(data={'SID':sids})
    df['EDFPath'] = ''
    df['AnnotPath'] = ''
    df['Channels'] = ''
    for i in tqdm(range(len(df))):
        sid = df.SID.iloc[i]
        edf_path = [x for x in edf_paths if sid.lower()+'.edf'==x.lower()]; assert len(edf_path)==1
        df.loc[i, 'EDFPath'] = os.path.join(root_dir, edf_path[0])
        ann_path = [x for x in ann_paths if sid.lower()+'_nsrr.xml'==x.lower()]; assert len(ann_path)==1
        df.loc[i, 'AnnotPath'] = os.path.join(root_dir, ann_path[0])
            
        warnings.filterwarnings("error")
        try:
            edf = mne.io.read_raw_edf(df.loc[i, 'EDFPath'], preload=False, verbose=False)
        except RuntimeWarning as e:
            msg = str(e)
            msg = msg[msg.index('(')+1:msg.rindex(')')]
            msg = msg.replace('-00', '-01')
            start_time = datetime.datetime.strptime(msg, '%Y-%m-%d %H:%M:%S')
        
        warnings.filterwarnings("ignore")
        edf = mne.io.read_raw_edf(df.loc[i, 'EDFPath'], preload=False, verbose=False)
        df.loc[i, 'Channels'] = str(edf.ch_names)
        df.loc[i, 'StartTime'] = start_time
        #df.loc[i, 'Fs'] = f.getSampleFrequncies()
    
    # add age
    dfa = pd.concat([
        pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/SOF-datasets/v8demogr.sas7bdat'),
        pd.read_sas('/data/haoqisun/dataset_MrOS_SOF_PSG/SOF-datasets/v8aademogr.sas7bdat'),
        ], axis=0, ignore_index=True)
    dfa = dfa.rename(columns={'V8AGE':'Age', 'ID':'SID'})
    dfa = dfa[dfa.SID.notna()].reset_index(drop=True)
    dfa['SID'] = dfa.SID.astype(int).astype(str)
    df = df.merge(dfa[['SID', 'Age']], on='SID', how='inner', validate='1:1')
    
    df = df[['SID', 'Age', 'StartTime', 'Channels', 'EDFPath', 'AnnotPath']]
    print(df)
    df.to_excel('SOF_mastersheet.xlsx', index=False)
    

if __name__=='__main__':
    main()

