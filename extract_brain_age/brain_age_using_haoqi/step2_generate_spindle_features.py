import datetime
import os
import subprocess
import numpy as np
import scipy.io as sio
from scipy.signal import resample
import pandas as pd
import pyedflib
from tqdm import tqdm
from myfunctions import delete_nan_start_end, load_mros_sof_data


def convert_to_edf(signal_path, annot_path, output_path, dataset):
    """
    """
    # load dataset
    if dataset=='MrOS_SOF':
        EEG, sleep_stages, df_event, params = load_mros_sof_data(signal_path, annot_path)
    else:
        raise NotImplementedError(f'Unknown dataset={dataset}')
    EEG_channels = params['ch_names']
    Fs = params['Fs']
    combined_EEG_channels = ['C'] #TODO
    sleep_stages = sleep_stages[np.arange(0, len(sleep_stages), int(round(Fs*30)))]

    # resample
    newFs = 200.
    EEG = resample(EEG, int(round(EEG.shape[1]/Fs*newFs)), axis=1)
    Fs = newFs
            
    channel_info = [
            {'label': EEG_channels[i],
             'dimension': 'uV',
             'sample_rate': Fs,
             'physical_max': 32767,
             'physical_min': -32768,
             'digital_max': 32767,
             'digital_min': -32768,
             'transducer': 'E',
             'prefilter': ''}
        for i in range(len(EEG_channels))]

    with pyedflib.EdfWriter(output_path, len(EEG), file_type=pyedflib.FILETYPE_EDFPLUS) as ff:
        ff.setSignalHeaders(channel_info)
        ff.writeSamples(EEG)
            
    return sleep_stages, EEG_channels, combined_EEG_channels


def convert_to_xml(sleep_stages, output_path):
    """
    """       
    sleep_stages[np.isnan(sleep_stages)] = -1
    sleep_stage_mapping = {-1:0, 0:0, 5:0, 4:5, 3:1, 2:2, 1:3}
    with open(output_path, 'w') as ff:
        ff.write('<CMPStudyConfig>\n')
        ff.write('<EpochLength>30</EpochLength>\n')
        ff.write('<SleepStages>\n')
        for ss in sleep_stages:
            ff.write('<SleepStage>%d</SleepStage>\n'%sleep_stage_mapping[ss])
        ff.write('</SleepStages>\n')
        ff.write('</CMPStudyConfig>')


def run_luna(output_dir, sids, edf_paths, xml_paths, channels, stage, cfreq=13.5, cycles=12):
    now = datetime.datetime.now().strftime('%H%M%S_%m%d%Y')
    # create the list file
    list_path = os.path.join(output_dir, f'luna_{now}.lst')
    df = pd.DataFrame(data={'sid': sids, 'edf':edf_paths, 'xml':xml_paths})
    df = df[['sid', 'edf', 'xml']]
    df.to_csv(list_path, sep='\t', index=False, header=False)

    #channels = [x.replace('-','_') for x in channels]
    
    # generate R code to convert luna output .db to .mat
    db_path = os.path.join(output_dir, f'luna_output_{now}.db')
    csv_path = os.path.join(output_dir, f'luna_output_{now}.csv')
    db_path = db_path.replace('\\', '\\\\')  # windows only
    csv_path = csv_path.replace('\\', '\\\\')  # windows only
    r_code = f"""library(luna)
k<-ldb("{db_path}")
d1<-lx(k,"SPINDLES", "CH_F")
d2<-lx(k,"SPINDLES", "CH")
d3 <- merge(d1,d2,by=c("ID","CH"))
write.csv(d3, "{csv_path}", row.names=F)
"""
    r_code_path = os.path.join(output_dir, f'convert_luna_output_db2mat_{now}.R')
    with open(r_code_path, 'w') as ff:
        ff.write(r_code)
    # maximum spindle duration 10s to account for pediatric subjects
    #with open(os.devnull, 'w') as FNULL:
    subprocess.check_call(['luna', list_path, '-o', db_path, '-s',
        'MASK ifnot=%s & RE & SPINDLES sig=%s max=10 fc=%g cycles=%d so mag=1.5'%(stage, ','.join(channels), cfreq, cycles)],)
    #    stdout=FNULL, stderr=subprocess.STDOUT)
    subprocess.check_call(['Rscript', r_code_path])#, stdout=FNULL, stderr=subprocess.STDOUT)
    
    df = pd.read_csv(csv_path)
    
    # delete intermediate files
    if os.path.exists(csv_path):
        os.remove(csv_path)
    if os.path.exists(r_code_path):
        os.remove(r_code_path)
    if os.path.exists(db_path):
        os.remove(db_path)
    if os.path.exists(list_path):
        os.remove(list_path)
        
    return df
    

if __name__=='__main__':
    dataset = 'MrOS_SOF'
    cfreq = 13.5  # [Hz]
    cycles = 12
    stage = 'NREM2'

    # get list of files
    df = pd.read_csv(f'../mastersheet_{dataset}.csv')
    assert not df.SID.str.contains(' ').any()
    
    cwd = os.getcwd()
    output_dir = os.path.join(cwd, f'features_{dataset}')
    os.makedirs(output_dir, exist_ok=True)
    xml_edf_dir = os.path.join(cwd, 'edf_xml_for_luna')
    os.makedirs(xml_edf_dir, exist_ok=True)

    # generate edf and xmls files and paths required by Luna
    for i in tqdm(range(len(df))):
        sid = df.SID.iloc[i]
        signal_path = df.EDFPath.iloc[i]
        annot_path = df.AnnotPath.iloc[i]
        edf_path = os.path.join(xml_edf_dir, sid+'.edf')
        xml_path = os.path.join(xml_edf_dir, sid+'.xml')

        try:
            sleep_stages, EEG_channels, combined_EEG_channels = convert_to_edf(signal_path, annot_path, edf_path, dataset)
            convert_to_xml(sleep_stages, xml_path)
        except Exception as ee:
            print(ee)
            continue

        df_res = run_luna(xml_edf_dir, [sid], [edf_path], [xml_path],
            EEG_channels, stage, cfreq=cfreq, cycles=cycles)

        # save space
        os.remove(edf_path)
        os.remove(xml_path)
        
        # make one row per subject
        for chn in combined_EEG_channels:
            df_res.loc[df_res.CH.str.startswith(chn),'CH2'] = chn
        df_res = df_res.iloc[:,3:].groupby('CH2').agg('mean').reset_index()

        df_res2 = []
        for chn in combined_EEG_channels:
            df_res_ = df_res[df_res.CH2==chn].drop(columns='CH2')
            df_res_ = df_res_.rename(columns={x:f'{x}_{chn}' for x in df_res_.columns}).reset_index(drop=True)
            df_res2.append(df_res_)
        df_res = pd.concat(df_res2, axis=1)
        df_res.insert(0, 'SID', sid)
        df_res.to_csv(os.path.join(output_dir, f'spindle_features_{sid}.csv'), index=False)

