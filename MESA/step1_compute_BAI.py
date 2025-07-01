import datetime, os, subprocess
from itertools import groupby
import numpy as np
import pandas as pd
from scipy.stats import mode
from tqdm import tqdm
import mne
import pyedflib
import xmltodict
from scipy.signal import peak_prominences
from neurokit2 import ecg_peaks


def load_mesa_data(sid, signal_path=None, annot_path=None, download=True, dropbox_file_list=None):
    if (signal_path is None or annot_path is None) and download:
        signal_path = f'mesa-sleep-{sid}.edf'
        annot_path  = f'mesa-sleep-{sid}-nsrr.xml'
        if not os.path.exists(signal_path):
            subprocess.run(['rclone', 'copy', 'dropbox_bidmc:/Datasets/zz_SLEEP/mesa/polysomnography/edfs/'+signal_path, '.'])
        if not os.path.exists(annot_path):
            subprocess.run(['rclone', 'copy', 'dropbox_bidmc:/Datasets/zz_SLEEP/mesa/polysomnography/annotations-events-nsrr/'+annot_path, '.'])
    
    edf = mne.io.read_raw_edf(signal_path, verbose=False, preload=False)
    Fs = edf.info['sfreq']
    signals = edf.get_data(picks=['EEG3', 'EKG'])*1e6
    assert len(signals)==2
    signals = np.array([signals[0], signals[0], signals[1]])
    
    with open(annot_path, 'r') as f:
        annot = xmltodict.parse(f.read())
    annot = pd.DataFrame(data=annot['PSGAnnotation']['ScoredEvents']['ScoredEvent'])
    annot['Start'] = annot.Start.astype(float)
    annot['Duration'] = annot.Duration.astype(float)
    
    txt2num = {'Wake|0':5, 'Stage 1 sleep|1':3, 'Stage 2 sleep|2':2, 'Stage 3 sleep|3':1, 'REM sleep|5':4}
    sleep_stages = np.zeros(signals.shape[1])+np.nan
    for i, r in annot[annot.EventType=='Stages|Stages'].iterrows():
        start = int(round(r.Start*Fs))
        end = start + int(round(r.Duration*Fs))
        start = max(0, start)
        end = min(signals.shape[1], end)
        if start<end and r.EventConcept in txt2num:
            sleep_stages[start:end] = txt2num[r.EventConcept]
    
    params = {'Fs':Fs,
        'start_time':edf.info['meas_date'].replace(tzinfo=None),
        'ch_names':['C3M2', 'C4M1', 'ECG'] }
    
    if download:
        if os.path.exists(signal_path):
            os.remove(signal_path)
        if os.path.exists(annot_path):
            os.remove(annot_path)
    return signals, sleep_stages, annot, params


def flip_ecg(ecg_, Fs, line_freq=60):
    # returns True/False to indicate if ecg should be flipped
    # returns None if (1) has <10 minutes non-nan signal or (2) overall heart rate is abnormal
    
    # find the longest non-nan and non-flat part
    nan_mask = np.isnan(ecg_)|np.isinf(ecg_)|(np.abs(ecg_)>5000)
    cc = 0; maxll = 0; start = -1
    for k,l in groupby(nan_mask):
        ll = len(list(l))
        if (not k) and ecg_[cc:cc+ll].std()>1 and ll>maxll:
            maxll = ll
            start = cc
        cc += ll
    if start<0:
        return None
    ecg_seg = ecg_[start:start+maxll]
    if len(ecg_seg)<Fs*30:
        return None
    # find middle two hours, so that the signal is short
    if len(ecg_seg)>Fs*7200:
        start = int(len(ecg_seg)//2-Fs*3600)
        end = int(len(ecg_seg)//2+Fs*3600)
        ecg_seg = ecg_seg[start:end]
    if line_freq<Fs/2:
        ecg_seg = mne.filter.notch_filter(ecg_seg, Fs, line_freq, verbose=False)
    ecg_seg = mne.filter.filter_data(ecg_seg, Fs, 5, 70 if 70<Fs/2 else None, verbose=False)
    try:
        #rpeaks1 = detect_heartbeats(ecg_seg, Fs)
        #rpeaks2 = detect_heartbeats(-ecg_seg, Fs)
        rpeaks1 = ecg_peaks(ecg_seg, sampling_rate=Fs)[1]['ECG_R_Peaks']
        rpeaks2 = ecg_peaks(-ecg_seg, sampling_rate=Fs)[1]['ECG_R_Peaks']
    except Exception as ee:
        print(str(ee))
        return None
    
    peakness1 = np.mean(peak_prominences(ecg_seg, rpeaks1)[0])
    peakness2 = np.mean(peak_prominences(-ecg_seg, rpeaks2)[0])
    if peakness1<peakness2:
        hr = len(rpeaks2)/(len(ecg_seg)/Fs/60)
    else:
        hr = len(rpeaks1)/(len(ecg_seg)/Fs/60)
    if 30<hr<100:
        return peakness1<peakness2
    else:
        return None


def remove_start_end_wake(signals, sleep_stages, Fs, start_time):
    sleep_stages2 = np.array(sleep_stages)
    sleep_stages2[np.isnan(sleep_stages2)] = 5
    start = 0
    for k,l in groupby(sleep_stages2):
        if k==5:
            start = len(list(l))
        break
    end = 0
    for k,l in groupby(sleep_stages2[::-1]):
        if k==5:
            end = len(list(l))
        break
    end = len(sleep_stages)-end
    
    return signals[...,start:end], sleep_stages[start:end], start_time+datetime.timedelta(seconds=start/Fs)


def get_brain_age_luna(model_path, sid, signals, sleep_stages, ch_names, start_time, Fs, age, exclude_start_stop_ids=None):
    cwd = os.getcwd()
    assert ' ' not in cwd
    
    # save edf
    edf_path = os.path.join(cwd, f'luna_{sid}.edf')
    headers = [
            {'label': ch_names[i],
             'dimension': 'uV',
             'sample_rate': Fs,
             'physical_max': 32767,
             'physical_min': -32768,
             'digital_max': 32767,
             'digital_min': -32768,
             'transducer': 'E',
             'prefilter': ''}
        for i in range(len(ch_names))]
    with pyedflib.EdfWriter(edf_path, len(signals), file_type=pyedflib.FILETYPE_EDF) as ff:
        ff.setSignalHeaders(headers)
        ff.setStartdatetime(start_time)
        ff.writeSamples(signals)
    
    # get per-epoch sleep stages
    sleep_stages2 = np.array(sleep_stages)
    if exclude_start_stop_ids is not None:
        sleep_stages2[exclude_start_stop_ids[0]:exclude_start_stop_ids[1]] = -1
    epoch_size = int(round(Fs*30))
    start_ids = np.arange(0, len(sleep_stages2)-epoch_size+1, epoch_size)
    sleep_stages2 = np.array([sleep_stages2[x:x+epoch_size] for x in start_ids])
    sleep_stages2[(sleep_stages2==-1).any(axis=1)] = -1
    sleep_stages2 = np.array([mode(x, keepdims=False).mode for x in sleep_stages2])
    
    # save xml
    xml_path = os.path.join(cwd, f'luna_{sid}.xml')
    sleep_stage_mapping = {5: 0, 4: 5, 3: 1, 2: 2, 1: 3}
    with open(xml_path, 'w') as ff:
        ff.write('<CMPStudyConfig>\n')
        ff.write('<EpochLength>30</EpochLength>\n')
        ff.write('<SleepStages>\n')
        for ss in sleep_stages2:
            ff.write('<SleepStage>%d</SleepStage>\n'%sleep_stage_mapping.get(ss,9))
        ff.write('</SleepStages>\n')
        ff.write('</CMPStudyConfig>')
    
    lst_path = os.path.join(cwd, f's_{sid}.lst')
    with open(lst_path, 'w') as f:
        f.write(f'{sid}\t{edf_path}\t{xml_path}')
    
    covar_path = os.path.join(cwd, 'covar.txt')
    ch_names = [x for x in ch_names if x.startswith('C')]
    ch_names2 = ','.join([x.replace('-','_') for x in ch_names])
    with open(covar_path, 'w') as f:
        f.write('ID\tage\tcen\tecg\tlinenoise\n')
        f.write(f'{sid}\t{age}\t{ch_names2}\tECG\t60')
    
    output_db_path = os.path.join(cwd, f'output_{sid}.db')
    with open(os.path.join(model_path, 'm1-adult-age-luna.txt'), 'r') as f:
        subprocess.run(['luna', lst_path, 'vars='+covar_path, 'th=3', 'mpath='+model_path, '-o', output_db_path], stdin=f)
    
    output_csv_path = os.path.join(cwd, f'output_{sid}.csv')
    output_csv_path2 = os.path.join(cwd, f'output2_{sid}.csv')
    rcode = f'''library(luna)
k <- ldb("{output_db_path}")
res <- lx(k, "PREDICT", "BL")
write.csv(res, "{output_csv_path}", row.names = F)
res <- lx(k,"PREDICT","FTR")[,c("FTR","X")]
write.csv(res, "{output_csv_path2}", row.names = F)
'''
    rcode_path = os.path.join(cwd, f'code_{sid}.R')
    with open(rcode_path,'w') as f:
        f.write(rcode)
    subprocess.run(['Rscript', rcode_path])
    
    df_res = pd.read_csv(output_csv_path)
    ba = float(df_res.Y1.iloc[0])
    ba_uncorrected = float(df_res.Y.iloc[0])
    
    df_features = pd.read_csv(output_csv_path2).T
    df_features.columns = df_features.iloc[0]
    df_features = df_features.iloc[1:].reset_index(drop=True)
    
    os.remove(rcode_path)
    os.remove(output_csv_path2)
    os.remove(output_csv_path)
    os.remove(output_db_path)
    os.remove(covar_path)
    os.remove(lst_path)
    os.remove(edf_path)
    os.remove(xml_path)
    
    return ba, ba_uncorrected, df_features


def main():
    """
    df = pd.read_csv('mesa-sleep-dataset-0.3.0.csv')
    df = df.rename(columns={'gender1':'sexM', 'sleepage5c':'age'})
    df2 = pd.read_csv('MESAdatasets/mesa_nsrr_bridge_ids.csv')
    df = df.merge(df2, on='mesaid', how='inner', validate='1:1')
    
    res1 = subprocess.run(['rclone', 'ls', 'dropbox_bidmc:/Datasets/zz_SLEEP/mesa/polysomnography/edfs'], capture_output=True, text=True)
    sids1 = [int(x.split('-')[-1][:-4]) for x in res1.stdout.strip().split('\n')]
    res2 = subprocess.run(['rclone', 'ls', 'dropbox_bidmc:/Datasets/zz_SLEEP/mesa/polysomnography/annotations-events-nsrr'], capture_output=True, text=True)
    sids2 = [int(x.split('-')[-2]) for x in res2.stdout.strip().split('\n')]
    exist_sids = sorted(set(sids1)&set(sids2))
    
    df = df[np.in1d(df.mesaid, exist_sids)].reset_index(drop=True)
    df.to_csv('mesa-sleep-dataset-haoqi.csv', index=False)
    """
    df = pd.read_csv('mesa-sleep-dataset-haoqi.csv')
    
    #model_path = '/data/haoqisun/BAI_dementia_community/extract_brain_age/BAI_model_luna/'
    model_path = '/data/haoqisun/BAI_dementia_community/MESA/BAI_model_luna/'
    cols = ['mesaid', 'sexM', 'age', 'BA', 'BAI', 'COUPL_OVERLAP_C', 'DENS_C', 'alpha_bandpower_kurtosis_C_N2',
       'alpha_bandpower_mean_C_N1', 'delta_alpha_mean_C_N3',
       'delta_bandpower_kurtosis_C_N2', 'delta_bandpower_mean_C_N3',
       'delta_theta_mean_C_N3', 'kurtosis_N2_C', 'kurtosis_N3_C',
       'sigma_bandpower_kurtosis_C_N2', 'theta_bandpower_kurtosis_C_N2',
       'theta_bandpower_kurtosis_C_N3']
    
    for i in tqdm(range(len(df))):
        sid = df.mesaid.iloc[i]
        age = df.age.iloc[i]
        if pd.isna(age):
            continue
        #if sid not in exist_sids:
        #    continue
        
        try:
            signals, sleep_stages, annot, params = load_mesa_data(f'{sid:04d}', download=True)
            Fs = params['Fs']
            start_time = params['start_time']
            ch_names = params['ch_names']
            
            # preprocess
            signals, sleep_stages, start_time = remove_start_end_wake(signals, sleep_stages, Fs, start_time)
            signals = mne.filter.notch_filter(signals, Fs, 60, verbose=False)
            eeg = signals[:2]
            #eeg = mne.filter.filter_data(eeg, Fs, 0.3, 35, verbose=False)
            ecg = signals[2]
            ecg = mne.filter.filter_data(ecg, Fs, 0.3, 70, verbose=False)
            
            flip = flip_ecg(ecg, Fs, line_freq=60)
            assert flip is not None
            if flip:
                ecg = -ecg
            
            # save to luna-compatible format
            ba, ba_uncorrected, df_feat = get_brain_age_luna(model_path, sid,
                np.vstack([eeg,ecg]), sleep_stages, ch_names, start_time, Fs, age)
            df.loc[i, 'BAI'] = ba-age
            df.loc[i, 'BA'] = ba
            #df.loc[i, 'BA_uncorrected'] = ba_uncorrected
            for col in df_feat.columns:
                df.loc[i, col] = df_feat[col].iloc[0]
            if i%10==1:
                df[cols].to_csv('BAI_MESA.csv', index=False)
        except Exception as ee:
            print(f'{sid}: {ee}')
            #continue

    df = df[cols]
    print(df)
    df.to_csv('BAI_MESA.csv', index=False)



if __name__=='__main__':
    main()

