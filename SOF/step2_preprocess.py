from itertools import groupby
import os
import numpy as np
from scipy.signal import resample, peak_prominences
import pandas as pd
from tqdm import tqdm
import pyedflib
import mne
import xml2dict
from neurokit2 import ecg_peaks


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
    #if line_freq<Fs/2:
    #    ecg_seg = mne.filter.notch_filter(ecg_seg, Fs, line_freq, verbose=False)
    #ecg_seg = mne.filter.filter_data(ecg_seg, Fs, 5, 70 if 70<Fs/2 else None, verbose=False)
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
        
        
def main():
    df = pd.read_excel('SOF_mastersheet.xlsx')
    
    output_dir = 'clean-edf'
    os.makedirs(output_dir, exist_ok=True)
    
    for i in tqdm(range(len(df))):
        sid = df.SID.iloc[i]
        output_path = os.path.join(output_dir, f'{sid}.edf')
        df.loc[i, 'CleanEDFPath'] = output_path
        if os.path.exists(output_path):
            continue
        
        edf_path = df.EDFPath.iloc[i]
        ann_path = df.AnnotPath.iloc[i]
        start_time = df.StartTime.iloc[i]
        ch_names = eval(df.Channels.iloc[i])
        
        edf_eeg = mne.io.read_raw_edf(edf_path, verbose=False, preload=False, exclude=[x for x in ch_names if x not in ['C3','C4','A1','A2']])
        edf_ecg = mne.io.read_raw_edf(edf_path, verbose=False, preload=False, exclude=[x for x in ch_names if x not in ['ECG1','ECG2']])
    
        eeg = edf_eeg.get_data(picks=['C3','C4','A1','A2'])*1e6
        eeg = np.array([eeg[0]-eeg[3], eeg[1]-eeg[2]])
        fs_eeg = edf_eeg.info['sfreq']
        df.loc[i, 'FsEEG'] = fs_eeg
        
        ecg = edf_ecg.get_data(picks=['ECG1','ECG2'])*1e6
        ecg = ecg[0]-ecg[1]
        fs_ecg = edf_ecg.info['sfreq']
        df.loc[i, 'FsECG'] = fs_ecg
        
        Fs = min(fs_eeg, fs_ecg)
        df.loc[i, 'Fs'] = Fs
        if fs_eeg<fs_ecg:
            ecg = resample(ecg, int(round(ecg.shape[-1]/fs_ecg*fs_eeg)))
        elif fs_eeg>fs_ecg:
            eeg = resample(eeg, int(round(eeg.shape[-1]/fs_eeg*fs_ecg)))
        
        if ecg.shape[-1]!=eeg.shape[-1]:
            L = min(ecg.shape[-1], eeg.shape[-1])
            ecg = ecg[...,:L]
            eeg = eeg[...,:L]
            
        eeg = eeg-eeg.mean(axis=1, keepdims=True)
        eeg = mne.filter.notch_filter(eeg, Fs, 60, verbose=False)
        eeg = mne.filter.filter_data(eeg, Fs, 0.3, 35, verbose=False)
        
        ecg = ecg-ecg.mean()
        ecg = mne.filter.notch_filter(ecg, Fs, 60, verbose=False)
        ecg = mne.filter.filter_data(ecg, Fs, 0.3, 70 if 70<Fs/2 else None, verbose=False)
        
        flip = flip_ecg(ecg, Fs)
        #assert flip is not None, sid
        if type(flip)==bool and flip:
            ecg = -ecg
        
        with open(ann_path, 'r') as ff:
            annot = xml2dict.parse(ff.read())
        annot = pd.DataFrame(annot['PSGAnnotation']['ScoredEvents']['ScoredEvent'])
        annot['Start'] = annot.Start.astype(float)
        annot['Duration'] = annot.Duration.astype(float)
        #annot = annot[annot.EventType=='Stages|Stages'].reset_index(drop=True)
        
        ch_names = ['C3A2', 'C4A1', 'ECG']
        channel_info = [
                {'label': ch,
                 'dimension': 'uV',
                 'sample_rate': Fs,
                 'physical_max': 32767,
                 'physical_min': -32768,
                 'digital_max': 32767,
                 'digital_min': -32768,
                 'transducer': 'E',
                 'prefilter': ''}
            for ch in ch_names]
        with pyedflib.EdfWriter(output_path, len(ch_names), file_type=pyedflib.FILETYPE_EDFPLUS) as ff:
            ff.setSignalHeaders(channel_info)
            ff.setStartdatetime(start_time)
            ff.writeSamples(np.vstack([eeg, ecg]))
            for ii,r in annot.iterrows():
                ff.writeAnnotation(r.Start, r.Duration, r.EventConcept)

    df = df[['SID', 'Age', 'Fs', 'FsEEG', 'FsECG', 'CleanEDFPath']]
    print(df)
    df.to_excel('SOF_mastersheet-after-step2.xlsx', index=False)


if __name__ == '__main__':
    main()

