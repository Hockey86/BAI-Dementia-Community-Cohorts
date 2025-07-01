import os
import sys
from collections import Counter
from itertools import groupby
import datetime
import numpy as np
import pandas as pd
import scipy.io as sio
from tqdm import tqdm
import joblib
import matplotlib
matplotlib.rc('pdf', fonttype=42)
matplotlib.rc('font', **{'size': 14})
import matplotlib.pyplot as plt
import seaborn
seaborn.set_style('ticks')
from myfunctions import delete_nan_start_end, load_mros_sof_data
from segment_EEG import segment_EEG
from extract_features_parallel import extract_features


def myprint(epoch_status):
    sm = Counter(epoch_status)
    for k, v in sm.items():
        print(f'{k}: {v}/{len(epoch_status)}, {v*100./len(epoch_status)}%')


if __name__=='__main__':
    epoch_length = 30 # [s]
    amplitude_thres = 500 # [uV]
    line_freq = 60.  # [Hz]
    bandpass_freq = [0.5, 20.]  # [Hz]
    n_jobs = 16
    clean_only = True
    newFs = 200.
    stages = ['W','N1','N2','N3','R']
    stage2num = {'W':5,'R':4,'N1':3,'N2':2,'N3':1}
    num2stage = {stage2num[x]:x for x in stage2num}
    minimum_epochs_per_stage = 2

    # get list of files
    dataset = 'MrOS_SOF'
    df = pd.read_csv(f'../mastersheet_{dataset}.csv')
            
    # define output folder
    output_feature_dir = f'features_{dataset}'
    os.makedirs(output_feature_dir, exist_ok=True)
    output_sleep_vis_dir = f'sleep_spectrogram_hypnogram_{dataset}'
    os.makedirs(output_sleep_vis_dir, exist_ok=True)

    # for each recording
    for si in tqdm(range(len(df))):
        sid = df.SID.iloc[si]
        signal_path = df.EDFPath.iloc[si]
        annot_path = df.AnnotPath.iloc[si]
        feature_path1 = os.path.join(output_feature_dir, f'features_{sid}.mat')
        feature_path2 = os.path.join(output_feature_dir, f'features_{sid}.csv')
        figure_path = os.path.join(output_sleep_vis_dir, sid+'.png')
        if os.path.exists(feature_path1) and os.path.exists(feature_path2):
            continue
        
        # load dataset
        try:
            if dataset=='MrOS_SOF':
                EEG, sleep_stages, df_event, params = load_mros_sof_data(signal_path, annot_path)
            else:
                raise NotImplementedError(f'Unknown dataset={dataset}')
        except Exception as ee:
            print(ee)
            continue
        Fs = params['Fs']
        start_time = params['start_time']
        ch_names = params['ch_names']
        EEG, sleep_stages, start_time = delete_nan_start_end(EEG, sleep_stages, start_time, Fs)
        combined_EEG_channels = ['C']
        combined_EEG_channels_ids = [[0,1]]

        # segment EEG
        epochs, sleep_stages, epoch_start_idx, epoch_status, specs, freq, qs = segment_EEG(EEG, sleep_stages, epoch_length, epoch_length, Fs, newFs, notch_freq=line_freq, bandpass_freq=bandpass_freq, amplitude_thres=amplitude_thres, n_jobs=n_jobs)
        if epochs.shape[0] <= 0:
            raise ValueError('Empty EEG segments')
        Fs = newFs
        specs_db = 10*np.log10(specs)

        # plot spectrogram and sleep stages
        plt.close()
        fig = plt.figure(figsize=(13.8,6.4))
        gs = fig.add_gridspec(1+len(combined_EEG_channels), 1, height_ratios=[1]+[3]*len(combined_EEG_channels))

        tt = np.arange(len(sleep_stages))*epoch_length/3600
        xticks = np.arange(0, np.floor(tt.max())+1) 
        xticklabels = [(start_time+datetime.timedelta(hours=x)).strftime('%H:%M') for x in xticks]
        
        ax_ss = fig.add_subplot(gs[0])
        ax_ss.step(tt, sleep_stages, color='k', where='post')
        ax_ss.yaxis.grid(True)
        ax_ss.set_ylim([0.7,5.3])
        ax_ss.set_yticks([1,2,3,4,5])
        ax_ss.set_yticklabels(['N3', 'N2', 'N1', 'R', 'W'])
        ax_ss.set_xlim([tt.min(), tt.max()])
        seaborn.despine()
        plt.setp(ax_ss.get_xlabel(), visible=False)
        plt.setp(ax_ss.get_xticklabels(), visible=False)
        
        for chi in range(len(combined_EEG_channels)):
            ax_spec = fig.add_subplot(gs[1+chi], sharex=ax_ss)
            specs_db_ch = (specs_db[:,chi*2]+specs_db[:,chi*2+1])/2
            ax_spec.imshow(
                    specs_db_ch.T, aspect='auto', origin='lower', cmap='turbo',
                    vmin=-5, vmax=15,
                    extent=(tt.min(), tt.max(), freq.min(), freq.max()))
            ax_spec.set_ylabel(f'Avg\n{combined_EEG_channels[chi]}', rotation=0, ha='right')
            #ax_spec.set_ylabel('freq (Hz)')
            ax_spec.set_xticks(xticks)
            ax_spec.set_xticklabels(xticklabels)
            if chi<len(combined_EEG_channels)-1:
                plt.setp(ax_spec.get_xlabel(), visible=False)
                plt.setp(ax_spec.get_xticklabels(), visible=False)
        
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.11)
        plt.savefig(figure_path, bbox_inches='tight', pad_inches=0.03)

        if clean_only:
            good_ids = np.where(epoch_status=='clean')[0]
            if len(good_ids)<=300:
                myprint(epoch_status)
                raise ValueError('<=300 clean epochs')
            epochs = epochs[good_ids]
            specs = specs[good_ids]
            sleep_stages = sleep_stages[good_ids]
            epoch_start_idx = epoch_start_idx[good_ids]

        # extract brain age features

        # normalize signal
        nch = epochs.shape[1]
        q1, q2, q3 = qs
        epochs = (epochs - q2.reshape(1,nch,1)) / (q3.reshape(1,nch,1)-q1.reshape(1,nch,1))

        features, feature_names = extract_features(
            epochs, Fs, ch_names, 2,
            2, 1, return_feature_names=True,
            combined_channel_names=combined_EEG_channels,
            n_jobs=n_jobs, verbose=False)
        artifact_ratio = 1-len(sleep_stages)/len(epoch_status)
        num_missing_stage = 5-len(set(sleep_stages[~np.isnan(sleep_stages)]))
        
        #myprint(epoch_status)
        sio.savemat(feature_path1, {
            'start_time':start_time.strftime('%Y-%m-%d %H:%M:%S'),
            #'EEG_feature_names':feature_names,
            #'EEG_features':features,
            'ch_names':ch_names,
            'combined_EEG_channels':combined_EEG_channels,
            'combined_EEG_channels_ids':combined_EEG_channels_ids,
            'EEG_specs':specs,
            'EEG_frequency':freq,
            'sleep_stages':sleep_stages,
            'epoch_start_idx':epoch_start_idx,
            #'age':age,
            #'gender':df.Gender.iloc[si],
            'epoch_status':epoch_status,
            'Fs':Fs,
            'artifact_ratio':artifact_ratio,
            'num_missing_stage':num_missing_stage,
            })
                
        # log-transform brain age features
        features_no_log = np.array(features)
        features = np.sign(features)*np.log1p(np.abs(features))
        
        # average features across sleep stages
        X = []
        X_no_log = []
        for stage in stages:
            ids = sleep_stages==stage2num[stage]
            if ids.sum()>=minimum_epochs_per_stage:
                X.append(np.nanmean(features[ids], axis=0))
                X_no_log.append(np.nanmean(features_no_log[ids], axis=0))
            else:
                X.append(np.zeros(features.shape[1])+np.nan)
                X_no_log.append(np.zeros(features.shape[1])+np.nan)
        X = np.concatenate(X)
        X_no_log = np.concatenate(X_no_log)
        
        cols = np.concatenate([[x.strip()+'_'+stage for x in feature_names] for stage in stages])
        df_feat = pd.DataFrame(data=X.reshape(1,-1), columns=cols)
        df_feat.insert(0, 'SID', sid)
        df_feat.to_csv(feature_path2, index=False)

        df_feat = pd.DataFrame(data=X_no_log.reshape(1,-1), columns=cols)
        df_feat.insert(0, 'SID', sid)
        df_feat.to_csv(feature_path2.replace('.csv','_no_log.csv'), index=False)
        
