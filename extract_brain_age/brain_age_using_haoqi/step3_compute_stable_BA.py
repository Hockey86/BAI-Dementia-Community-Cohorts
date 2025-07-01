import os
import sys
import pickle
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.special import logsumexp
from tqdm import tqdm
from sklearn.neighbors import KernelDensity


def get_perc(BAI, age, sex, bmi, df_ref, ref_bai_col, plot=False, fig_path=None):
    if pd.notna(sex):
        df_ref = df_ref[df_ref.Sex==sex].reset_index(drop=True)
    else:
        df_ref = df_ref.copy()

    # similarity = exp(-(x-y)^2/sigma)

    # find sigma so that at bound, similarity=0.05
    sigma_age = -(10**2)/np.log(0.05)
    sigma_bmi = -(5**2)/np.log(0.05)

    cols = []
    log_sim = []
    if pd.notna(age):
        log_sim.append(-(df_ref.Age.values-age)**2/sigma_age)
    if pd.notna(bmi):
        log_sim.append(-(df_ref.BMI.values-bmi)**2/sigma_bmi)
    sim = np.exp(sum(log_sim)/len(log_sim))

    bai_ref = df_ref[ref_bai_col].values
    age_bound = np.abs(bai_ref).max()

    std = np.median(np.abs(bai_ref-np.median(bai_ref)))*1.4826
    kde = KernelDensity(bandwidth=std*0.4)
    kde.fit(bai_ref.reshape(-1,1), sample_weight=sim)
    bai_grid = np.linspace(-age_bound-10,age_bound+10,10000)
    bai_log_prob = kde.score_samples(bai_grid.reshape(-1,1))

    perc = np.exp(logsumexp(bai_log_prob[bai_grid<BAI]) - logsumexp(bai_log_prob))*100

    if plot:
        plt.close()
        fig = plt.figure(figsize=(10,8))

        ax = fig.add_subplot(111)
        res = ax.hist(bai_ref, bins=30, density=True, weights=sim)
        ax.plot(bai_grid, np.exp(bai_log_prob), c='k', lw=2)
        ax.plot([BAI]*2, [0,res[0].max()*1.01], ls='--', c='r', lw=2)
        ax.text(BAI, res[0].max()*1.01, f'{perc:.0f}%', ha='center', va='bottom', color='r', fontweight='bold')
        ax.set_xlabel('BAI (year)')
        ax.set_ylabel('Probability density')
        ax.set_xlim(-age_bound-1, age_bound+1)
        seaborn.despine()

        plt.tight_layout()
        if fig_path is None:
            plt.show()
        else:
            plt.savefig(fig_path)

    return perc, (bai_grid, bai_log_prob)


def main():
    # get list of subjects
    df = pd.read_csv('../mastersheet_MrOS_SOF.csv')
    dataset = 'MrOS'
    df = df[df.Dataset==dataset].reset_index(drop=True)
    df2 = pd.read_csv(f'../dataset_{dataset}.csv')
    df2['SID'] = df2.id.str.lower()
    df2 = df2.rename(columns={'vsage1':'Age', 'hwbmi':'BMI'})
    df2['Sex'] = 1
    df = df[['SID']].merge(df2[['SID', 'Age', 'Sex', 'BMI']], on='SID', how='inner', validate='1:1')

    # load brain age model
    brain_age_dir = 'brain_age_model_c'
    feature_names = list(pd.read_csv(os.path.join(brain_age_dir, 'BA_features_used.csv')).feature.values)
    sys.path.insert(0, brain_age_dir)
    with open(os.path.join(brain_age_dir, 'stable_brain_age_model.pickle'), 'rb') as ff:
        model = pickle.load(ff)
    df_ref = pd.read_csv('BAI_MGH_healthy.csv')

    # read features
    feature_dir = 'features_MrOS_SOF'
    for i in tqdm(range(len(df))):
        sid = df.SID.iloc[i]
        try:
            df_feat = pd.read_csv(os.path.join(feature_dir, f'features_{sid}_no_log.csv'))
            df_sp = pd.read_csv(os.path.join(feature_dir, f'spindle_features_{sid}.csv'))
            mat = sio.loadmat(os.path.join(feature_dir, f'features_{sid}.mat'),
                    variable_names=['ch_names', 'combined_EEG_channels', 'combined_EEG_channels_ids',
                    'artifact_ratio', 'num_missing_stage'])
        except Exception as ee:
            print(str(ee))
            continue
        combined_EEG_channels = np.char.strip(mat['combined_EEG_channels'].flatten())
        combined_EEG_channels_ids = mat['combined_EEG_channels_ids']
        EEG_channels = np.char.strip(mat['ch_names'].flatten())
        artifact_ratio = float(mat['artifact_ratio'])
        num_missing_stage = int(mat['num_missing_stage'])
        
        # match the feature names in model
        df_feat = df_feat.rename(columns={x:x.replace('/','_') for x in df_feat.columns if '/' in x})
        
        stages = ['W','N1','N2','N3','R']
        for c1, c2 in zip(combined_EEG_channels_ids, combined_EEG_channels):
            for stage in stages:
                df_feat[f'kurtosis_{stage}_{c2}'] = (df_feat[f'kurtosis_{EEG_channels[c1[0]]}_{stage}']+df_feat[f'kurtosis_{EEG_channels[c1[1]]}_{stage}'])/2
            
        
        # compute brain age
        # model contains all preprocessing and adjustment steps
        X = df_feat.merge(df_sp, on='SID', how='inner')[feature_names].values

        # COUPL_OVERLAP has different mean, remove that feature
        #for chn in combined_EEG_channels:
        #    if 'COUPL_OVERLAP_'+chn in feature_names:
        #        idx = feature_names.index('COUPL_OVERLAP_'+chn)
        #        X[:,idx] = model.steps[0][1].mean_[idx]

        age = df.Age.iloc[i]
        sex = df.Sex.iloc[i]
        bmi = df.BMI.iloc[i]
        if type(sex)==str:
            if sex.lower()=='m':
                sex = 1
            elif sex.lower()=='f':
                sex = 0
        if sex not in [0,1]:
            print(f'Unknown sex encoding (M or 1, F or 0): {sex}\nIgnore it.')
            sex = np.nan

        age2 = 70 if pd.isna(age) else age
        BA = model.predict(X, y=age2)[0]

        # get BAI
        BAI = BA-age

        # get BAI percent
        BAI_perc = np.nan
        if pd.notna(BAI) and (pd.notna(age)|pd.notna(sex)|pd.notna(bmi)):
            try:
                BAI_perc, hist = get_perc(BAI, age, sex, bmi, df_ref, 'BAI')
            except Exception as ee:
                BAI_perc = np.nan
        
        df.loc[i, 'RobustBA'] = BA
        df.loc[i, 'RobustBAI'] = BAI
        df.loc[i, 'RobustBAIPerc'] = BAI_perc
        df.loc[i, 'artifact_ratio'] = artifact_ratio
        df.loc[i, 'num_missing_stage'] = num_missing_stage

        if i%10==9 or i==len(df)-1:
            print(df)
            df.to_csv(f'robust_BA_{dataset}.csv', index=False)
    

if __name__ == '__main__':
    main()

