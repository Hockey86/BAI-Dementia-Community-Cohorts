"""
For each individual, compute BAI percentile among 'similar' healthy reference people.

Similar person definition:
    - diff(age)+diff(BMI)+diff(AHI) <= 10
    - same sex (sexM)
    - same race (race_Asian, race_Black, race_White all match)
    - same APOE4 (binary)  (only when --apoe4 flag is used)

Output:
    (no flag)                       dataset_all_percentile_BAI_noAPOE4.csv
    --apoe4                         dataset_all_percentile_BAI_withAPOE4.csv
    --without-apoe4-same-dataset    dataset_all_percentile_BAI_noAPOE4_sameDataset.csv
    Columns: dataset, BAIPercentile, n_similar, event, time2event

The default (no flag) run uses a larger reference pool and keeps cohorts without
APOE4 genotyping (e.g. ARIC), so it is not comparable to --apoe4. Use
--without-apoe4-same-dataset to isolate the effect of APOE4 matching alone.
"""

import argparse, pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from step1_generate_healthy_ref_dataset import *


def compute_percentiles(df, df_ref, min_n_match, apoe4=False):
    datasets  = df_ref['dataset'].values
    ages      = df_ref['age'].values.astype(float)
    sexs      = df_ref['sexM'].values.astype(float)
    races_A   = df_ref['race_Asian'].values.astype(float)
    races_B   = df_ref['race_Black'].values.astype(float)
    races_W   = df_ref['race_White'].values.astype(float)
    bmis      = df_ref['BMI'].values.astype(float)
    log_ahis  = df_ref['log_AHI'].values.astype(float)
    bais      = df_ref['BAI'].values.astype(float)
    if apoe4:
        apoe4v = df_ref['APOE4'].values.astype(float)

    n = len(df)
    percentiles = np.full(n, np.nan)
    n_similars  = np.zeros(n, dtype=int)

    print(f"Computing percentiles for {n} individuals...")
    for i in tqdm(range(n)):
        #    (np.abs(ages - df.age.iloc[i]) <= 5) &
        #    (np.abs(bmis - df.BMI.iloc[i]) <= 5) &
        #    (np.abs(ahis - df.AHI.iloc[i]) <= 5) &
        mask = (
            (np.abs(ages - df.age.iloc[i]) + np.abs(bmis - df.BMI.iloc[i]) + np.abs(log_ahis - df.log_AHI.iloc[i]) <= 3) &
            (sexs == df.sexM.iloc[i]) &
            (races_A == df.race_Asian.iloc[i]) &
            (races_B == df.race_Black.iloc[i]) &
            (races_W == df.race_White.iloc[i]) &
            (datasets== df.dataset.iloc[i])
        )
        if apoe4:
            mask &= (apoe4v == df.APOE4.iloc[i])

        diffs = (np.abs(ages - df.age.iloc[i]) +
                 np.abs(bmis - df.BMI.iloc[i]) +
                 np.abs(log_ahis - df.log_AHI.iloc[i]))

        similar_bais = bais[mask]
        similar_diffs = diffs[mask]
        n_sim = len(similar_bais)
        n_similars[i] = n_sim

        # Gaussian weights: exp(-diff^2/sigma), sigma chosen so weight=0.01 when diff=X
        # exp(-X^2/sigma) = 0.01  =>  sigma = X^2 / ln(100)
        sigma = 3**2 / np.log(100.0)
        #print(f'{sigma = }')

        if n_sim >= min_n_match:
            weights = np.exp(-similar_diffs**2 / sigma)
            w_total = weights.sum()
            w_below = weights[similar_bais <= df.BAI.iloc[i]].sum()
            percentiles[i] = 100.0 * w_below / w_total

    return percentiles, n_similars


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--apoe4', action='store_true',
                        help='Use APOE4 in similarity matching')
    group.add_argument('--without-apoe4-same-dataset', action='store_true',
                        help='Do not match on APOE4, but use the same reference pool and '
                             'the same individuals as --apoe4, so the two runs differ only '
                             'in whether APOE4 is a matching variable')
    args = parser.parse_args()

    same_dataset = args.without_apoe4_same_dataset
    # same_dataset borrows the withAPOE4 reference pool but does not match on APOE4
    ref_suffix = '_withAPOE4' if (args.apoe4 or same_dataset) else '_noAPOE4'
    out_suffix = '_noAPOE4_sameDataset' if same_dataset else ref_suffix
    ref_path = f'dataset_ref_original_BAI{ref_suffix}.csv'
    out_path = f'dataset_all_percentile_BAI{out_suffix}.csv'

    df_ref = pd.read_csv(ref_path)

    dfs = [
        load_mesa(),
        load_mros(),
        load_sof(),
        load_fhs(),
        load_aric(),
    ]
    df = pd.concat(dfs, ignore_index=True, axis=0)
    df = df[df.time2event>0].reset_index(drop=True)

    # Restrict to individuals with known APOE4 so both APOE4 runs start from the same set.
    # Under --apoe4 these rows can never match anyway (NaN==x is always False), so this
    # only makes the exclusion explicit; under same_dataset it is what makes it fair.
    if args.apoe4 or same_dataset:
        n_before = len(df)
        df = df[df.APOE4.notna()].reset_index(drop=True)
        print(f"Restricted to known APOE4: {n_before} -> {len(df)}")
    print(df.shape)

    # standardize based on the healthy reference subset, and be cohort specific
    cohorts = df.dataset.unique()
    df_ref['log_AHI'] = np.log1p(df_ref.AHI)
    df['log_AHI'] = np.log1p(df.AHI)
    cols = ['age', 'BMI', 'log_AHI']
    means = {}; stds = {}
    for coh in cohorts:
        means[coh] = df_ref.loc[df_ref.dataset==coh, cols].mean().values
        stds[coh] = df_ref.loc[df_ref.dataset==coh, cols].std().values
        for col in cols:
            df[col] = df[col].astype(float)
            df_ref[col] = df_ref[col].astype(float)
        df.loc[df.dataset==coh, cols] = (df.loc[df.dataset==coh, cols].values-means[coh])/stds[coh]
        df_ref.loc[df_ref.dataset==coh, cols] = (df_ref.loc[df_ref.dataset==coh, cols].values-means[coh])/stds[coh]

    min_n_match = 0
    percentiles, n_similars = compute_percentiles(df, df_ref, min_n_match, apoe4=args.apoe4)

    out = pd.DataFrame({
        'dataset':       df['dataset'].values,
        'BAI':       df['BAI'].values,
        'BAIPercentile': percentiles,
        'n_similar':     n_similars,
        'event':         df['event'].values,
        'time2event':    df['time2event'].values,
    })
    print(f"n_similar summary:\n{pd.Series(n_similars).describe().to_string()}")

    min_n_match = 10
    print(f"Rows with n_similar >= {min_n_match}: {(n_similars >= min_n_match).sum()} ({(n_similars>=min_n_match).mean()*100.:.1f}%)")
    print(f"Rows with n_similar < {min_n_match}: {(n_similars < min_n_match).sum()} ({(n_similars<min_n_match).mean()*100.:.1f}%)")
    out = out[out.n_similar>=min_n_match].reset_index(drop=True)
    out.to_csv(out_path, index=False)

    with open(f'means_stds_ref_cohorts{out_suffix}.pickle', 'wb') as ff:
        pickle.dump({'means':means, 'stds':stds}, ff)


if __name__ == '__main__':
    main()
