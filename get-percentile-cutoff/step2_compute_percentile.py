#!/usr/bin/env python3
"""
Step 2: For each individual, compute BAI percentile among 'similar' people.

Similar person definition:
    - age ± 5
    - same sex (sexM)
    - same race (race_Asian, race_Black, race_White all match)
    - BMI ± 3
    - AHI ± 3
    - same APOE4 (binary)  (only when --apoe4 flag is used)

Usage:
    python step2_compute_percentile.py              # reads dataset_all_original_BAI_noAPOE4.csv
    python step2_compute_percentile.py --apoe4      # reads dataset_all_original_BAI_withAPOE4.csv

Output:
    dataset_all_percentile_BAI_noAPOE4.csv   or   dataset_all_percentile_BAI_withAPOE4.csv
    Columns: dataset, BAIPercentile, n_similar, event, time2event
"""

import argparse
import numpy as np
import pandas as pd


def compute_percentiles(df, apoe4=False):
    ages      = df['age'].values.astype(float)
    sex       = df['sexM'].values.astype(float)
    race_A    = df['race_Asian'].values.astype(float)
    race_B    = df['race_Black'].values.astype(float)
    race_W    = df['race_White'].values.astype(float)
    bmis      = df['BMI'].values.astype(float)
    ahis      = df['AHI'].values.astype(float)
    bais      = df['BAI'].values.astype(float)
    if apoe4:
        apoe4v = df['APOE4'].values.astype(float)

    n = len(df)
    percentiles = np.full(n, np.nan)
    n_similars  = np.zeros(n, dtype=int)

    print(f"Computing percentiles for {n} individuals...")
    for i in range(n):
        if i % 500 == 0:
            print(f"  {i}/{n}", flush=True)

        mask = (
            (np.abs(ages - ages[i]) <= 5) &
            (sex == sex[i]) &
            (race_A == race_A[i]) &
            (race_B == race_B[i]) &
            (race_W == race_W[i]) &
            (np.abs(bmis - bmis[i]) <= 3) &
            (np.abs(ahis - ahis[i]) <= 3)
        )
        if apoe4:
            mask &= (apoe4v == apoe4v[i])

        similar_bais = bais[mask]
        n_sim = len(similar_bais)
        n_similars[i] = n_sim

        if n_sim > 0:
            # Percentile rank: fraction of similar people with BAI <= this person's BAI
            percentiles[i] = 100.0 * np.mean(similar_bais <= bais[i])

    return percentiles, n_similars


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--apoe4', action='store_true',
                        help='Use APOE4 in similarity matching')
    args = parser.parse_args()

    suffix = '_withAPOE4' if args.apoe4 else '_noAPOE4'
    in_path  = f'dataset_all_original_BAI{suffix}.csv'
    out_path = f'dataset_all_percentile_BAI{suffix}.csv'

    print(f"Reading {in_path}...")
    df = pd.read_csv(in_path)
    print(f"  {len(df)} rows, {len(df.columns)} columns")

    percentiles, n_similars = compute_percentiles(df, apoe4=args.apoe4)

    out = pd.DataFrame({
        'dataset':       df['dataset'].values,
        'BAIPercentile': percentiles,
        'n_similar':     n_similars,
        'event':         df['event'].values,
        'time2event':    df['time2event'].values,
    })

    out.to_csv(out_path, index=False)
    print(f"\nSaved {len(out)} rows to {out_path}")
    print(f"n_similar summary:\n{pd.Series(n_similars).describe().to_string()}")
    print(f"Rows with n_similar >= 10: {(n_similars >= 10).sum()}")


if __name__ == '__main__':
    main()
