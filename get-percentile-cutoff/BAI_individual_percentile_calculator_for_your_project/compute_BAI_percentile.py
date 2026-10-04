"""
Compute each person's brain age index (BAI) percentile among similar healthy
reference people.

Similar person:
    - same sex and race
    - |age| + |BMI| + |log(1+AHI)| distance <= 3, each standardized by the reference SD
    - closer reference people get more weight (Gaussian weight on the distance)

Usage:
    python compute_BAI_percentile.py my_data.csv my_data_with_percentile.csv

Input CSV columns (one row per person):
    BAI          brain age index in years (brain age minus chronological age)
    age          years
    sexM         1 = male, 0 = female
    race_Asian, race_Black, race_White   1/0 each (all 0 = other race)
    BMI          kg/m^2
    AHI          apnea-hypopnea index, events/hour
Other columns are kept unchanged in the output.

Output adds two columns:
    BAIPercentile   0-100; NaN if fewer than MIN_N_SIMILAR similar reference people
    n_similar       number of similar reference people
"""

import argparse, os
import numpy as np
import pandas as pd

MAX_DIST      = 3.0   # distance threshold, in reference SD units
MIN_N_SIMILAR = 10    # minimum number of similar reference people to report a percentile


def main():
    parser = argparse.ArgumentParser(description='Compute BAI percentiles among similar healthy reference people.')
    parser.add_argument('input_csv')
    parser.add_argument('output_csv')
    args = parser.parse_args()

    ref = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dataset_healthy_reference.csv'))
    df  = pd.read_csv(args.input_csv)

    exact_cols = ['sexM', 'race_Asian', 'race_Black', 'race_White']
    missing = [c for c in ['BAI', 'age', 'BMI', 'AHI'] + exact_cols if c not in df.columns]
    if missing:
        raise ValueError(f'Input is missing columns: {missing}')

    # standardize age, BMI, log(1+AHI) using the reference mean and SD
    dist_cols = ['age', 'BMI', 'log_AHI']
    ref['log_AHI'] = np.log1p(ref.AHI)
    x_new = np.column_stack([df.age, df.BMI, np.log1p(df.AHI)]).astype(float)
    mean, std = ref[dist_cols].mean().values, ref[dist_cols].std().values
    x_ref = (ref[dist_cols].values - mean) / std
    x_new = (x_new - mean) / std
    exact_ref = ref[exact_cols].values.astype(float)
    exact_new = df[exact_cols].values.astype(float)
    bai_ref   = ref.BAI.values
    bai_new   = df.BAI.values.astype(float)

    # Gaussian weight, equal to 0.01 at the distance threshold
    sigma = MAX_DIST**2 / np.log(100.0)

    percentiles = np.full(len(df), np.nan)
    n_similar   = np.zeros(len(df), dtype=int)
    for i in range(len(df)):
        dist = np.abs(x_ref - x_new[i]).sum(axis=1)
        mask = (dist <= MAX_DIST) & (exact_ref == exact_new[i]).all(axis=1)
        n_similar[i] = mask.sum()
        if n_similar[i] >= MIN_N_SIMILAR and not np.isnan(bai_new[i]):
            w = np.exp(-dist[mask]**2 / sigma)
            percentiles[i] = 100.0 * w[bai_ref[mask] <= bai_new[i]].sum() / w.sum()

    df['BAIPercentile'] = percentiles
    df['n_similar']     = n_similar
    df.to_csv(args.output_csv, index=False)
    print(f'{np.isfinite(percentiles).sum()} of {len(df)} people got a percentile '
          f'(the rest had < {MIN_N_SIMILAR} similar reference people or missing values).')
    print(f'Saved to {args.output_csv}')


if __name__ == '__main__':
    main()
