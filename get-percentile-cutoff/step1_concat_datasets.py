"""
Concatenate all datasets, harmonize columns, remove missing rows.

Usage:
    python step1_concat_datasets.py              # no APOE4, all 5 datasets
    python step1_concat_datasets.py --apoe4      # with APOE4 (binary), excludes ARIC

Output:
    dataset_all_original_BAI_noAPOE4.csv   or   dataset_all_original_BAI_withAPOE4.csv
"""

import argparse
import pandas as pd
import numpy as np

DATA_DIR = '..'


def load_mesa(apoe4=False):
    df = pd.read_csv(f'{DATA_DIR}/MESA/dataset_MESA.csv')
    df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)].reset_index(drop=True)
    df = df[df.race_Hispanic!=1].reset_index(drop=True)
    df['dataset'] = 'MESA'
    df['race_Asian'] = df['race_Chinese']   # map Chinese -> Asian
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'event', 'time2event', 'BAI']
    if apoe4:
        df['APOE4'] = (df['APOE4Count'] >= 1).astype(int)
        cols.append('APOE4')
    return df[cols]


def load_mros(apoe4=False):
    df = pd.read_csv(f'{DATA_DIR}/MrOS/dataset_MrOS.csv')
    df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)].reset_index(drop=True)
    df['dataset'] = 'MrOS'
    df['sexM'] = 1                                  # all-male study
    df['race_White'] = (1 - df['race_NonWhite']).astype(int)
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'event', 'time2event', 'BAI']
    if apoe4:
        df['APOE4'] = (df['APOE4Count'] >= 1).astype(int)
        cols.append('APOE4')
    return df[cols]


def load_sof(apoe4=False):
    df = pd.read_csv(f'{DATA_DIR}/SOF/dataset_SOF.csv')
    df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)].reset_index(drop=True)
    df['dataset'] = 'SOF'
    df['sexM'] = 0                                  # all-female study
    # Race encoding in SOF: 1 = White, 2 = Black
    df['race_White'] = (df['Race'] == 1).astype(int)
    df['race_Black'] = (df['Race'] == 2).astype(int)
    df['race_Asian'] = 0
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'event', 'time2event', 'BAI']
    if apoe4:
        cols.append('APOE4')   # SOF already has a binary APOE4 column
    return df[cols]


def load_aric():
    df = pd.read_csv(f'{DATA_DIR}/SHHS/ARIC/dataset_ARIC.csv')
    df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)].reset_index(drop=True)
    df['dataset'] = 'ARIC'
    df['race_Asian'] = 0    # no Asian participants in ARIC
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'event', 'time2event', 'BAI']
    return df[cols]


def load_fhs(apoe4=False):
    df = pd.read_csv(f'{DATA_DIR}/SHHS/FHS/dataset_FHS.csv')
    df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)].reset_index(drop=True)
    df['dataset'] = 'FHS'
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'event', 'time2event', 'BAI']
    if apoe4:
        df['APOE4'] = (df['APOE4Count'] >= 1).astype(int)
        cols.append('APOE4')
    return df[cols]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--apoe4', action='store_true',
                        help='Include APOE4 (binary) in output (excludes ARIC which lacks it)')
    args = parser.parse_args()

    dfs = [
        load_mesa(args.apoe4),
        load_mros(args.apoe4),
        load_sof(args.apoe4),
        load_fhs(args.apoe4),
    ]
    #if not args.apoe4:
    #    dfs.append(load_aric())

    combined = pd.concat(dfs, ignore_index=True, axis=0)
    n_before = len(combined)
    combined = combined.dropna()
    n_after = len(combined)

    suffix = '_withAPOE4' if args.apoe4 else '_noAPOE4'
    out_path = f'dataset_all_original_BAI{suffix}.csv'
    combined.to_csv(out_path, index=False)

    print(f"Rows before dropna: {n_before}, after: {n_after}")
    print(f"Dataset counts:\n{combined['dataset'].value_counts().to_string()}")
    print(f"\nSaved to: {out_path}")


if __name__ == '__main__':
    main()
