"""
Concatenate all datasets, harmonize columns, remove missing rows.

Output:
    dataset_all_original_BAI_noAPOE4.csv   or   dataset_all_original_BAI_withAPOE4.csv
"""

import pandas as pd
import numpy as np

DATA_DIR = '..'


def load_mesa(healthy=False):
    df = pd.read_csv(f'{DATA_DIR}/MESA/dataset_MESA_table1.csv')
    if healthy:
        df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)&(df.smoke_current==0)].reset_index(drop=True)
    df = df[df.race_Hispanic!=1].reset_index(drop=True)
    df['dataset'] = 'MESA'
    df['race_Asian'] = df['race_Chinese']   # map Chinese -> Asian
    df['APOE4'] = (df['APOE4Count'] >= 1).astype(int)
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'APOE4', 'event', 'time2event', 'BAI']
    return df[cols]


def load_mros(healthy=False):
    df = pd.read_csv(f'{DATA_DIR}/MrOS/dataset_MrOS_table1.csv')
    if healthy:
        df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)&(df.smoke_current==0)].reset_index(drop=True)
    df['dataset'] = 'MrOS'
    df['sexM'] = 1                                  # all-male study
    df['race_White'] = (1 - df['race_NonWhite']).astype(int)
    df['APOE4'] = (df['APOE4Count'] >= 1).astype(int)
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'APOE4', 'event', 'time2event', 'BAI']
    return df[cols]


def load_sof(healthy=False):
    df = pd.read_csv(f'{DATA_DIR}/SOF/dataset_SOF_table1.csv')
    if healthy:
        df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)&(df.smoke_current==0)].reset_index(drop=True)
    df['dataset'] = 'SOF'
    df['sexM'] = 0                                  # all-female study
    # Race encoding in SOF: 1 = White, 2 = Black
    df['race_White'] = (df['Race'] == 1).astype(int)
    df['race_Black'] = (df['Race'] == 2).astype(int)
    df['race_Asian'] = 0
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'APOE4', 'event', 'time2event', 'BAI']
    return df[cols]


def load_aric(healthy=False):
    df = pd.read_csv(f'{DATA_DIR}/SHHS/ARIC/dataset_ARIC_table1.csv')
    if healthy:
        df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)&(df.smoke_current==0)].reset_index(drop=True)
    df['dataset'] = 'ARIC'
    df['race_Asian'] = 0    # no Asian participants in ARIC
    df['APOE4'] = np.nan
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'APOE4', 'event', 'time2event', 'BAI']
    return df[cols]


def load_fhs(healthy=False):
    df = pd.read_csv(f'{DATA_DIR}/SHHS/FHS/dataset_FHS_table1.csv')
    if healthy:
        df = df[(df.sleepmed==0)&(df.stroke==0)&(df.heartattack==0)&(df.smoke_current==0)].reset_index(drop=True)
    df['dataset'] = 'FHS'
    df['APOE4'] = (df['APOE4Count'] >= 1).astype(int)
    cols = ['dataset', 'age', 'sexM', 'race_Asian', 'race_Black', 'race_White',
            'BMI', 'AHI', 'APOE4', 'event', 'time2event', 'BAI']
    return df[cols]


def main():

    dfs = [
        load_mesa(healthy=True),
        load_mros(healthy=True),
        load_sof(healthy=True),
        load_fhs(healthy=True),
        load_aric(healthy=True),
    ]
    combined = pd.concat(dfs, ignore_index=True, axis=0)

    suffix = '_noAPOE4'
    n_before = len(combined)
    combined = combined.dropna(subset=[x for x in combined.columns if x !='APOE4'])
    n_after = len(combined)
    out_path = f'dataset_ref_original_BAI{suffix}.csv'
    combined.to_csv(out_path, index=False)
    print(suffix)
    print(f"Rows before dropna: {n_before}, after: {n_after}")
    print(f"Dataset counts:\n{combined['dataset'].value_counts().to_string()}")
    print(f"Saved to: {out_path}")

    suffix = '_withAPOE4'
    n_before = len(combined)
    combined = combined.dropna()
    n_after = len(combined)
    out_path = f'dataset_ref_original_BAI{suffix}.csv'
    combined.to_csv(out_path, index=False)
    print(suffix)
    print(f"Rows before dropna: {n_before}, after: {n_after}")
    print(f"Dataset counts:\n{combined['dataset'].value_counts().to_string()}")
    print(f"Saved to: {out_path}")

if __name__ == '__main__':
    main()
