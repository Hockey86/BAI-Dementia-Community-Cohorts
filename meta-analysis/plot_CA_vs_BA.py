import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import seaborn
seaborn.set_style('ticks')


def main():
    paths = {
    'MESA': '/data/haoqisun/BAI_dementia_community/MESA/dataset_MESA.csv',
    'ARIC': '/data/haoqisun/BAI_dementia_community/SHHS/ARIC/dataset_ARIC.csv',
    'FHS-OFF1': '/data/haoqisun/BAI_dementia_community/SHHS/FHS/dataset_FHS.csv',
    'MrOS': '/data/haoqisun/BAI_dementia_community/MrOS/dataset_MrOS.csv',
    'SOF': '/data/haoqisun/BAI_dementia_community/SOF/dataset_SOF.csv',
    }

    plt.close()
    fig = plt.figure(figsize=(8.1,5))

    lim = [30,100]
    ticks = [40,60,80,100]
    datasets = list(paths.keys())
    for i,dataset in enumerate(datasets):
        print(dataset)
        df = pd.read_csv(paths[dataset])
        df = df.dropna(subset=['age','BAI'])
        age = df.age.values
        ba = df.age.values+df.BAI.values
        mae = np.mean(np.abs(ba-age))
        corr = pearsonr(age,ba)[0]
        if i==0:
            ax = fig.add_subplot(2,3,i+1)
            ax0 = ax
        else:
            ax = fig.add_subplot(2,3,i+1,sharex=ax0,sharey=ax0)
        ax.plot(lim, lim, c='r', ls='--')
        ax.scatter(df.age, df.age+df.BAI, s=2, c='k')
        ax.text(0.05, 0.99, dataset, ha='left', va='top', transform=ax.transAxes)
        ax.text(0.98, 0.01, f'MAE {mae:.1f}y\nCorr {corr:.2f}', ha='right', va='bottom', transform=ax.transAxes, fontsize=11)
        ax.text(-0.35, 1.06, chr(ord('a')+i)+'.', ha='right', va='top', transform=ax.transAxes, fontweight='bold')

        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)
        ax.set_xlabel('Age (Year)')
        ax.set_ylabel('Brain Age (Year)')
        seaborn.despine()

    plt.tight_layout()
    #plt.show()
    plt.savefig('CA_vs_BA.png', bbox_inches='tight')


if __name__=='__main__':
    main()

