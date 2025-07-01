import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from sklearn.metrics import auc
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import seaborn
seaborn.set_style('ticks')


def aggregate_roc(fprs, tprs, method=np.mean):
    x = np.concatenate(fprs)
    fpr = np.linspace(x.min(), x.max(), 100)
    assert len(fprs)==len(tprs)
    tpr_ = []
    for f,t in zip(fprs, tprs):
        assert np.isnan(f).sum()==0 and np.isnan(t).sum()==0
        func = interp1d(f,t)
        tpr_.append( func(fpr) )
    tpr = method(np.array(tpr_), axis=0)
    return fpr, tpr


def main():
    paths = {
    'MESA': '/data/haoqisun/BAI_dementia_community/MESA/ROC_results_MESA_nbt100.xlsx',
    'ARIC': '/data/haoqisun/BAI_dementia_community/SHHS/ARIC/ROC_results_ARIC_nbt100.xlsx',
    'FHS': '/data/haoqisun/BAI_dementia_community/SHHS/FHS/ROC_results_FHS_nbt100.xlsx',
    'MrOS': '/data/haoqisun/BAI_dementia_community/MrOS/ROC_results_MrOS_nbt100.xlsx',
    'SOF': '/data/haoqisun/BAI_dementia_community/SOF/ROC_results_SOF_nbt100.xlsx',
    }

    """
    dfs = []
    for k,p in paths.items():
        df = pd.read_csv(p)
        df.insert(0, 'Dataset', k)
        dfs.append(df)
    df = pd.concat(dfs, axis=0, ignore_index=True)
    print(df)
    df.to_csv('sens_spec.csv', index=False)
    """
    nbt = 100

    plt.close()
    fig = plt.figure(figsize=(6,6))
    ax = fig.add_subplot(111)
    ax.plot([0,1],[0,1],c='k',lw=1,ls='--')
    alpha = 0.3

    fprs_cov = []; tprs_cov = []; aucs_cov = []
    fprs_bai_cov = []; tprs_bai_cov = []; aucs_bai_cov = []
    for k,p in paths.items():
        with pd.ExcelFile(p) as f:
            df_auc = pd.read_excel(f, 'AUC', index_col=0)
            #df_bai_roc = pd.read_excel(f, 'BAI-only-ROC')
            #df_bai_roc_bt = pd.read_excel(f, f'BAI-only-ROC-bt{nbt}')
            df_cov_roc = pd.read_excel(f, 'cov-only-ROC')
            #df_cov_roc_bt = pd.read_excel(f, f'cov-only-ROC-bt{nbt}')
            df_bai_cov_roc = pd.read_excel(f, 'BAI_cov-ROC')
            #df_bai_cov_roc_bt = pd.read_excel(f, f'BAI_cov-ROC-bt{nbt}')
        df_cov_roc = df_cov_roc[['sens','specA']].dropna().drop_duplicates(ignore_index=True)
        df_bai_cov_roc = df_bai_cov_roc[['sens','specA']].dropna().drop_duplicates(ignore_index=True)
        ax.plot(1-df_cov_roc.specA, df_cov_roc.sens, c='r', alpha=alpha, lw=1)
        ax.plot(1-df_bai_cov_roc.specA, df_bai_cov_roc.sens, c='b', alpha=alpha, lw=1)
        aucs_cov.append(df_auc.loc['perf1','AUC'])
        fprs_cov.append(1-df_cov_roc.specA)
        tprs_cov.append(df_cov_roc.sens)
        aucs_bai_cov.append(df_auc.loc['perf2', 'AUC'])
        fprs_bai_cov.append(1-df_bai_cov_roc.specA)
        tprs_bai_cov.append(df_bai_cov_roc.sens)

    method = np.median
    fpr_cov, tpr_cov = aggregate_roc(fprs_cov, tprs_cov, method=method)
    fpr_bai_cov, tpr_bai_cov = aggregate_roc(fprs_bai_cov, tprs_bai_cov, method=method)
    #auc_cov = method(aucs_cov)
    auc_cov = auc(fpr_cov, tpr_cov)
    auc_cov_min = np.min(aucs_cov)
    auc_cov_max = np.max(aucs_cov)
    #auc_bai_cov = method(aucs_bai_cov)
    auc_bai_cov = auc(fpr_bai_cov, tpr_bai_cov)
    auc_bai_cov_min = np.min(aucs_bai_cov)
    auc_bai_cov_max = np.max(aucs_bai_cov)
    ax.plot(fpr_cov, tpr_cov, c='r', alpha=1, lw=2, label=f'Covariates only: AUC= {auc_cov:.3f}, range= {auc_cov_min:.3f}-{auc_cov_max:.3f}')
    ax.plot(fpr_bai_cov, tpr_bai_cov, c='b', alpha=1, lw=2, label=f'BAI + Covariates: AUC= {auc_bai_cov:.3f}, range= {auc_bai_cov_min:.3f}-{auc_bai_cov_max:.3f}')
    ax.legend(frameon=True, loc='lower right', fontsize=11)

    ax.set_xlim(-0.01, 1.01)
    ax.set_ylim(-0.01, 1.01)
    ax.set_xlabel('FPR (1-specificity)')
    ax.set_ylabel('TPR (sensitivity)')
    ax.grid(True)
    seaborn.despine()

    plt.tight_layout()
    #plt.show()
    plt.savefig('ROC_results.png', bbox_inches='tight')


if __name__=='__main__':
    main()


