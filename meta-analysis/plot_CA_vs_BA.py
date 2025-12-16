import numpy as np
import pandas as pd
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn

plt.rcParams.update({
    #'font.family': 'Arial',
    'font.size': 9,
    'axes.linewidth': 0.5,
    #'axes.spines.top': False,
    #'axes.spines.right': False,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size': 2.5,
    'ytick.major.size': 2.5,
    #'legend.frameon': False,
    #'legend.fontsize': 6,
    #'figure.dpi': 300
})
seaborn.set_style('ticks')


def main():
    paths = {
    'MESA': '/data/haoqisun/BAI_dementia_community/MESA/dataset_MESA.csv',
    'ARIC': '/data/haoqisun/BAI_dementia_community/SHHS/ARIC/dataset_ARIC.csv',
    'FHS-OS': '/data/haoqisun/BAI_dementia_community/SHHS/FHS/dataset_FHS.csv',
    'MrOS': '/data/haoqisun/BAI_dementia_community/MrOS/dataset_MrOS.csv',
    'SOF': '/data/haoqisun/BAI_dementia_community/SOF/dataset_SOF.csv',
    }
    datasets = list(paths.keys())

    #bais = []
    #for i,dataset in enumerate(datasets):
    #    df = pd.read_csv(paths[dataset])
    #    bais.extend(df.BAI)
    #print(np.nanpercentile(bais, (0,1,2.5,5,10,90,95,97.5,99,100)))
    #[-34.75991851 -17.01709961 -14.17059486 -12.0287224   -9.65446022  5.56231041   7.50739092   9.45098299  11.33536682  19.16720285]
    clim = [-10,10]
    dataset2color = {'MESA':'m', 'ARIC':'g', 'FHS-OS':'b', 'MrOS':'gray', 'SOF':'orange'}

    cmap = matplotlib.cm.get_cmap('coolwarm')
    norm = matplotlib.colors.Normalize(vmin=clim[0], vmax=clim[1])
    plt.close()
    fig = plt.figure(figsize=(8.1*0.8,5*0.8))
    ax = fig.add_subplot(2,3,1)
    ax0 = ax

    lim = [30,100]
    lim2 = [35,95]
    all_ages = []
    all_bas = []
    all_datasets = []
    for i,dataset in enumerate(datasets+['All']):
        print(dataset)
        if i>0:
            ax = fig.add_subplot(2,3,i+1,sharex=ax0,sharey=ax0)
        ax.plot(lim2, lim2, c='k', ls='--', lw=0.5)
        if dataset!='All':
            df = pd.read_csv(paths[dataset])
            df = df.dropna(subset=['age','BAI'])
            age = df.age.values
            ba = df.age.values+df.BAI.values
            mae = np.mean(np.abs(ba-age))
            corr = pearsonr(age,ba)[0]
            c = [cmap(norm(x)) for x in df.BAI.values]
            ax.scatter(df.age, df.age+df.BAI, s=4, facecolor=c, edgecolor='none', alpha=0.5)
            ax.text(0.05, 0.97, dataset, ha='left', va='top', transform=ax.transAxes)
            ax.text(0.985, 0.03, f'N = {len(df):,}\nMAE = {mae:.1f} y\nCorr = {corr:.2f}', ha='right', va='bottom', transform=ax.transAxes)#, fontsize=11)
            
            # Add color bar to SOF panel (panel E, i=4)
            if dataset == 'MESA':
                cbar_ax = inset_axes(ax, width="3%", height="40%", loc='upper left',
                                   bbox_to_anchor=(0.018,-0.17,1,1), bbox_transform=ax.transAxes)
                cbar = plt.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
                                  cax=cbar_ax, orientation='vertical')
                cbar.set_label('BAI (y)', fontsize=7, labelpad=18, rotation=0, ha='right', va='center')
                cbar.ax.tick_params(labelsize=5, pad=1, width=0.5, length=2)
                cbar.set_ticks([-10, 0, 10], labels=['-10','0','+10'])
                # Fix PDF rendering: rasterize the colorbar content to prevent floating
                #cbar.solids.set_rasterized(True)
                #for axx in fig.get_axes():
                #    axx.set_rasterized(True)
            all_ages.extend(age)
            all_bas.extend(ba)
            all_datasets.extend([dataset]*len(age))
        else:
            all_ages = np.array(all_ages)
            all_bas = np.array(all_bas)
            all_datasets = np.array(all_datasets)
            for x in datasets:
                ids = all_datasets==x
                ax.scatter(all_ages[ids], all_bas[ids], s=4, facecolor=dataset2color[x], edgecolor='none', alpha=0.15)
                ax.scatter([-100],[0], s=18, facecolor=dataset2color[x], edgecolor='none', label=x)
            ax.legend(loc='upper left', fontsize=8, frameon=False, labelspacing=0.2, bbox_to_anchor=(-0.065,1.05), handletextpad=0.17)

        ax.text(-0.22+0.12*int(i%3>0), 1.04, chr(ord('A')+i), ha='right', va='top', transform=ax.transAxes, fontweight='bold')

        ax.set_xlim(lim)
        ax.set_ylim(lim)
        ticks = [40,60,80,100]
        ax.set_xticks(ticks)
        ax.set_yticks(ticks, labels=['40','60','80',''])

        if i//3==2-1:
            ax.set_xlabel('Age (Year)')
        if i%3==0:
            ax.set_ylabel('Brain Age (Year)')
        seaborn.despine()

    plt.tight_layout()
    #plt.show()
    plt.subplots_adjust(wspace=0.27)
    #plt.savefig('CA_vs_BA.png', bbox_inches='tight', dpi=300, facecolor='w')
    plt.savefig('CA_vs_BA.pdf')#, bbox_inches='tight', dpi=300, facecolor='w')


if __name__=='__main__':
    main()

