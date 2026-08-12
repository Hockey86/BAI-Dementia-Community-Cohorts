import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 20,
    'axes.titlesize': 20,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 13,
    'legend.title_fontsize': 12,
    'axes.linewidth': 1.5,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'lines.linewidth': 2.2,
})
import seaborn as sns
sns.set_style('ticks')
from sksurv.nonparametric import cumulative_incidence_competing_risks


def fit_cif(durations, events):
    """Aalen-Johansen CIF for event=1 (dementia) using scikit-survival.

    Returns arrays prepended with (t=0, CIF=0) so the step plot starts at origin.
    events: 0=censored, 1=dementia, 2=death.

    sksurv returns conf[cause, bound, time] where for cause 1 (dementia):
      conf[1, 0] = lower CI bound, conf[1, 1] = upper CI bound.
    """
    times, probs, conf = cumulative_incidence_competing_risks(
        events.astype(np.int32), durations, conf_type='log-log')
    t     = np.concatenate([[0.0], times])
    cif   = np.concatenate([[0.0], probs[1]])
    ci_lo = np.concatenate([[0.0], conf[1, 0]])
    ci_hi = np.concatenate([[0.0], conf[1, 1]])
    return t, cif, ci_lo, ci_hi


suffix    = '_noAPOE4_finegray'
apoe4     = False
recalc    = False
in_path   = 'dataset_all_percentile_BAI_noAPOE4.csv'

df = pd.read_csv(f'cutoff_results_lococv{suffix}.csv')

cutoff    = 90
out_plot  = f'figure1{suffix}_cut{cutoff}.png'
out_plot_pdf = f'figure1{suffix}_cut{cutoff}.pdf'
out_csv   = f'cic_dementia{suffix}_cut{cutoff}.csv'
apoe4_lbl = 'with APOE4' if apoe4 else 'no APOE4'

COLORS = {'high': '#d62728', 'low': '#1f77b4'}

df.loc[df.held_cohort=='FHS', 'held_cohort'] = 'FHS-OS'

cohort2color = {
    'MESA': '#CC79A7',
    'ARIC': '#009E73',
    'FHS-OS': '#0072B2',
    'MrOS': 'gray',
    'SOF': '#E69F00'}

plt.close()
fig = plt.figure(figsize=(11,6))
gs = fig.add_gridspec(2,2,height_ratios=[1,1], width_ratios=[5,6])
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0], sharex=ax1)

for cohort, color in cohort2color.items():
    mask = df.held_cohort==cohort
    cutoff_vals = df.cutoff[mask].values
    hr = df.HR[mask].values
    log_pval = -np.log10(df.pvalue[mask].values)
    best_id = np.argmax(log_pval)

    ax1.plot(cutoff_vals, hr, c=color, label=cohort)
    ax1.scatter(cutoff_vals[[best_id]], hr[[best_id]], c=color, marker='*', s=80)
    ax2.plot(cutoff_vals, log_pval, c=color, label=cohort)
    ax2.scatter(cutoff_vals[[best_id]], log_pval[[best_id]], c=color, marker='*', s=80)
ax1.yaxis.grid(True)
ax1.set_ylabel('Hazard ratio')
ax1.set_xlim(49, 96)
ax1.text(-0.14, 1.02, 'A', ha='right', va='top', transform=ax1.transAxes, fontweight='bold')
ax1.set_xlabel('')
plt.setp(ax1.get_xticklabels(), visible=False)
sns.despine(ax=ax1)

ax1.legend(title='CV Fold:', ncols=2, loc='upper left', alignment='left', framealpha=0.5)
ax2.text(49.2, -np.log10(0.001)+0.003, 'p = 0.001', ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
ax2.axhline(-np.log10(0.001), c='k', ls='--')
ax2.text(49.2, -np.log10(0.01)+0.003, 'p = 0.01', ha='left', va='bottom',
          bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
ax2.axhline(-np.log10(0.01), c='k', ls='--')
ax2.yaxis.grid(True)
ax2.set_ylabel('-log10(p)')
ax2.set_xlabel('BAI percentile cutoff (%)')
ax2.text(-0.14, 1.0, 'B', ha='right', va='top', transform=ax2.transAxes, fontweight='bold')
sns.despine(ax=ax2)



if os.path.exists(out_csv) and not recalc:
    print(f'Loading CIF values from {out_csv}')
    cif_df = pd.read_csv(out_csv)
else:
    print(f'Reading {in_path}')
    df = pd.read_csv(in_path)
    print(f'  rows: {len(df)},  events: {df.event.value_counts().to_dict()}')

    df['event_int'] = df['event'].map({'censor': 0, 'dementia': 1, 'death': 2})
    high_mask = df['BAIPercentile'] >= cutoff

    groups = {
        'high': (df[high_mask],  f'BAI ≥ {cutoff}% percentile'),
        'low':  (df[~high_mask], f'BAI < {cutoff}% percentile'),
    }

    records = []
    for key, (sub, legend_label) in groups.items():
        print(f'  Fitting CIF for {legend_label} (n={len(sub):,}) ...')
        t, cif, ci_lo, ci_hi = fit_cif(
            sub['time2event'].values, sub['event_int'].values)
        for ti, ci, lo, hi in zip(t, cif, ci_lo, ci_hi):
            records.append({
                'group':  legend_label,
                'time':   round(float(ti), 6),
                'CIF':    round(float(ci), 8),
                'CI_lo':  round(float(lo), 8),
                'CI_hi':  round(float(hi), 8),
            })

    cif_df = pd.DataFrame(records)
    cif_df.to_csv(out_csv, index=False)
    print(f'CIF values saved to: {out_csv}')

# ---- Plot ---------------------------------------------------------------
ax3 = fig.add_subplot(gs[:,1])

for key, color in COLORS.items():
    grp_label = (f'BAI ≥ {cutoff}% percentile' if key == 'high'
                    else f'BAI < {cutoff}% percentile')
    sub = cif_df[cif_df['group'] == grp_label]
    if sub.empty:
        continue
    ax3.step(sub['time'], sub['CIF'] * 100,
            where='post', color=color, linewidth=2, label=grp_label)
    ax3.fill_between(sub['time'], sub['CI_lo'] * 100, sub['CI_hi'] * 100,
                    step='post', alpha=0.15, color=color)

ax3.set_xlabel('Years since sleep study')
ax3.set_ylabel('Cumulative dementia incidence (%)')
#ax3.set_title(
#    f'Cumulative Incidence of Dementia (Aalen-Johansen, {apoe4_lbl})\n'
#    f'{cutoff}% percentile cutoff, 95% CI')
ax3.set_xlim(0,10)
ax3.set_ylim(0,16)
ax3.legend(loc='upper left', framealpha=0.5, ncols=1)
ax3.grid(True)#, alpha=0.3, linestyle='--')
ax3.text(-0.1, 1.01, 'C', ha='right', va='top', transform=ax3.transAxes, fontweight='bold')
sns.despine(ax=ax3)

plt.tight_layout()
plt.savefig(out_plot, dpi=300, bbox_inches='tight')
print(f'Plot saved to:       {out_plot}')
plt.savefig(out_plot_pdf, bbox_inches='tight')
print(f'Plot saved to:       {out_plot_pdf}')
plt.close()

print(f'\n{"Group":<35} {"5yr":>12} {"10yr":>12}')# {"15yr":>12}')
print('-' * 73)
for grp in cif_df['group'].unique():
    sub = cif_df[cif_df['group'] == grp]
    row = [f'{grp:<35}']
    for yr in [5, 10]:#, 15]:
        # Last step value at or before the requested year
        before = sub[sub['time'] <= yr]
        if before.empty:
            row.append(f'{"N/A":>12}')
        else:
            r = before.iloc[-1]
            row.append(f'{r.CIF*100:.1f}% [{r.CI_lo*100:.1f}%, {r.CI_hi*100:.1f}%]')
    print('  '.join(row))

