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
df_test  = pd.read_csv(f'cutoff_results_lococv_test{suffix}.csv')
df_boot_hr = pd.read_csv(f'cutoff_bootstrap_HR{suffix}.csv')

cutoff    = 90
out_plot  = f'figure1{suffix}_cut{cutoff}.png'
out_plot_pdf = f'figure1{suffix}_cut{cutoff}.pdf'
out_csv   = f'cic_dementia{suffix}_cut{cutoff}.csv'
apoe4_lbl = 'with APOE4' if apoe4 else 'no APOE4'

COLORS = {'high': '#d62728', 'low': '#1f77b4'}

for d in [df, df_test]:
    d.loc[d.held_cohort=='FHS', 'held_cohort'] = 'FHS-OS'

cohort2color = {
    'MESA': '#CC79A7',
    'ARIC': '#009E73',
    'FHS-OS': '#0072B2',
    'MrOS': 'gray',
    'SOF': '#E69F00'}

plt.close()
fig = plt.figure(figsize=(14,6))
gs = fig.add_gridspec(2,3,height_ratios=[1,1], width_ratios=[4,4,6])


def plot_sweep(ax_hr, ax_p, lines, panel_hr, panel_p):
    """HR and -log10(p) vs cutoff. lines: list of (df, color, label, mark_selected);
    star = cutoff selected on the training cohorts."""
    for d, color, label, mark_selected in lines:
        cutoff_vals = d.cutoff.values
        hr = d.HR.values
        log_pval = -np.log10(d.pvalue.values)

        ax_hr.plot(cutoff_vals, hr, c=color, label=label)
        ax_p.plot(cutoff_vals, log_pval, c=color, label=label)
        if mark_selected:
            sel = d.selected.values.astype(bool)
            ax_hr.scatter(cutoff_vals[sel], hr[sel], color=color, marker='*', s=80)
            ax_p.scatter(cutoff_vals[sel], log_pval[sel], color=color, marker='*', s=80)
    ax_hr.yaxis.grid(True)
    ax_hr.set_ylabel('Hazard ratio')
    ax_hr.set_xlim(49, 96)
    ax_hr.set_xticks([50, 60, 70, 80, 90])
    ax_hr.text(-0.17, 1.02, panel_hr, ha='right', va='top', transform=ax_hr.transAxes, fontweight='bold')
    ax_hr.set_xlabel('')
    plt.setp(ax_hr.get_xticklabels(), visible=False)
    sns.despine(ax=ax_hr)

    for p in [0.001, 0.01]:
        ax_p.text(49.2, -np.log10(p)+0.003, f'p = {p:g}', ha='left', va='bottom',
                  bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
        ax_p.axhline(-np.log10(p), c='k', ls='--')
    ax_p.yaxis.grid(True)
    ax_p.set_ylabel('-log10(p)')
    ax_p.set_xlabel('BAI percentile threshold (%)')
    ax_p.text(-0.17, 1.0, panel_p, ha='right', va='top', transform=ax_p.transAxes, fontweight='bold')
    sns.despine(ax=ax_p)


# A, B: training cohorts of each fold (used to select the cutoff)
ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
plot_sweep(ax1, ax2,
           [(df[df.held_cohort==c], color, c, True) for c, color in cohort2color.items()],
           'A', 'B')
ax1.set_ylim(0.95, 2.0)  # leave room for the legend, and keep HR = 1 off the bottom axis
leg1 = ax1.legend(title='LOCO-CV Training Folds:', ncols=2, loc='upper left', alignment='left',
                  frameon=False, columnspacing=0.8)
leg1.get_title().set_fontweight('bold')

# C, D: held-out cohorts of all folds pooled together (evaluation only)
ax3 = fig.add_subplot(gs[0,1])
ax4 = fig.add_subplot(gs[1,1], sharex=ax3)
plot_sweep(ax3, ax4, [(df_test[df_test.held_cohort=='Pooled'], 'k', None, False)], 'C', 'D')
ax3.fill_between(df_boot_hr.cutoff, df_boot_hr.HR_boot_lo, df_boot_hr.HR_boot_hi,
                 color='gray', alpha=0.25, lw=0)   # bootstrap 95% CI
txt3 = ax3.text(0.03, 0.97, 'Pooled LOCO-CV Testing Folds:', ha='left', va='top', transform=ax3.transAxes,
         fontsize=plt.rcParams['legend.title_fontsize'], fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.5, edgecolor='none'))
ax3.axhline(1, c='k', lw=1)

# same y-axis range for A & C and for B & D
for ax_a, ax_b in [(ax1, ax3), (ax2, ax4)]:
    lo = min(ax_a.get_ylim()[0], ax_b.get_ylim()[0])
    hi = max(ax_a.get_ylim()[1], ax_b.get_ylim()[1])
    ax_a.set_ylim(lo, hi)
    ax_b.set_ylim(lo, hi)
# HR ticks every 0.2, skipping one near the top that would collide with the panel label
hr_ticks = np.arange(1.0, ax1.get_ylim()[1] - 0.1, 0.2)
ax1.set_yticks(hr_ticks)
ax3.set_yticks(hr_ticks)


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
ax5 = fig.add_subplot(gs[:,2])

for key, color in COLORS.items():
    grp_label = (f'BAI ≥ {cutoff}% percentile' if key == 'high'
                    else f'BAI < {cutoff}% percentile')
    sub = cif_df[cif_df['group'] == grp_label]
    if sub.empty:
        continue
    ax5.step(sub['time'], sub['CIF'] * 100,
            where='post', color=color, linewidth=2, label=grp_label.replace('%','th'))
    ax5.fill_between(sub['time'], sub['CI_lo'] * 100, sub['CI_hi'] * 100,
                    step='post', alpha=0.15, color=color)

ax5.set_xlabel('Years since sleep study')
ax5.set_ylabel('Cumulative dementia incidence (%)')
ax5.set_xlim(0,10)
ax5.set_ylim(0,16)
ax5.legend(loc='upper left', framealpha=0.5, ncols=1)
ax5.grid(True)#, alpha=0.3, linestyle='--')
ax5.text(-0.15, 1.01, 'E', ha='right', va='top', transform=ax5.transAxes, fontweight='bold')
sns.despine(ax=ax5)

plt.tight_layout(w_pad=0.2)

# top-align the panel C text with the legend title in panel A (same row, so same axes-fraction height)
fig.canvas.draw()
title_top = leg1.get_title().get_window_extent().y1
txt3.set_y(ax1.transAxes.inverted().transform((0, title_top))[1])
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

