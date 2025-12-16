import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 12})
import seaborn as sns
sns.set_style('ticks')


feat_names = [
'BAI',

'delta_bandpower_mean_C_N3',
'delta_alpha_mean_C_N3',
'delta_theta_mean_C_N3',
'theta_bandpower_kurtosis_C_N3',
'kurtosis_N3_C',

'COUPL_OVERLAP_C',
'DENS_C',
'sigma_bandpower_kurtosis_C_N2',
'delta_bandpower_kurtosis_C_N2',
'theta_bandpower_kurtosis_C_N2',
'alpha_bandpower_kurtosis_C_N2',
'kurtosis_N2_C',

'alpha_bandpower_mean_C_N1',
]
feat_names2 = [
'BAI',

'Delta power in N3',
'Delta-alpha ratio in N3',
'Delta-theta ratio in N3',
'Theta kurtosis in N3',
'Waveform kurtosis in N3',

'Spindle-SO overlap% in N2',
'Spindle density in N2',
'Sigma kurtosis in N2',
'Delta kurtosis in N2',
'Theta kurtosis in N2',
'Alpha kurtosis in N2',
'Waveform kurtosis in N2',

'Alpha power in N1',
]

#"""
df = {'Name':[], 'HR':[], 'LB':[], 'UB':[], 'P':[]}
for fn in feat_names:
    with open(os.path.join('withoutAPOE', f'meta_analysis_result-survival-{fn}-Intermediate.txt'), 'r') as ff:
        for l in ff:
            if l.startswith('Random effects model'):
                #try:
                xx = l[len('Random effects model '):].split()
                df['Name'].append(fn)
                df['HR'].append( float(xx[0]) )
                df['LB'].append( float(xx[1][1:-1]) )
                df['UB'].append( float(xx[2][:-1]) )
                df['P'].append(  ''.join(xx[4:]) )
                if df['P'][-1]=='<0.0001':
                    df['P'][-1] = 0.
                else:
                    df['P'][-1] = float(df['P'][-1])
                #except Exception as ee:
                #    print(fn, l, str(ee))
                break

df = pd.DataFrame(data=df)
print(df)
df.to_csv('BAI_feature_coefs.csv', index=False)
#"""
#df = pd.read_csv('BAI_feature_coefs.csv')

plt.close()
fig, ax = plt.subplots(figsize=(10*0.7, 8*0.7))

# Define colors for different sleep stages
def get_color(feature_name):
    if 'N3' in feature_name:
        return  ((149-10)/255,(220-5)/255,(221+5)/255)
    elif 'N2' in feature_name:
        return (141/255,(188+20)/255,133/255)
    elif 'N1' in feature_name:
        return (252/1.1/255,233/1.1/255,142/1.1/255)
    elif feature_name == 'BAI':
        return (107/255,107/255,107/255)

# Create y-positions for features (reverse order for top-to-bottom display)
y_positions = range(len(df), 0, -1)

# Add vertical dashed line at HR = 1
ax.axvline(x=1, color='red', linestyle='--', alpha=0.7, linewidth=2)

# Plot horizontal lines with confidence intervals
for i, (idx, row) in enumerate(df.iterrows()):
    y_pos = y_positions[i]
    color = get_color(feat_names2[i])
    
    # Plot horizontal line from LB to UB
    ax.plot([row['LB'], row['UB']], [y_pos, y_pos], 
            color=color, linewidth=2, alpha=0.7)
    
    # Plot center point (HR)
    ax.plot(row['HR'], y_pos, 'o', color=color, markersize=8, 
            markeredgecolor='black', markeredgewidth=0.5)
    
    # Add significance markers
    if row['P'] < 0.05:
        if row['HR'] < 1:  # Significantly protective
            ax.text(row['LB'] - 0.02, y_pos-0.21, '*', fontsize=18, 
                   color='k', ha='right', va='center', weight='bold')
        else:  # Significantly harmful
            ax.text(row['UB'] + 0.02, y_pos-0.21, '*', fontsize=18, 
                   color='k', ha='left', va='center', weight='bold')

# Set y-axis labels with feature names
ax.set_yticks(y_positions)
ax.set_yticklabels([feat_names2[i] for i in range(len(feat_names2))])

# Set axis labels
ax.set_xlabel('Hazard Ratio (HR)', weight='bold')
#ax.set_ylabel('Features', fontsize=12, weight='bold')

# Set x-axis limits with some padding
x_min = df[['LB', 'UB', 'HR']].min().min() * 0.9
x_max = df[['LB', 'UB', 'HR']].max().max() * 1.1
ax.set_xlim(x_min, x_max)
ax.set_xticks(np.arange(0.7,1.8,0.1))

# Add grid
ax.grid(True, alpha=0.6, axis='x')

# Remove top and right spines
sns.despine()

plt.tight_layout()
#plt.savefig('BAI_features_forest_plot.png', dpi=300, bbox_inches='tight')
plt.savefig('BAI_features_forest_plot.pdf', bbox_inches='tight')
#plt.show()

