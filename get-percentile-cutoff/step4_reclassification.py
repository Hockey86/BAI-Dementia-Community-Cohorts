import numpy as np
import pandas as pd


df = pd.read_csv('dataset_all_percentile_BAI_noAPOE4.csv')
#BAI_median = df.BAI.median()
#print(BAI_median)
BAI_median = 0.
cutoff = 90

rows = []

groups = [
    (f'BAI>0 & Percentile>{cutoff}', (df.BAI>BAI_median) & (df.BAIPercentile>cutoff)),
    (f'BAI<=0 & Percentile>{cutoff}',  (df.BAI<=BAI_median)  & (df.BAIPercentile>cutoff)),
    (f'BAI>0 & Percentile<={cutoff}',  (df.BAI>BAI_median) & (df.BAIPercentile<=cutoff)),
    (f'BAI<=0 & Percentile<={cutoff}',   (df.BAI<=BAI_median)  & (df.BAIPercentile<=cutoff)),
    (f'BAI>0',   (df.BAI>BAI_median) ),
    (f'BAI<=0',   (df.BAI<=BAI_median) ),
    (f'Percentile>{cutoff}', (df.BAIPercentile>cutoff)),
    (f'Percentile<={cutoff}', (df.BAIPercentile<=cutoff)),
]

for label, mask in groups:
    n = mask.sum()
    pct = mask.mean() * 100
    dementia_pct = (df.event[mask] == "dementia").mean() * 100
    death_pct    = (df.event[mask] == "death").mean()    * 100
    print(f'{n} ({pct:.1f}%), dementia = {dementia_pct:.1f}%, death = {death_pct:.1f}%')
    rows.append({
        'group':        label,
        'n':            n,
        'pct_of_total': round(pct, 1),
        'dementia_pct': round(dementia_pct, 1),
        'death_pct':    round(death_pct, 1),
    })

output_file = f'reclassification_results_noAPOE4_cutoff{cutoff}.csv'
results = pd.DataFrame(rows)
results.to_csv(output_file, index=False)
print(f'Saved to {output_file}')