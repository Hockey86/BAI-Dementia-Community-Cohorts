import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.proportion import proportions_ztest


def format_pvalue(p_value):
    """Format p-value according to specified rules"""
    if np.isnan(p_value):
        return 'N/A'
    elif p_value < 0.001:
        return 'p<0.001'
    elif p_value <= 0.01:
        return f'{p_value:.3f}'
    elif p_value <= 0.99:
        return f'{p_value:.2f}'
    else:
        return 'p>0.99'


def main():
    cohort_paths = {
        'MESA': '/data/haoqisun/BAI_dementia_community2/MESA/dataset_MESA.csv',
        'ARIC': '/data/haoqisun/BAI_dementia_community2/SHHS/ARIC/dataset_ARIC.csv',
        'FHS': '/data/haoqisun/BAI_dementia_community2/SHHS/FHS/dataset_FHS.csv',
        'MrOS': '/data/haoqisun/BAI_dementia_community2/MrOS/dataset_MrOS.csv',
        'SOF': '/data/haoqisun/BAI_dementia_community2/SOF/dataset_SOF.csv'
    }
    
    # Ordered variables list with updated names and column mappings
    ordered_variables = [
        ('Age at sleep study, year', 'age', 'continuous'),
        ('Sex (Female)', 'sexM', 'categorical'),
        ('BMI, kg/m2', 'BMI', 'continuous'),
        ('College degree', 'educcollege', 'categorical'),
        ('APOE4 carrier', 'APOE4Count', 'categorical'),  # Will handle SOF separately
        ('AHI, /hour', 'AHI', 'continuous'),
        ('Sleep medication use', 'sleepmed', 'categorical'),
        ('Hypertension', 'hypertension', 'categorical'),
        ('Diabetes', 'diabetes', 'categorical'),
        ('Myocardial infarction', 'heartattack', 'categorical'),  # Updated for most cohorts
        ('Stroke', 'stroke', 'categorical'),
        ('Depression', 'depression', 'categorical')
    ]
    
    results = []
    
    for cohort_name, file_path in cohort_paths.items():
        if not os.path.exists(file_path):
            print(f"Warning: {file_path} not found, skipping {cohort_name}")
            continue
            
        try:
            df = pd.read_csv(file_path)
            print(f"Processing {cohort_name}: {len(df)} participants")
            
            # Updated BAI groups: <= -3 vs >= 3
            df_young = df[df['BAI'] <= -3].copy()
            df_old = df[df['BAI'] >= 3].copy()
            
            n_young = len(df_young)
            n_old = len(df_old)
            
            # Add cohort header
            results.append({
                'Variable': cohort_name,
                'Young BAI (<=-3y)': '',
                'Old BAI (>=+3y)': '',
                'P-value': ''
            })
            
            # Add N row
            results.append({
                'Variable': 'N',
                'Young BAI (<=-3y)': str(n_young),
                'Old BAI (>=+3y)': str(n_old),
                'P-value': ''
            })
            
            # Process variables in the specified order
            for var_info in ordered_variables:
                var_name, col_name, var_type = var_info
                
                # Handle special column mappings
                actual_col_name = col_name
                
                # Handle MESA-specific myocardial infarction column
                if var_name == 'Myocardial infarction' and cohort_name == 'MESA':
                    actual_col_name = 'mi'
                
                # Handle SOF-specific APOE4 column
                if var_name == 'APOE4 carrier' and cohort_name == 'SOF':
                    actual_col_name = 'APOE4'
                
                if var_type == 'continuous':
                    # Handle continuous variables
                    if actual_col_name in df.columns:
                        young_data = df_young[actual_col_name].dropna()
                        old_data = df_old[actual_col_name].dropna()
                        
                        if len(young_data) > 0 and len(old_data) > 0:
                            young_mean = young_data.mean()
                            young_std = young_data.std()
                            old_mean = old_data.mean()
                            old_std = old_data.std()
                            
                            try:
                                t_stat, p_value = ttest_ind(young_data, old_data)
                                p_str = format_pvalue(p_value)
                            except:
                                p_str = 'N/A'
                            
                            results.append({
                                'Variable': var_name,
                                'Young BAI (<=-3y)': f'{young_mean:.1f} ({young_std:.1f})',
                                'Old BAI (>=+3y)': f'{old_mean:.1f} ({old_std:.1f})',
                                'P-value': p_str
                            })
                        else:
                            results.append({
                                'Variable': var_name,
                                'Young BAI (<=-3y)': 'N/A',
                                'Old BAI (>=+3y)': 'N/A',
                                'P-value': 'N/A'
                            })
                    else:
                        print(f"Warning: {actual_col_name} column not found in {cohort_name}")
                        results.append({
                            'Variable': var_name,
                            'Young BAI (<=-3y)': 'N/A',
                            'Old BAI (>=+3y)': 'N/A',
                            'P-value': 'N/A'
                        })
                
                elif var_type == 'categorical':
                    # Handle categorical variables
                    if actual_col_name not in df.columns:
                        print(f"Warning: {actual_col_name} column not found in {cohort_name}")
                        results.append({
                            'Variable': var_name,
                            'Young BAI (<=-3y)': 'N/A',
                            'Old BAI (>=+3y)': 'N/A',
                            'P-value': 'N/A'
                        })
                        continue
                    
                    # Handle APOE4 carrier
                    if var_name == 'APOE4 carrier':
                        young_data = df_young[actual_col_name].dropna()
                        old_data = df_old[actual_col_name].dropna()
                        
                        if cohort_name == 'SOF':
                            # SOF APOE4 is already binary
                            young_pos = (young_data == 1).sum()
                            old_pos = (old_data == 1).sum()
                        else:
                            # Other cohorts: convert count to binary
                            young_binary = (young_data > 0).astype(int)
                            old_binary = (old_data > 0).astype(int)
                            young_pos = young_binary.sum()
                            old_pos = old_binary.sum()
                        
                        young_total = len(young_data)
                        old_total = len(old_data)
                    else:
                        # Regular binary variables
                        young_data = df_young[actual_col_name].dropna()
                        old_data = df_old[actual_col_name].dropna()
                        
                        # For Sex (Female), count sexM=0 (females), for others count =1
                        if var_name == 'Sex (Female)':
                            young_pos = (young_data == 0).sum()
                            old_pos = (old_data == 0).sum()
                        else:
                            young_pos = (young_data == 1).sum()
                            old_pos = (old_data == 1).sum()
                        young_total = len(young_data)
                        old_total = len(old_data)
                    
                    if young_total > 0 and old_total > 0:
                        young_pct = (young_pos / young_total) * 100
                        old_pct = (old_pos / old_total) * 100
                        
                        # Calculate p-value using statsmodels proportions z-test
                        try:
                            if young_pos + old_pos > 0 and young_total + old_total > young_pos + old_pos:
                                z_stat, p_value = proportions_ztest([young_pos, old_pos], [young_total, old_total])
                                p_str = format_pvalue(p_value)
                            else:
                                p_str = 'N/A'
                        except:
                            p_str = 'N/A'
                        
                        results.append({
                            'Variable': var_name,
                            'Young BAI (<=-3y)': f'{young_pos} ({young_pct:.1f}%)',
                            'Old BAI (>=+3y)': f'{old_pos} ({old_pct:.1f}%)',
                            'P-value': p_str
                        })
                    else:
                        results.append({
                            'Variable': var_name,
                            'Young BAI (<=-3y)': 'N/A',
                            'Old BAI (>=+3y)': 'N/A',
                            'P-value': 'N/A'
                        })
                    
        except Exception as e:
            print(f"Error processing {cohort_name}: {e}")
            continue
    
    df_results = pd.DataFrame(results)
    print(df_results)
    df_results.to_excel('BAI_descriptive_results.xlsx', index=False)
    


if __name__ == "__main__":
    main()

