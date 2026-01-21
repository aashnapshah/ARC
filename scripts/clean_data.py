import warnings
warnings.filterwarnings('ignore')

import os
import sys
import json
import re
import pandas as pd
import numpy as np

sys.path.append("../")
from data.dicts.col_dict import *
from sklearn.preprocessing import StandardScaler

def get_reference(cohort):
    date = '2023-07-18' if cohort.startswith('nh') else '2023-12-04'
    path = f'../data/raw/{cohort}/{cohort}_{date}.csv'
    sep = '\t' if cohort.startswith('nh4') else ','
    df = pd.read_csv(path, sep=sep)
    df = df.replace('South Asian', 'Asian')

    def print_patient_info(df, label, verbose=False):
        race_values = df['RIDRETH'].value_counts().sort_index().to_dict()
        print(f"{label}: {df[id_col].nunique()} unique patients")
        if verbose:
            print("Patients by Race:", ", ".join([f"{race}: {count}" for race, count in race_values.items()]))
        print()

    id_col = 'SEQN' if cohort.startswith('nh') else 'EID'
    
    df = df.dropna(subset=['FEV1','FVC'])
    df['FEV1/FVC'] = df['FEV1']/df['FVC']

    df = df.dropna(subset=['RIDAGEYR'])
    df = df[(df['RIDAGEYR']>=20) & (df['RIDAGEYR']<=80)]

    admin_cols = list(set(admin_dict.keys()).intersection(df.columns))
    anthro_cols = list(set(anthropometrics_dict.keys()).intersection(df.columns))
    demo_cols = list(set(demographics_dict.keys()).intersection(df.columns))
    socio_cols = list(set(sociodemographics_dict.keys()).intersection(df.columns))
    bmx_cols = set(anthropometrics_dict.keys()).intersection(df.columns)
    smoke_cols = list(set(smoking_dict.keys()).intersection(df.columns))
    condition_cols = list(set(condition_dict.keys()).intersection(df.columns))
    symp_cols = list(set(symptoms_dict.keys()).intersection(df.columns))
    gen_cols = list(set(genetics_dict.keys()).intersection(df.columns))
    exp_cols = list(set(exposures_dict.keys()).intersection(df.columns))
    keep_cols = list(set(data_dict.keys()).intersection(df.columns))
    # remove 
    remove_cols = ['DMD_YEARS_IN_UK', 'born_subregion', 'poverty_ratio']
    df = df.dropna(subset=demo_cols)
    print_patient_info(df, "After dropping demo cols:", verbose=True)

    df = df.dropna(subset=anthro_cols)
    df = df[~df[smoke_cols].any(axis=1)]
    print_patient_info(df, "Smokers:", verbose=True)

    df = df[~df[condition_cols].any(axis=1)]
    print_patient_info(df, "No conditions:", verbose=True)

    df = df[df['FEV1/FVC'] >= 0.7]
    print_patient_info(df, "FEV1/FVC >= 0.7:", verbose=True)
    
    #df_gen = df_excl_symp.dropna(subset=gen_cols)
    #print_patient_info(df_gen, "Genetics")
    df = df[~df[symp_cols].any(axis=1)]
    print_patient_info(df, "No symptoms:", verbose=True)
    
    print(df[keep_cols].isna().sum().sort_values(ascending=False))
    df = df.dropna(subset=socio_cols)
    print_patient_info(df, "After dropping socio cols:", verbose=True)
    
    
    print(df[keep_cols].isna().sum().sort_values(ascending=False))
    # remove f.eid from keep_cols
    keep_cols = [col for col in keep_cols if col != 'f.eid']
    df = df[keep_cols].dropna().set_index(id_col).reset_index()
    print_patient_info(df, "Reference:", verbose=True)

    return df


def scale(df, numeric_cols=None):
    scale_cols = [col for col in numeric_cols if col in df.columns]
    scaler = StandardScaler()
    scaler.fit(df[scale_cols])
    df[scale_cols] = scaler.transform(df[scale_cols])
    return df

def combine_smoke_exposure(df):
    cols = df.columns[df.columns.str.contains('smoke', case=False)]
    df['smoke_exposure'] = df[cols].any(axis=1).astype(int)
    return df

def impute_sh(df):
    df.rename(columns={'sit_height': 'sit_height_old'}, inplace=True)
    imp_sh_path = '../results/imputed_sh/tables/imputed_results.csv'
    imp_sh = pd.read_csv(imp_sh_path).query('model == "_xgboost" and race_adj == False')
    imp_sh = imp_sh.applymap(lambda x: x.lower() if isinstance(x, str) else x)
    df = df.merge(imp_sh[['seqn', 'imputed_value']], on='seqn', how='left')
    df.rename(columns={'imputed_value': 'sit_height'}, inplace=True)
    print(df.head())
    return df

def unit_convert(df, targets=['fev1', 'fvc', 'fev1_fvc']):
    for col in targets:
        df[col] = df[col] * 1000
    return df

def get_cols(df):
    categorical_cols = []
    numeric_cols = []
    other_cols = []
    
    name_exclusions = ['eid', 'seqn', 'idx', 'fev1', 'fvc', 'gpc', 'age', 'strata', 'psu', 'WTMEC6YR']
    excluded_by_name = {col for col in df.columns if any(ex in col.lower() for ex in name_exclusions)}

    for col in df.columns:
        col_series = df[col]
        col_lower = col.lower()

        if col_lower in name_exclusions:
            continue

        if col_lower == 'race':
            categorical_cols.append(col)
            continue

        if col in excluded_by_name:
            other_cols.append(col)
            continue

        try:
            if pd.api.types.is_bool_dtype(col_series) or col_series.nunique(dropna=True) <= 2:
                other_cols.append(col)
                continue

            if pd.api.types.is_object_dtype(col_series) or pd.api.types.is_categorical_dtype(col_series) or pd.api.types.is_string_dtype(col_series):
                categorical_cols.append(col)
                continue

            if pd.api.types.is_numeric_dtype(col_series):
                if pd.api.types.is_integer_dtype(col_series):
                    unique_vals = col_series.nunique(dropna=True)
                    if col_lower == 'hh_size':
                        numeric_cols.append(col)
                    elif  unique_vals <= max(10, int(0.02 * max(len(col_series), 1))):
                        categorical_cols.append(col)
                    else:
                        numeric_cols.append(col)
                else:
                    numeric_cols.append(col)
                continue
        except:
            print(f"Error processing column {col}")

            continue
        # Fallback: unknown dtypes
        other_cols.append(col)   

    return categorical_cols, numeric_cols, other_cols

def bin_income(df, cohort):
    if 'nh' in cohort:
        df['inc_pov_bin'] = (df['poverty_ratio'] > 1).astype(int)
    else: 
        df['inc_pov_bin'] = (df['income_category'] > 17760).astype(int)
    return df

def encode(df, categorical_cols=None):
    for col in list(categorical_cols):
        if col not in df.columns:
            continue   
        if col == 'race':
            df['race'] = pd.Categorical(
            df['race'], 
            categories=[0, 1, 2, 3, 4], 
            ordered=True
            )
            df = pd.get_dummies(df, columns=[col], drop_first=True)
        else:
            df = pd.get_dummies(df, columns=[col], drop_first=False)
    return df

def combine_nh(cohort):
    path = f'../data/processed/{cohort}/{cohort}_ref.csv'
    nh3_ref = pd.read_csv(path.replace('nh', 'nh3'))
    nh4_ref = pd.read_csv(path.replace('nh', 'nh4'))
    nhanes_cols = list(set(nh3_ref.columns).intersection(nh4_ref.columns))
    df = pd.concat([nh3_ref, nh4_ref])[nhanes_cols]
    return df.rename(columns=names_dict) 

def get_sex(df):
    intersection = list(set(df.columns).intersection(['RIAGENDR', 'IS_FEMALE', 'sex']))
    if 'RIAGENDR' in df.columns and 'IS_FEMALE' in df.columns:
        df = df.drop(columns=['RIAGENDR'])
    if 'sex' in df.columns:
        df = df.drop(columns=['sex'])
    return df

def immigrant(df, cohort):
    born_col = 'born_uk' if 'ukb' in cohort else 'born_us'
    if born_col == 'born_uk':
        df['immigrant'] = df['born_uk'] != 'uk'
    else:
        df['immigrant'] = df['born_us'] == 0
    return df

def clean_df(cohort, savePath=None):
    if cohort != 'nh':
        df = get_reference(cohort)
        df = get_sex(df)
        df = df.rename(columns=names_dict) 
        df = bin_income(df, cohort)
        df = df.replace(' ', '_', regex=True)
        df = df[df['age'] >= 40]
        df = df.applymap(lambda x: x.lower() if isinstance(x, str) else x)
        df = df.replace('sub-', 'sub_', regex=True)
        df = impute_sh(df) if 'nh4' in cohort else df
        df = immigrant(df, cohort)
        df['race'] = df['race'].map(race_dict)
        df = combine_smoke_exposure(df)
        df = unit_convert(df, targets=['fev1', 'fvc', 'fev1_fvc'])
    else:
        df = combine_nh(cohort)
        
    categorical_cols, numeric_cols, other_cols = get_cols(df)
    print()
    print(f"Categorical: {' '.join(categorical_cols)}")
    print(f"Numeric: {' '.join(numeric_cols)}")
    print(f"Other: {' '.join(other_cols)}")
    print()
    
    df_encoded = encode(df, categorical_cols=categorical_cols)
    df_scaled = scale(df_encoded, numeric_cols=numeric_cols)
    df_scaled_encoded = scale(df_encoded, numeric_cols=numeric_cols)
    
    if savePath:
        df.to_csv(savePath + '_ref.csv', index=False)
        df_encoded.to_csv(savePath + '_ref_d.csv', index=False)
        df_scaled.to_csv(savePath + '_ref_std.csv', index=False)
        df_scaled_encoded.to_csv(savePath + '_ref_std_d.csv', index=False)
        print(f"Saved to {savePath}")
        
    # df.to_csv(path.replace('.csv', '_imp.csv'), index=False)
    # df_encoded = encode(df, path.replace('.csv', '_encoded.csv'), categorical_cols=categorical_cols)
    # df_scaled = scale(df, path.replace('.csv', '_scaled.csv'), numeric_cols=numeric_cols)
    # df_scaled_encoded = scale(df_encoded, path.replace('.csv', '_scaled_encoded.csv'), numeric_cols=numeric_cols)
    return df

def main():
    cohorts = ['ukb', 'nh3', 'nh4', 'nh'] 
    for cohort in cohorts:
        print('-'*100)
        print(f"Processing {cohort}...")
        out_path = f'../data/processed/{cohort}/{cohort}'
        
        df = clean_df(cohort, savePath=out_path)
        df['inc_pov_bin'] = df['inc_pov_bin'].astype(int)
        print('-'*100)
        
    # data_dir = "../data/processed"
    # cohorts = ['nh3', 'nh4', 'nh', 'ukb']
    # # for cohort in cohorts:
    #     # path = os.path.join(data_dir, f"{cohort}_ref.csv")
    #     # df = clean_df(path)
    #   #  df.to_csv(path.replace('.csv', '_imp.csv'), index=False)

if __name__ == "__main__":
    main()