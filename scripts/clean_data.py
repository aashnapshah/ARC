import pandas as pd
import numpy as np

import os
import sys
from sklearn.preprocessing import StandardScaler
import json
import re

sys.path.append("../")
from data.dicts.col_dict import col_names_dict

def scale(df, path=None, numeric_cols=None):
    """
    Scale provided numeric columns (if given). If not provided, infer numeric columns.
    """
    scale_cols = [col for col in numeric_cols if col in df.columns]
    scaler = StandardScaler()
    scaler.fit(df[scale_cols])
    df[scale_cols] = scaler.transform(df[scale_cols])
    if path:
        df.to_csv(path, index=False)
    return df

def impute_sh(df):
    df.rename(columns={'sit_height': 'sit_height_old'}, inplace=True)
    imp_sh_path = '/n/data1/hms/dbmi/manrai/aashna/ARC/results/imputed_sh/imputed_results_ref.csv'
    imp_sh = pd.read_csv(imp_sh_path).set_index('IDX', drop=True).query('model == "_xgboost" and race_adj == False')
    df = df.reset_index(drop=True).merge(imp_sh, right_index=True, left_index=True, how='left')
    df.rename(columns={'imputed_value': 'sit_height'}, inplace=True)
    return df

def read_and_prepare(df, targets=['fev1', 'fvc', 'fev1_fvc']):
    for col in targets:
        df[col] = df[col] * 1000
    return df

def get_cols(df):
    """
    Determine which columns are categorical (to encode), numeric (to scale),
    and other (ids, targets, booleans, already-binary, excluded by name).
    """
    categorical_cols = []
    numeric_cols = []
    other_cols = []

    # Columns that should never be scaled (ids, indexes, targets, etc.)
    name_exclusions = ['eid', 'seqn', 'idx', 'fev1', 'fvc', 'gpc', 'age']
    excluded_by_name = {col for col in df.columns if any(ex in col.lower() for ex in name_exclusions)}

    for col in df.columns:
        col_series = df[col]
        col_lower = col.lower()

        if col_lower == 'race':
            categorical_cols.append(col)
            continue

        # Exclude certain columns by name from scaling/encoding
        if col in excluded_by_name:
            other_cols.append(col)
            continue

        # Booleans and binary flags (already encoded)
        if pd.api.types.is_bool_dtype(col_series) or col_series.nunique(dropna=True) <= 2:
            other_cols.append(col)
            continue

        # Object/string/categorical -> categorical
        if pd.api.types.is_object_dtype(col_series) or pd.api.types.is_categorical_dtype(col_series) or pd.api.types.is_string_dtype(col_series):
            categorical_cols.append(col)
            continue

        # Numeric
        if pd.api.types.is_numeric_dtype(col_series):
            # Integer with low cardinality (heuristic) -> categorical
            if pd.api.types.is_integer_dtype(col_series):
                unique_vals = col_series.nunique(dropna=True)
                if unique_vals <= max(10, int(0.02 * max(len(col_series), 1))):
                    categorical_cols.append(col)
                else:
                    numeric_cols.append(col)
            else:
                numeric_cols.append(col)
            continue

        # Fallback: unknown dtypes
        other_cols.append(col)

    return categorical_cols, numeric_cols, other_cols

def encode(df, path=None, categorical_cols=None):
    """
    One-hot encode provided categorical columns (if given). If not provided,
    fallback to original automatic detection.
    """
    for col in list(categorical_cols):
        if col not in df.columns:
            continue
        if col == 'race':
            df[col] = df[col].map({'white': 0, 'black': 1, 'asian': 2, 'hispanic': 3, 'other': 4})
            df = pd.get_dummies(df, columns=[col], drop_first=False)
        else:
            df = pd.get_dummies(df, columns=[col], drop_first=False)

    if path:
        df.to_csv(path, index=False)
    return df

def clean_df(path):
    df = pd.read_csv(path)
    df = df.rename(columns=col_names_dict) 
    
    df = df.replace(' ', '_', regex=True)
    df = df.applymap(lambda x: x.lower() if isinstance(x, str) else x)
    df = df.replace('sub-', 'sub_', regex=True)
    
    born_col = 'born_uk' if 'ukb' in path else 'born_us'
    df['immigrant'] = df[born_col] == 0
    
    if 'nh4' in path:
        df = impute_sh(df)

    # Determine columns for encoding/scaling and print them
    categorical_cols, numeric_cols, other_cols = get_cols(df)
    print(f"Categorical columns to encode ({len(categorical_cols)}): {categorical_cols}")
    print(f"Numeric columns to scale ({len(numeric_cols)}): {numeric_cols}")
    print(f"Other columns ({len(other_cols)}): {other_cols}")
    
    df = read_and_prepare(df, targets=['fev1', 'fvc', 'fev1_fvc'])
    df.to_csv(path.replace('.csv', '_imp.csv'), index=False)
    df_encoded = encode(df, path.replace('.csv', '_encoded.csv'), categorical_cols=categorical_cols)
    df_scaled = scale(df, path.replace('.csv', '_scaled.csv'), numeric_cols=numeric_cols)
    df_scaled_encoded = scale(df_encoded, path.replace('.csv', '_scaled_encoded.csv'), numeric_cols=numeric_cols)
    print(path.replace('.csv', '_scaled_encoded.csv'))
    return df

def main():
    data_dir = "../data/processed"
    nh3_ref_path = os.path.join(data_dir, "nh3_ref.csv")
    nh4_ref_path = os.path.join(data_dir, "nh4_ref.csv")
    ukb_ref_path = os.path.join(data_dir, "ukb_ref.csv")

    nh3_ref = clean_df(nh3_ref_path)
    nh4_ref = clean_df(nh4_ref_path)
    ukb_ref = clean_df(ukb_ref_path)
    

if __name__ == "__main__":
    main()