# Descriptive summary statistics for cohorts
# Author: Aashna Shah
# Date: 2023-06-09

import pandas as pd
import numpy as np
import sys

sys.path.append("../")
from data.dicts.col_dict import int_race_dict

def format_mean_std(series):
    vals = series.dropna()
    return f"{vals.mean():.1f} ({vals.std():.1f})"

def format_percent(series, value=1):
    pct = 100 * np.mean(series == value)
    return f"{pct:.1f}%"

def clean_binaries(df, binaries):
    for col in binaries:
        df[col] = df[col].map({'yes': 1, 'no': 0, True: 1, False: 0}).fillna(0).astype(int)
    return df

def summary_by_group(
    df,
    group_col,
    numeric_cols=None,
    binary_cols=None,
    group_order=None
):
    out = {}
    grouped = df.groupby(group_col)
    out['Sample Size'] = grouped.size()

    for label, col in numeric_cols.items():
        if col in df.columns:
            out[label] = grouped[col].apply(format_mean_std)

    for label, col in binary_cols.items():
        if col in df.columns:
            out[label] = grouped[col].apply(format_percent)
            
    summary = pd.DataFrame(out)
    summary = summary.T

    # Reorder race columns if group_order provided and group_col is 'race'
    if group_col == 'race' and group_order is not None:
        columns_to_keep = [g for g in group_order if g in summary.columns]
        extra_columns = [c for c in summary.columns if c not in columns_to_keep]
        summary = summary[columns_to_keep + extra_columns]

    return summary

NUMERIC_MAP = {
    'Age (SD) - years': 'age',
    "FEV1 (SD) - mL": 'fev1',
    "FVC (SD) - mL": 'fvc',
    "FEV1/FVC (SD)": 'fev1_fvc',
    "Height (SD) - cm": 'height',
    "Sitting Height (SD) - cm": 'sit_height',
    "Waist Circ. (SD) - cm": 'waist_circ',
}

BINARY_MAP = {
    "Female (%)": 'sex',             
    "Immigrant (%)": 'immigrant',
    "Education (%)": 'hs_grad',
    "Smoke Exposure (%)": 'smoke_exposure',
    "Income:Poverty (%)": 'inc_pov_bin',
}

ALL_COHORTS = ['ukb', 'nh', 'nh3', 'nh4']
RACE_ORDER = ['White', 'Black', 'Asian', 'Other', 'Hispanic']

def process_cohort(cohort):
    df = pd.read_csv(f"../data/processed/{cohort}/{cohort}_ref.csv")
    cols = ['sex', 'race', 'age', 'fev1', 'fvc', 'fev1_fvc', 'height', 'sit_height', 'waist_circ', 'immigrant', 'hs_grad', 'inc_pov_bin', 'smoke_exposure']
    df = df[cols]
    df = clean_binaries(df, BINARY_MAP.values())   
    return df

def main():
    for cohort in ALL_COHORTS:
        print(f"\n========== {cohort.upper()} ==========")
        df = process_cohort(cohort)
        # turn into int if it is  
        if cohort == 'ukb':
            df['race'] = pd.to_numeric(df['race'], errors='coerce').astype('Int64')
        df['race'] = df['race'].astype(str).map(int_race_dict).str.title()
        summary = summary_by_group(
            df,
            group_col='race',
            numeric_cols=NUMERIC_MAP,
            binary_cols=BINARY_MAP,
            group_order=RACE_ORDER,
        )
        print(summary.to_string())
        summary.to_csv(f"../results/stats/{cohort}_stats.csv", index=False)

if __name__ == "__main__":
    main()
