import warnings
warnings.filterwarnings('ignore')

import os
import sys
import argparse
from itertools import product
import numpy as np
import pandas as pd
from patsy import dmatrices
import statsmodels.api as sm
from statsmodels.stats.weightstats import DescrStatsW

sys.path.append('../')
from data.dicts.col_dict import *

def get_params(model, col, outcome):
    cols = ['index', 'coef', 'std_err', 't', 'p', 'ci_lower', 'ci_upper']
    res = model.summary2().tables[1].reset_index() 
    res.columns = cols
    pred = [pred for pred in predictors if col.startswith(pred)][0]
    value = col.replace(pred + "_", "") if '_' in col else None
    res['target'] = outcome
    res['group'] = groups_dict[pred]
    res['cov'] = pred
    res['val'] = value
    return res[['target', 'group', 'cov', 'val'] + cols]

def run_ewas(df, targets):
    results = []
    relevant_predictor_names = set()
    for d in (anthropometrics_dict, demographics_dict, sociodemographics_dict, exposures_dict):
        relevant_predictor_names.update([v[0] for v in d.values()])

    cols = []
    for c in df.columns:
        for base in relevant_predictor_names:
            if c == base or c.startswith(base + "_"):
                cols.append(c)
                break

    pairs = product(targets, cols)
    for outcome, col in pairs:
        try:
            formula = f"{outcome} ~ age + sex + {col}"
            y, X = dmatrices(formula, data=df, return_type='dataframe')
            model = sm.OLS(y, X).fit()
            res = get_params(model, col, outcome)
            results.append(res)
        except Exception as e:
            continue
    return pd.concat(results, ignore_index=True)

def process_cohort(df, cohort):
    df = df[df.apply(lambda row: str(row['cov']) in str(row['index']), axis=1)]
    df['title'] = df['cov'].map(titles_dict)
    df = df[(df['group'].notna()) & (df['group'] != 'NA')]
    for i, row in df.iterrows():
        if "race" in row['cov'] and row['val'] is not np.nan:
            val = row['val']
            df.at[i, 'title'] = int_race_dict[val].title()
            df.at[i, 'group'] = 'Race'
    round_cols_2dp = ['coef', 'std_err', 'ci_lower', 'ci_upper']
    for col in round_cols_2dp:
        df[col] = df[col].round(2)
    df['p'] = df['p'].round(3)
    df['cohort'] = cohort.upper()
    df['group'] = df['group'].str.title()
    df['target'] = df['target'].str.upper()
    df.sort_values(by=['cohort', 'target', 'group', 'title'], inplace=True)
    cols = ['cohort', 'target', 'group', 'title', 'coef', 'std_err', 'p']
    return df[cols]

def main(args):
    for cohort in args.cohorts:
        df = pd.read_csv(f"{args.data_dir}/{cohort}/{cohort}_ref_std_d.csv")
        result = run_ewas(df, args.targets)
        result.to_csv(f"{args.output_dir}/tables/raw/{cohort}_ref_ewas.csv", index=False)
        result = process_cohort(result, cohort)
        print(result.head())
        result.to_csv(f"{args.output_dir}/tables/processed/{cohort}_ref_ewas.csv", index=False)
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ewas workflow.")
    parser.add_argument('--data_dir', default='../data/processed')
    parser.add_argument('--output_dir', default='../results/ewas/')
    parser.add_argument('--cohorts', nargs='+', default=['ukb', 'nh', 'nh3', 'nh4'])
    parser.add_argument('--targets', nargs='+', default=['fev1', 'fvc'])
    args = parser.parse_args()
    main(args)