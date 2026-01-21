# Cleaned and Dynamic Refactor of race_adjusted.R Logic (Python)

import os
import pandas as pd
import numpy as np
import argparse
import sys
from itertools import product

sys.path.append('../')
from data.dicts.col_dict import *
from statsmodels.formula.api import ols

def mse(model):
    return np.mean(model.resid**2)

def rss(model):
    return np.sum(model.resid**2)

# Re-write get_predictors to use only vars in the dict for each category
def get_predictors(df, groups=['anthropometrics', 'sociodemographics', 'exposures']):
    predictors = []
    for group in groups:
        valid_cols = groups_list_dict[group]
        intersection = set(valid_cols) & set(df.columns)
        predictors.extend([col for col in intersection if not (col.endswith('adj') or col.endswith('old'))])
    return predictors

def get_params(model, col, outcome, formula=None):
    cols = ['index', 'coef', 'std_err', 't', 'p', 'ci_lower', 'ci_upper']
    res = model.summary2().tables[1].reset_index() 
    res.columns = cols
    if formula is None:
        pred = 'baseline'
        group = 'baseline'
        n_vars = 1
    else:   
        col = col.replace('bs(', '').replace(')', '')
        pred = col
        group = None
        # Determine which group pred belongs to
        for g in ['anthropometrics', 'exposures', 'sociodemographics']:
            if pred in groups_list_dict[g]:
                group = g
                break
        
        # get the last part of the formula after ' ~ ... + '
        n_vars = len(formula.split(' ~ ')[1].split(' + ')) if ' ~ ' in formula else 1
    res['group'] = group
    res['cov'] = pred
    res['val'] = None
    res['target'] = outcome
    res['mse'] = mse(model)
    res['rss'] = rss(model)
    res['n_vars'] = n_vars
    return res[['target', 'group', 'cov', 'val'] + cols + ['mse', 'rss', 'n_vars']]
    
def get_base_formula(df, outcome):
    variables = ["sex", "age", "sex*age"] 
    race_cols = [col for col in df.columns if col.startswith('race_') and not col.endswith('_adj')]
    variables += race_cols
    return f"{outcome} ~ {' + '.join(variables)}"

def stepwise_select(df, base_formula, predictors, loss=mse):
    base_formula_vars = base_formula.split(' ~ ')[1].split(' + ')
    predictors = [pred for pred in predictors if pred not in base_formula_vars]
    predictors = [pred for pred in predictors if not pred.endswith('_adj')]
    predictors = [pred for pred in predictors if not pred.endswith('_old')]
    selected, remaining, results = [], set(predictors), []
    best_score = np.inf

    while remaining:
        candidate, candidate_score = None, np.inf
        for var in remaining:
            new_vars = selected + [var]
            formula = base_formula + ' + ' + ' + '.join(new_vars)
            model = ols(formula, data=df).fit()
            score = loss(model)
            if score < candidate_score:
                candidate, candidate_score = var, score
        if candidate is None:
            break
        selected.append(candidate)
        remaining.remove(candidate)
        results.append({'formula_vars': selected.copy(), 'metric': candidate_score})
    return selected, results
    
def run_stepwise_df(cohort_name, df, outcome, groups, base_formula, reset=False):
    model = ols(base_formula, data=df).fit()
    result = get_params(model, 'Baseline', outcome)
    print(f"{'Base Formula':<20}: {base_formula}")

    for group in groups:
        print(f"{'Group':<20}: {group}")
        print('-'*100)
        if reset == 'cont':
            cols = get_predictors(df, groups)
        else:
            cols = get_predictors(df, [group])
        sel_vars, _ = stepwise_select(df, base_formula, cols, loss=mse)
        for idx in range(len(sel_vars)):
            use_vars = sel_vars[:idx + 1]
            # if var is continuous, can u add spline bs function?
            if df[use_vars[-1]].dtype == 'float64':
                use_vars[-1] = f'bs({use_vars[-1]}, df=3)'
            formula = base_formula + ' + ' + ' + '.join(use_vars)
            model = ols(formula, data=df).fit()
            try:
                formula_rhs = formula.split(' ~ ')[1]
                # try to isolate only after 'race_4 + ' if present
                if 'race_4 + ' in formula_rhs:
                    formula_rhs = formula_rhs.split('race_4 + ')[1]
            except Exception:
                formula_rhs = formula
            print(f"{'Covariates':<20}: {formula_rhs}")
            res = get_params(model, use_vars[-1], outcome, formula)
            result = pd.concat([result, res], ignore_index=True)
        if not reset:
            base_formula = formula
            break
        print('-'*100)
    return result

def fraction_calc(df, reset=None, num=10):
    def frac_calc(group):
        # Set the order for how the group appears - here, by 'n_vars' ascending, so group is in stepwise order
        race, target = group['index'].iloc[0], group['target'].iloc[0]
        baseline_row = df.query('index == @race and target == @target and group == "baseline"')
        baseline_coef = baseline_row['coef'].iloc[0]
        group = group.sort_values('n_vars', ascending=True).reset_index(drop=True)
        group['fraction'] = 100 * (baseline_coef - group['coef']) / baseline_coef
        rest_rows = group.iloc[1:].sort_values('n_vars', ascending=True).reset_index(drop=True)
        out = pd.concat([baseline_row, rest_rows], ignore_index=True)
        return out.fillna(0)
    cols = ['index', 'target', 'group'] if reset == 'reset' else ['target', 'index']
    df = df.groupby(cols, as_index=False).apply(lambda group: frac_calc(group).reset_index(drop=True)).reset_index(drop=True)
    return df

def process_df(df, ds, reset=None):
    df['cov'] = df['cov'].map(titles_dict)
    race_cols = [var for var in df['index'].unique() if var.startswith('race_') and not '_adj' in var]
    df = df.query('index in @race_cols')
    df['index'] = df['index'].str.replace('race_', '').str.split('[').str[0].map(int_race_dict)
    cols_to_round = [col for col in df.columns if col != 'p']
    df[cols_to_round] = df[cols_to_round].round(2)
    df = df.sort_values(['index', 'target', 'group'])
    if ds == 'ukb': df = df.query('index != "hispanic"')
    if ds == 'nh3': df = df.query('index != "asian"')
    df = fraction_calc(df, reset)
    df['cov'] = df['cov']
    df = df.rename(columns={'[0.025': 'ci_low', '0.975]': 'ci_high'})
    df = df.drop_duplicates(subset=['target', 'index', 'cov'], keep='first')
    return df[['n_vars', 'group', 'target', 'index', 'cov', 'fraction', 'coef', 'std_err', 't', 'p', 'mse']].reset_index(drop=True)

def main(args):
    output_types = [ 'reset', 'cont' ]
    for cohort in args.cohorts:
        results_dict = {suffix: pd.DataFrame() for suffix in output_types}
        
        for target in args.target_pfts:
            print('\n\n------------------------------------------------------------------------------------------------')
            print(f"{'Cohort':<20}: {cohort}")
            print(f"{'Target':<20}: {target}")
            df = pd.read_csv(f"{args.data_dir}/{cohort}/{cohort}_ref_d.csv")
            base_formula = get_base_formula(df, target)

            for suffix in output_types:
                res = run_stepwise_df(cohort, df, target, args.groups, base_formula, reset=suffix)
                results_dict[suffix] = pd.concat([results_dict[suffix], res], ignore_index=True)
        
        for suffix in output_types:
            results_dict[suffix].to_csv(f"{args.out_dir}/raw/{cohort}_{suffix}.csv", index=False)
        
        for suffix in output_types:
            processed = process_df(results_dict[suffix], cohort, reset=suffix)
            processed.to_csv(f"{args.out_dir}/processed/{cohort}_{suffix}.csv", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run stepwise regression explanation workflow.")
    parser.add_argument('--data_dir', default='../data/processed')
    parser.add_argument('--out_dir', default='../results/explain/tables')
    parser.add_argument('--cohorts', nargs='+', default=['nh3', 'nh4', 'nh', 'ukb'])
    parser.add_argument('--groups', nargs='+', default=['anthropometrics', 'sociodemographics', 'exposures'])
    parser.add_argument('--target_pfts', nargs='+', default=["fev1", "fvc"])
    args = parser.parse_args()
    main(args)
