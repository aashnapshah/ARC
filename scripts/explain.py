# Cleaned and Dynamic Refactor of race_adjusted.R Logic (Python)

import os
import pandas as pd
import numpy as np
import argparse
import sys

sys.path.append('../')
from data.dicts.col_dict import cat_col_dict

from statsmodels.formula.api import ols

def mse(model):
    return np.mean(model.resid**2)

def rss(model):
    return np.sum(model.resid**2)

def ensure_dir(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)

def get_pred_map(csv_path, avail_cols, cats=None):
    pm = pd.read_csv(csv_path)
    pm['Predictor.Name.R'] = pm['Predictor.Name'].str.replace(" ", ".", regex=False)
    pm = pm[pm['Predictor.Name.R'].isin(avail_cols)]
    if cats:
        pm = pm[pm['Category'].isin(cats)]
    return pm

def stepwise_select(df, base_formula, predictors, fit_metric=mse):
    selected, remaining, results = [], set(predictors), []
    best_score = np.inf

    while remaining:
        candidate, candidate_score = None, np.inf
        for var in remaining:
            new_vars = selected + [var]
            formula = base_formula + ' + ' + ' + '.join(new_vars)
            model = ols(formula, data=df).fit()
            score = fit_metric(model)
            if score < candidate_score:
                candidate, candidate_score = var, score
        if candidate is None:
            break
        selected.append(candidate)
        remaining.remove(candidate)
        results.append({'formula_vars': selected.copy(), 'metric': candidate_score})
    return selected, results

def save_coef_summary(model, model_desc, cat, extra=None, fit_metric=None):
    summary = model.summary2().tables[1].reset_index().rename(columns={'index': 'Covariate'})
    summary['Model'] = model_desc
    summary['Category'] = cat
    summary['R2'] = model.rsquared
    summary['MSE'] = mse(model)
    summary['RSS'] = rss(model)
    if fit_metric:
        summary['MSE'] = fit_metric
    if extra:
        for k, v in extra.items():
            summary[k] = v
    return summary

def get_base_formula(df, pft):
    variables = ["sex", "age", "sex*age"] 
    race_cols = [col for col in df.columns if col.startswith('race_') and not col.endswith('_adj')]
    variables += race_cols
    return f"{pft} ~ {' + '.join(variables)}"
    
def run_stepwise_df(cohort_name, df, pft, categories, base_formula, exclude_pat=None):
    """Stepwise by adding predictors"""
    output_df = pd.DataFrame()
    model = ols(base_formula, data=df).fit()
    output_df = pd.concat([output_df, save_coef_summary(model, 'Baseline', 'Baseline', fit_metric=mse(model))], ignore_index=True)

    for cat in categories:
        preds = cat_col_dict[cat] 
        preds = [col for col in df.columns if col in preds]
        sel_vars, _ = stepwise_select(df, base_formula, preds, fit_metric=mse)
        
        for idx in range(len(sel_vars)):
            use_vars = sel_vars[:idx + 1]
            base_formula = base_formula + ' + ' + ' + '.join(use_vars)
            print(base_formula)
            model = ols(base_formula, data=df).fit()
            output_df = pd.concat(
                [output_df, save_coef_summary(model, '+ ' + use_vars[-1], cat, fit_metric=mse(model))],
                ignore_index=True,
            )
    return output_df

def main(args):
    for cohort in args.cohorts:
        if cohort == 'nh':
            nh3_path = f"{args.data_dir}/nh3_ref_encoded.csv"
            nh4_path = f"{args.data_dir}/nh4_ref_encoded.csv"
            nh3_ref = pd.read_csv(nh3_path)
            nh4_ref = pd.read_csv(nh4_path)
            nhanes_cols = list(set(nh3_ref.columns).intersection(nh4_ref.columns))
            df = pd.concat([nh3_ref, nh4_ref])[nhanes_cols]
        else:   
            file_path = f"{args.data_dir}/{cohort}_ref_encoded.csv"
            df = pd.read_csv(file_path)
        print(df.columns)
        for pft in args.target_pfts:
            base_formula = get_base_formula(df, pft)
            print(f"Stepwise: {base_formula}")
            _ = run_stepwise_df(cohort, df, pft, args.categories, base_formula)
            print(_)
            
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run stepwise regression explanation workflow.")
    parser.add_argument('--data_dir', default='../data/processed')
    parser.add_argument('--cohorts', nargs='+', default=['nh3', 'nh4', 'nh', 'ukb'])
    parser.add_argument('--categories', nargs='+', default=['Anthropometrics', 'Sociodemographics', 'Exposures'])
    parser.add_argument('--target_pfts', nargs='+', default=["fev1", "fvc", "fev1_fvc"])
    args = parser.parse_args()
    main(args)


# # ---- Section 3: Stepwise with all variables (category All), sorted, with race covariates ----
# for cohort in COHORTS:
#     file_path = f"../../processed/{cohort}/{cohort}_healthy_encoded_2023-07-18.csv" if cohort == "nhanes" else f"../../processed/{cohort}/{cohort}_healthy_encoded_2023-12-04.csv"
#     df = pd.read_csv(file_path)
#     pred_map_path = "../../data/predictor_category_map.csv"
#     for pft in ["FEV1", "FVC"]:
#         pft_df = pd.DataFrame()
#         base_formula = base_formula_by_cohort(cohort, pft, race=True)
#         # Add all eligible variables in sorted order, not greedy stepwise (mimic R)
#         categories = ['Anthropometrics', 'Exposures', 'Sociodemographics', 'PRS']
#         pred_map = get_pred_map(pred_map_path, df.columns, categories)
#         all_vars = clean_predictors(pred_map, exclude_pattern="BMXARML|(al)?$|composite?$")
#         selected_vars = []
#         pft_df = pd.concat([pft_df, save_coef_summary(ols(base_formula, data=df).fit(), 'Baseline', 'Baseline')], ignore_index=True)
#         for var in all_vars:
#             selected_vars.append(var)
#             formula = base_formula + ' + ' + ' + '.join(selected_vars)
#             model = ols(formula, data=df).fit()
#             pft_df = pd.concat(
#                 [pft_df, save_coef_summary(model, '+ ' + var, 'All')], ignore_index=True
#             )
#         out_path = f'data/{cohort}/{cohort}_race_est_nc_sorted_w_race_{pft.lower()}_2023-12-15.csv'
#         ensure_dir(out_path)
#         pft_df.to_csv(out_path, index=False)
#         print(f"Stepwise-sorted-w-race: {out_path}")
