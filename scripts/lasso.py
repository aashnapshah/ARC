import pandas as pd
import numpy as np
import sys
import argparse
sys.path.append('../')
from data.dicts.col_dict import *
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from explain import get_predictors

def lasso_regression(df, target, predictors, alpha=0.01, random_state=42):
    X = df[predictors].copy()
    y = df[target].copy()
    X = X.fillna(0)
    non_na_idx = ~y.isna()
    X = X.loc[non_na_idx]
    y = y.loc[non_na_idx]
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=random_state)
    lasso = Lasso(alpha=alpha, random_state=random_state, max_iter=10000)
    lasso.fit(X_train, y_train)
    y_pred = lasso.predict(X_test)
    feature_importance = pd.Series(lasso.coef_, index=predictors)
    r2 = r2_score(y_test, y_pred)
    return feature_importance, r2

def main(args):
    cohorts = ['nh3', 'nh4', 'nh', 'ukb']
    targets = ['fev1', 'fvc']
    groups = args.groups if hasattr(args, 'groups') else ['anthropometrics', 'sociodemographics', 'exposures']
    
    for cohort in cohorts:
        results = []
        df = pd.read_csv(f'../data/processed/{cohort}/{cohort}_ref_d.csv')
        predictors = get_predictors(df, groups)
        for target in targets:
            feature_importance, r2 = lasso_regression(df, target, predictors)
            out_df = pd.DataFrame({
                'cohort': cohort,
                'target': target,
                'feature': feature_importance.index,
                'importance': feature_importance.values,
                'r2': r2
            })
            results.append(out_df)
        results_df = pd.concat(results, ignore_index=True)
        results_df.to_csv(f'../results/lasso/feat_imp_{cohort}.csv', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run lasso regression feature importance.")
    parser.add_argument('--groups', nargs='+', default=['anthropometrics', 'sociodemographics', 'exposures'])
    args = parser.parse_args()
    main(args)