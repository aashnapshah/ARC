import os
import argparse
import pandas as pd
import numpy as np
import joblib
from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
import sys

sys.path.append('../')
from data.dicts.col_dict import race_dict, int_race_dict, titles_dict

def recode_race(df):
    his_val = race_dict['hispanic']
    other_val = race_dict['other']
    df['race_c'] = df['race'].where(df['race'] != his_val, other_val)
    return df

def combine_smoke_exposure(df):
    cols = df.columns[df.columns.str.contains('smoke', case=False)]
    df['smoke_exposure'] = df[cols].any(axis=1).astype(int)
    return df

def _xgboost_grid(X, y):
    """XGBoost with grid search."""
    param_grid = {
        "n_estimators": [50, 100, 200],
        "max_depth": [3, 5, 7],
        "learning_rate": [0.1, 0.01],
        "subsample": [1]
    }
    xgb = XGBRegressor(objective='reg:squarederror', verbosity=1, n_jobs=-1)
    cv = GridSearchCV(xgb, param_grid, cv=3, scoring='neg_root_mean_squared_error', verbose=0)
    cv.fit(X, y)
    return cv.best_estimator_

def _xgboost(X_train, y_train, sample_weights=None):
    model = XGBRegressor(objective='reg:squarederror', n_jobs=-1)
    model.fit(X_train, y_train, sample_weight=sample_weights)
    return model

def metrics(y_true, y_pred):
    return {
        "n": len(y_true),
        "mae": f"{mean_absolute_error(y_true, y_pred):.2f}",
        "mse": f"{mean_squared_error(y_true, y_pred):.2f}",
        "r2": f"{r2_score(y_true, y_pred):.2f}"
    }

def evaluate(model, X_test, y_test, group_col=None, groups=None):
    y_pred = model.predict(X_test)
    results = []
    overall_metrics = metrics(y_test, y_pred)
    overall_metrics.update({"group": "All"})
    results.append(overall_metrics)
    unique_groups = np.unique(groups)
    for val in unique_groups:
        mask = groups == val
        group_metrics = metrics(y_test[mask], y_pred[mask])
        group_metrics.update({"group": int_race_dict[str(val)]})
        results.append(group_metrics)
    return pd.DataFrame(results)[['group', 'n', 'mae', 'mse', 'r2']]

def process_data(data_dir, cohort):
    path = os.path.join(data_dir, f'{cohort}/{cohort}_ref.csv')
    df = pd.read_csv(path)
    print('-'*100)
    cols = ['sex', 'race', 'age', 'fev1', 'fvc', 'fev1_fvc', 'height', 'sit_height', 'waist_circ', 'immigrant', 'hs_grad', 'smoke_exposure']
    print(df[cols].head())
    print('-'*100)
    df = recode_race(df)
    df = combine_smoke_exposure(df)
    return df

def get_data_dict(data_dir, features, targets, race_col, train_cohort='ukb', test_size=0.3, seed=42):
    df_dict = {k: process_data(data_dir, k) for k in ['ukb', 'nh', 'nh3', 'nh4']}
    all_cohorts = list(df_dict.keys())
    train_df_full = df_dict[train_cohort]
    test_cohorts = [c for c in all_cohorts if c != train_cohort]
    train_df_full = train_df_full.dropna(subset=features + targets)
    for cohort in test_cohorts:
        df_dict[cohort] = df_dict[cohort].dropna(subset=features + targets)

    stratify_col = train_df_full[race_col] if race_col in train_df_full.columns else None
    train_df, val_df = train_test_split(
        train_df_full,
        test_size=test_size,
        random_state=seed,
        stratify=stratify_col
    )
    tests = {f'{train_cohort}_val': val_df}
    for cohort in test_cohorts:
        tests[cohort] = df_dict[cohort]
    return train_df, tests

def run_experiment(train_df, test_dfs, features, target, race_col, save_path=None, seed=42):
    all_results = pd.DataFrame()

    for i in range(len(features)):
        for race_adj in [True, False]:
            base = ['age', 'sex', 'race_c'] if race_adj else ['age', 'sex']
            curr_feats = base + features[:i+1] 
            X_train, y_train = train_df[curr_feats].values, train_df[target].values 

            model = _xgboost(X_train, y_train)

            for name, df in test_dfs.items():
                X_test = df[curr_feats].values
                y_test = df[target].values
                groups = df[race_col].values
                eval_df = evaluate(model, X_test, y_test, groups=groups)
                eval_df["set"] = name
                eval_df["target"] = target
                eval_df["n_cov"] = len(curr_feats)
                eval_df["cov"] = titles_dict[curr_feats[-1]]
                eval_df["race_adj"] = race_adj
                all_results = pd.concat([all_results, eval_df], ignore_index=True)

    return all_results


def main(args):
    train_cohort = ['ukb', 'nh', 'nh3', 'nh4']
    feats = ['height', 'sit_height', 'waist_circ', 'immigrant', 'inc_pov_bin', 'hs_grad', 'smoke_exposure']
    for train_cohort in train_cohort:
        train_df, test_dfs = get_data_dict(args.data_dir, features=feats, targets=args.targets, race_col=args.race_col, train_cohort=train_cohort)
        all_targets_results = []
        for target in args.targets:
            results = run_experiment(train_df, test_dfs, features=feats, target=target,
                                    race_col=args.race_col)
            all_targets_results.append(results)

        summary_df = pd.concat(all_targets_results, ignore_index=True)
        print(summary_df)
        save_path = os.path.join(args.save_path, f'metrics_train_{train_cohort}.csv') 
        summary_df.to_csv(save_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_dir", type=str, default="../data/processed/")
    parser.add_argument("--save_path", type=str, default="../results/train/tables")
    parser.add_argument("--targets", nargs="+", default=["fvc", "fev1"])
    parser.add_argument("--race_col", type=str, default="race")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    main(args)
