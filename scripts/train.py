import os
import argparse
import pandas as pd
import numpy as np
import joblib
from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression

import statsmodels.api as sm
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

def train_gli_global(train_df, features, target, include_race):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, Formula
    from rpy2.robjects.packages import importr
    from rpy2.robjects.conversion import localconverter
    
    pandas2ri.activate()
    base = importr('base')
    gamlss = importr('gamlss')
    df = train_df.copy()
    df['log_age'] = np.log(df['age'])
    df['log_height'] = np.log(df['height'])
    df['sex'] = df['sex'].astype(int)
    if include_race:
        race_levels = list(int_race_dict.keys())
        df['race_c'] = df['race_c'].astype(int).astype(str)
        df['race_c'] = pd.Categorical(df['race_c'], categories=race_levels)
    if not include_race:
        weight_df = (
            df.groupby(['sex', 'race_c'])
            .size()
            .div(len(df))
            .reset_index(name='prob')
        )
        weight_df['ipw'] = 1 / weight_df['prob']
        df = df.merge(weight_df, on=['sex', 'race_c'], how='left')
        df['ipw'] = df['ipw'].fillna(1.0)
        print(df)
    else:
        df['ipw'] = df.get('ipw', 1.0)
    with localconverter(robjects.default_converter + pandas2ri.converter):
        rdf = robjects.conversion.py2rpy(df)
    rhs_terms = ["sex", "log_height", "scs(log_age, penalty=4)", "sex:scs(log_age, penalty=4)"]
    if include_race:
        rhs_terms.insert(2, "race_c")
    formula_str = "{} ~ {}".format(target, " + ".join(rhs_terms))
    sigma_formula_str = "~ 1 + cs(log_age, 2)"
    nu_formula_str = "~1"
    robjects.r('library(gamlss)')
    robjects.r('library(splines)')
    gamlss_func = robjects.r['gamlss']
    w = np.asarray(df["ipw"], dtype=float)
    w_r = robjects.FloatVector(w)

    fitted_model = gamlss_func(
        Formula(formula_str),
        data=rdf,
        sigma_formula=Formula(sigma_formula_str),
        nu_formula=Formula(nu_formula_str),
        weights=w_r,
        family="BCCGo",
    )
    return fitted_model

def gli_global_predict(model, df, include_race):
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter
    import numpy as np

    X = df[['age', 'sex', 'height']].copy()
    X['log_age'] = np.log(X['age'])
    X['log_height'] = np.log(X['height'])
    X['sex'] = X['sex'].astype(int)
    if include_race:
        race_levels = list(int_race_dict.keys())
        X['race_c'] = df['race_c'].astype(int).astype(str)
        X['race_c'] = pd.Categorical(X['race_c'], categories=race_levels)

    pandas2ri.activate()
    with localconverter(ro.default_converter + pandas2ri.converter):
        rX = ro.conversion.py2rpy(X)

    r_predict = ro.r['predict']
    yhat = r_predict(model, newdata=rX, what="mu", type="response")

    return np.asarray(yhat, dtype=float)

def _xgboost_grid(X, y):
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
        "mae": mean_absolute_error(y_true, y_pred),
        "mse": mean_squared_error(y_true, y_pred),
        "r2": r2_score(y_true, y_pred)
    }

def ci_low_high(arr, alpha=0.05):
    lower = np.percentile(arr, 100*alpha/2)
    upper = np.percentile(arr, 100*(1-alpha/2))
    return lower, upper

def bootstrap_metrics(y_true, y_pred, n_bootstraps=500, seed=42):
    n = len(y_true)
    rng = np.random.RandomState(seed)
    maes, mses, r2s = [], [], []
    for _ in range(n_bootstraps):
        idx = rng.choice(np.arange(n), size=n, replace=True)
        y_true_bs = y_true[idx]
        y_pred_bs = y_pred[idx]
        maes.append(mean_absolute_error(y_true_bs, y_pred_bs))
        mses.append(mean_squared_error(y_true_bs, y_pred_bs))
        try:
            r2s.append(r2_score(y_true_bs, y_pred_bs))
        except Exception:
            r2s.append(np.nan)
    mae_mean, mae_ci = np.mean(maes), ci_low_high(maes)
    mse_mean, mse_ci = np.mean(mses), ci_low_high(mses)
    r2_mean, r2_ci  = np.nanmean(r2s), ci_low_high([x for x in r2s if not np.isnan(x)])
    return {
        "mae": f"{mae_mean:.2f}",
        "mae_ci_lower": f"{mae_ci[0]:.2f}",
        "mae_ci_upper": f"{mae_ci[1]:.2f}",
        "mse": f"{mse_mean:.2f}",
        "mse_ci_lower": f"{mse_ci[0]:.2f}",
        "mse_ci_upper": f"{mse_ci[1]:.2f}",
        "r2": f"{r2_mean:.2f}",
        "r2_ci_lower": f"{r2_ci[0]:.2f}",
        "r2_ci_upper": f"{r2_ci[1]:.2f}",
        "n": n
    }

def evaluate(model, X_test, y_test, group_col=None, groups=None, n_bootstraps=500, seed=42):
    y_pred = model.predict(X_test)
    results = []
    overall_stats = bootstrap_metrics(y_test, y_pred, n_bootstraps=n_bootstraps, seed=seed)
    overall_stats.update({"group": "All"})
    results.append(overall_stats)
    unique_groups = np.unique(groups)
    for val in unique_groups:
        mask = groups == val
        stats = bootstrap_metrics(y_test[mask], y_pred[mask], n_bootstraps=n_bootstraps, seed=seed)
        stats.update({"group": int_race_dict[str(val)]})
        results.append(stats)
    return pd.DataFrame(results)[[
        'group', 'n',
        'mae', 'mae_ci_lower', 'mae_ci_upper',
        'mse', 'mse_ci_lower', 'mse_ci_upper',
        'r2', 'r2_ci_lower', 'r2_ci_upper'
    ]]

def process_data(data_dir, cohort):
    path = os.path.join(data_dir, f'{cohort}/{cohort}_ref.csv')
    df = pd.read_csv(path)
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

def evaluate_gli(model, test_df, target, race_col, n_bootstraps=500, seed=42, include_race=True):
    y_true = test_df[target].values
    y_pred = gli_global_predict(model, test_df, include_race)
    groups = test_df[race_col].values
    results = []
    overall_stats = bootstrap_metrics(y_true, y_pred, n_bootstraps=n_bootstraps, seed=seed)
    overall_stats.update({"group": "All"})
    results.append(overall_stats)
    unique_groups = np.unique(groups)
    for val in unique_groups:
        mask = groups == val
        stats = bootstrap_metrics(y_true[mask], y_pred[mask], n_bootstraps=n_bootstraps, seed=seed)
        stats.update({"group": int_race_dict[str(val)]})
        results.append(stats)
    return pd.DataFrame(results)[[
        'group', 'n',
        'mae', 'mae_ci_lower', 'mae_ci_upper',
        'mse', 'mse_ci_lower', 'mse_ci_upper',
        'r2', 'r2_ci_lower', 'r2_ci_upper'
    ]]

def run_experiment(train_df, test_dfs, features, target, race_col, save_path=None, seed=42, n_bootstraps=500):
    all_results = pd.DataFrame()
    # Only fit GLI global for the "canonical" covariate sets: age, sex, height [+ race]
    # (do not waste compute fitting GLI on other sets)
    for race_adj in [True, False]:
        # GLI global as before
        base_feats = ['age', 'sex', 'height']
        if race_adj:
            base_feats += ['race_c']
        # XGBoost: incrementally add features and re-evaluate
        # e.g., ["age"], ["age", "sex"], ["age", "sex", "height"], etc. (+ optional extra features passed in 'features')
        feat_names = ['age', 'sex', 'height'] + [f for f in features if f not in ['age', 'sex', 'height', 'race_c']]
        if race_adj and 'race_c' not in feat_names:
            feat_names.append('race_c')
        # For each incrementally larger feature set
        for i in range(1, len(feat_names) + 1):
            xgb_feats = feat_names[:i]
            # Always ensure 'race_c' is last if race_adj
            if race_adj and 'race_c' not in xgb_feats:
                xgb_feats.append('race_c')
            # Skip if only race_c in features (nonsense)
            if xgb_feats == ['race_c']:
                continue
            X_train = train_df[xgb_feats].values
            y_train = train_df[target].values
            xgb_model = _xgboost(X_train, y_train)
            for name, df in test_dfs.items():
                # Only test on data that has required features available
                if not all([f in df.columns for f in xgb_feats]):
                    continue
                X_test = df[xgb_feats].values
                y_test = df[target].values
                groups = df[race_col].values
                eval_df_xgb = evaluate(xgb_model, X_test, y_test, groups=groups, n_bootstraps=n_bootstraps, seed=seed)
                eval_df_xgb["set"] = name
                eval_df_xgb["target"] = target
                eval_df_xgb["n_cov"] = len(xgb_feats)
                eval_df_xgb["cov"] = feat_names[i-1]
                eval_df_xgb["race_adj"] = race_adj
                eval_df_xgb["model"] = "xgboost"
                all_results = pd.concat([all_results, eval_df_xgb], ignore_index=True)
        # Fit/fetch GLI global model and evaluate on all test sets (as before)
        model = train_gli_global(train_df, ['age', 'sex', 'height'], target, include_race=race_adj)
        for name, df in test_dfs.items():
            eval_df_gli = evaluate_gli(model, df, target, race_col, n_bootstraps=n_bootstraps, seed=seed, include_race=race_adj)
            eval_df_gli["set"] = name
            eval_df_gli["target"] = target
            eval_df_gli["n_cov"] = 3 + int(race_adj)
            eval_df_gli["cov"] = "GLI global"
            eval_df_gli["race_adj"] = race_adj
            eval_df_gli["model"] = "gli_global"
            all_results = pd.concat([all_results, eval_df_gli], ignore_index=True)
    return all_results

def main(args):
    train_cohorts = ['ukb', 'nh', 'nh3', 'nh4']
    feats = ['height', 'sit_height', 'waist_circ', 'immigrant', 'inc_pov_bin', 'hs_grad', 'smoke_exposure']
    for train_cohort in train_cohorts:
        train_df, test_dfs = get_data_dict(args.data_dir, features=feats, targets=args.targets, race_col=args.race_col, train_cohort=train_cohort)
        all_targets_results = []
        for target in args.targets:
            results = run_experiment(train_df, test_dfs, features=feats, target=target,
                                    race_col=args.race_col, seed=args.seed, n_bootstraps=args.n_bootstraps)
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
    parser.add_argument("--n_bootstraps", type=int, default=1000)
    args = parser.parse_args()
    main(args)
