import warnings
warnings.filterwarnings("ignore")

import argparse
import pandas as pd
import numpy as np
import xgboost as xgb
import matplotlib.pyplot as plt
import json
import joblib
import os

from collections import defaultdict
from itertools import product
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression  
from sklearn.experimental import enable_iterative_imputer  
from sklearn.impute import IterativeImputer   

DMARETHN_MAP = {
    1: "Non-Hispanic white",
    2: "Non-Hispanic black",
    3: "Mexican American",
    8: "Other"
}

RIDRETH_MAP = {
    1: "Non-Hispanic white",
    2: "Non-Hispanic black",
    3: "Mexican American",
    4: "Asian",
    5: "Other",
} 

NHANES_III_MAP = {
        'HSSEX': 'RIAGENDR',
        'HSAGEIR': 'RIDAGEYR',
        'BMPLEG': 'BMXLEG',
        'BMPHTIN': 'BMXHT',
        'DMARACER': 'RIDRETH',
        'RIDRETH1': 'RIDRETH',
        'BMPARML': 'BMXARML',
        
    }

def to_numeric(X):
    for col in X:
        if col.lower() == 'seqn':
            continue
        col_lower = col.lower()
        if col_lower in ["riagendr", "is_female", "sex"]:
            if set(X[col].dropna().unique()).intersection({1, 2}):  # e.g., NHANES: 1=Male, 2=Female
                X[col] = (X[col] == 2).astype(int)
            else:
                X[col] = X[col].astype(str).str.strip().str.lower().map({'female': 1, 'male': 0})
        elif X[col].dtype == "object" or pd.api.types.is_categorical_dtype(X[col]):
            X[col] = pd.Categorical(X[col]).codes
        elif pd.api.types.is_bool_dtype(X[col]):
            X[col] = X[col].astype(int)
    return X

    
def nhanes4(path, features, target, race_col):
    all_cols = set(features) | {target, race_col, "SEQN"}
    # bmx_df = pd.read_csv(f"{path}/BMX_E.csv", usecols=lambda c: c in all_cols)
    # demo_df = pd.read_csv(f"{path}/DEMO_E.csv", usecols=lambda c: c in all_cols)
    # df = demo_df.merge(bmx_df, on="SEQN", how="inner").rename(columns={"RIDRETH1": "RIDRETH"})
    path = f"../data/raw/nh4/nh4_2023-07-18.csv"
    nh_ref = pd.read_csv(path, sep="\t", usecols=lambda c: c in all_cols).rename(columns={"SEQN": "seqn"})
    return nh_ref

def process_data(path, features, target, race_col=None, include_race=False):
    if 'nhanes3' in path:
        df = pd.read_csv(path, na_values=["", " ", "  ", "   ", "    ", "     "], keep_default_na=True)
        df = df.rename(columns=NHANES_III_MAP)
    elif 'nhanes4' in path:
        df = nhanes4(path, features, target, race_col)
    else:
        raise ValueError(f"Invalid path: {path}")
    df = df.replace([-9999, -99999, -999999, -9999999, -99999999, -999999999, 8888, 88888, 9999, 99999], np.nan)
    df = df.rename(columns=NHANES_III_MAP)
    if include_race and race_col is not None:
        features = features + [race_col]
    df = to_numeric(df)
    cols = features + [target]
    cols = [col for col in cols if col in df.columns]
    df = df.dropna(subset=cols)
    return df

def compute_ipw(race_id):
    vc = race_id.value_counts(normalize=True)
    weights = race_id.map(lambda x: 1.0 / vc[x])
    return weights

def evaluate(y_true, y_pred, groups=None, params=None):
    results = {}
    results["All"] = {
        "mae": mean_absolute_error(y_true, y_pred), 
        "mse": mean_squared_error(y_true, y_pred), 
        "r2": r2_score(y_true, y_pred), 
        "n": len(y_true)}
    
    for grp in np.unique(groups):
        mask = groups == grp
        group_size = int(np.sum(mask))
        results[DMARETHN_MAP.get(grp, str(grp))] = {
            "mae": mean_absolute_error(y_true[mask], y_pred[mask]),
            "mse": mean_squared_error(y_true[mask], y_pred[mask]),
            "r2": r2_score(y_true[mask], y_pred[mask]),
            "n": group_size,
        }
    return results

def _linreg(X_train, y_train, X_test, sample_weights=None):
    model = LinearRegression()
    model.fit(X_train, y_train, sample_weight=sample_weights)
    return model

def _xgboost(X_train, y_train, X_test, sample_weights=None):
    model = xgb.XGBRegressor(n_jobs=-1)
    model.fit(X_train, y_train, sample_weight=sample_weights)
    return model

def _mice(X_train, y_train, X_test, sample_weights=None):
    imp = IterativeImputer(sample_posterior=False, max_iter=10)
    X_train_imp = imp.fit_transform(X_train)
    X_test_imp = imp.transform(X_test)
    model = LinearRegression()
    model.fit(X_train_imp, y_train, sample_weight=sample_weights)
    return model

def run_experiment(
    nh3, nh4, 
    features, target, race_col, race_adj=False, seed=42,
    test_size=0.3,
    balance_mode="ipw",
    mfunc=_xgboost   
):
    features = features + [race_col] if race_adj else features
    # set age to 18-65
    nh3['bin_sh'] = pd.qcut(nh3['BMPSITHT'], q=4, labels=False, duplicates='drop')
    nh3['black'] = nh3['RIDRETH'] == 2
    print(nh3['bin_sh'].value_counts())
    print(nh3['black'].value_counts())
    strat_cols = ['bin_sh', 'black']
    df_train, df_test = train_test_split(nh3, test_size=test_size, random_state=seed, stratify=nh3[strat_cols])
    sample_weights = compute_ipw(df_train[race_col])
    X_train, y_train = df_train[features], df_train[target]
    X_test, y_test = df_test[features], df_test[target]
    model = mfunc(X_train, y_train, X_test, sample_weights)
    y_pred = model.predict(X_test)
    metrics = evaluate(y_test, y_pred, groups=df_test[race_col])    
    imputed_y = model.predict(nh4[features])
    pred_df = nh4[['seqn',race_col]] 
    pred_df['imputed_value'] = imputed_y
    return metrics, pred_df, model  

def summ_stats(df):
    df = df.rename(columns=NHANES_III_MAP)
    print(f"{'Variable':<18} {'Unique/Min':<15} {'Max':<10} {'Mean':<12} {'Std':<12} {'Median':<12}")
    print("-" * 75)
    print(f"{'Race':<18} {str(df['RIDRETH'].unique())}")
    print(f"{'Sex':<18} {str(df['RIAGENDR'].unique())}")
    #print(f"{'Age':<18} {df['RIDAGEYR'].min():<15.2f} {df['RIDAGEYR'].max():<10.2f} {df['RIDAGEYR'].mean():<12.2f} {df['RIDAGEYR'].std():<12.2f}")
    print(f"{'Leg Length':<18} {df['BMXLEG'].min():<15.2f} {df['BMXLEG'].max():<10.2f} {df['BMXLEG'].mean():<12.2f} {df['BMXLEG'].std():<12.2f}")
    print(f"{'Height':<18} {df['BMXHT'].min():<15.2f} {df['BMXHT'].max():<10.2f} {df['BMXHT'].mean():<12.2f} {df['BMXHT'].std():<12.2f}")
    if "BMPSITHT" in df.columns:
        print(f"{'Sitting Height':<18} {df['BMPSITHT'].min():<15.2f} {df['BMPSITHT'].max():<10.2f} {df['BMPSITHT'].mean():<12.2f} {df['BMPSITHT'].std():<12.2f} {df['BMPSITHT'].median():<12.2f}")
    
def main(args):      
    print('IMPUTE PATH: ', args.impute_path)
    nh3 = process_data(args.train_path, args.features, args.target, args.race_col)    
    nh4 = process_data(args.impute_path, args.features, args.target, args.race_col)
    # nh4 = nh4[nh4['RIDAGEYR'] >= 18]
    # nh3 = nh3[nh3['RIDAGEYR'] >= 18]
    summ_stats(nh3)
    summ_stats(nh4)
    
    m_funcs = [eval(m) for m in args.models]
    pairs = list(product(m_funcs, args.race_adj))

    metrics_rows = []
    imputed_rows = []
    
    models_dir = os.path.join(args.save_path, "models")
    os.makedirs(models_dir, exist_ok=True)

    for m_func, race_adj in pairs:
        metrics, pred_df, model = run_experiment(
            nh3, nh4, args.features, args.target, args.race_col, 
            race_adj=race_adj, balance_mode="ipw", mfunc=m_func
        )
        model_name = m_func.__name__
 
        pred_df['model'] = model_name
        pred_df['race_adj'] = race_adj
        imputed_rows.append(pred_df)
        if not race_adj:
            preds = pred_df['imputed_value']
            print('-'*100)
            print('Model: ', model_name, 'Race Adj: ', race_adj)
            print(preds.mean(), preds.std(), preds.min(), preds.max())
            print(preds.groupby(pred_df['RIDRETH']).mean())
         #   print(preds.groupby(pred_df['RIDRETH']).std())
            # print(preds.groupby(pred_df['RIDRETH']).min()
            # print(preds.groupby(pred_df['RIDRETH']).max())
        for key, value in metrics.items():
            metric_row = {"model": model_name, "race_adj": race_adj, "race": key}
            for k, v in value.items():
                metric_row.update({k: v})
            
            metrics_rows.append(metric_row)

    # Save CSVs
    df = pd.concat(imputed_rows)
    df.to_csv(f"{args.save_path}/tables/imputed_results.csv", index=False)
    pd.DataFrame(metrics_rows).to_csv(f"{args.save_path}/tables/metrics.csv", index=False)
    print(pd.DataFrame(metrics_rows))
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--train_path", type=str, default="../data/raw/nhanes/nh3_adult_youth_exam/nhanes3_exam_processed.csv")
    parser.add_argument("--impute_path", type=str, default="../data/raw/nhanes/nhanes4")
    parser.add_argument("--save_path", type=str, default="../results/imputed_sh")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--features", type=list, default=['RIDAGEYR', 'RIAGENDR', 'BMXLEG', 'BMXHT', 'BMXARML'])
    parser.add_argument("--target", type=str, default="BMPSITHT")
    parser.add_argument("--race_col", type=str, default="RIDRETH")
    parser.add_argument("--models", type=list, default=["_linreg", "_xgboost"])
    parser.add_argument("--race_adj", type=list, default=[False, True])
    args = parser.parse_args()
    main(args)