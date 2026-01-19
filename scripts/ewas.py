import pandas as pd
import numpy as np
import os
from patsy import dmatrices
import statsmodels.api as sm
from statsmodels.stats.weightstats import DescrStatsW


def run_ewas(df, targets):
    results = []
    for outcome in targets:
        for col in df.columns:
            if col in ['age', 'sex']+targets:
                continue
            formula = f"{outcome} ~ age + sex + {col}"
            y, X = dmatrices(formula, data=df, return_type='dataframe')
            model = sm.OLS(y, X).fit()
            coef = model.summary2().tables[1].reset_index()  # columns: index, Coef., Std.Err., t, P>|t|, [0.025, 0.975]
            coef.rename(columns={'index':'Covariate'}, inplace=True)
            coef['Target'] = outcome
            coef['Model'] = col
            results.append(coef)
    
    result_df = pd.concat(results, ignore_index=True)
    print(result_df)
    return result_df

def main():
    # Data and variable settings
    targets = ['fev1', 'fvc', 'fev1_fvc']
    data_dir = "../data/processed"
    output_dir = "../results/ewas/"

    nh3_ref = pd.read_csv(f"{data_dir}/nhanes3/nh3_ref_scaled.csv")
    nh4_ref = pd.read_csv(f"{data_dir}/nhanes4/nh4_ref_scaled.csv")
    nhanes_cols = list(set(nh3_ref.columns).intersection(nh4_ref.columns))

    nh_ref = pd.concat([nh3_ref, nh4_ref])[nhanes_cols]
    ukb_ref = pd.read_csv(f"{data_dir}/ukb/ukb_ref_scaled_encoded.csv")

    # Run standard EWAS for NHANES and UKB
    cohort_names = {'nh3_ref': nh3_ref, 'nh4_ref': nh4_ref, 'nh_ref': nh_ref, 'ukb_ref': ukb_ref}
    for name, cohort in cohort_names.items():
        df = read_and_prepare(cohort, targets)
        result = run_ewas(df, targets)
        result.to_csv(f"{output_dir}/{name}_ewas.csv", index=False)

if __name__ == "__main__":
    main()