import os
import pandas as pd
import numpy as np

raw_dir = "../data/raw"
processed_dir = "../data/processed"

# Choose the cohort to analyze
# For example: cohort = 'nhanes3'

# Load data (add missing os import and ensure variables exist)
exam_df = pd.read_csv(os.path.join(raw_dir, "nhanes/nh3_adult_youth_exam/nhanes3_exam_processed.csv"))
adult_df = pd.read_csv(os.path.join(raw_dir, "nhanes/nh3_adult_youth_exam/nhanes3_adult_processed.csv"))

cohort = "nh3"  
processed = pd.read_csv(os.path.join(processed_dir, f'{cohort}_ref_imp.csv'))
ids = processed['SEQN'].unique()

# Filter each DataFrame by SEQN in processed set
exam_filt = exam_df[exam_df['SEQN'].isin(ids)]
adult_filt = adult_df[adult_df['SEQN'].isin(ids)]

# Count NA in each column for filtered raw datasets
print("NA counts for exam data:")
print(exam_filt.isna().sum().reset_index().rename(columns={'index': 'column', 0: 'num_NA'}))

print("\nNA counts for adult data:")
print(adult_filt.isna().sum().reset_index().rename(columns={'index': 'column', 0: 'num_NA'}))

print("\nNA counts for spirometry data:")
print(spiro_filt.isna().sum().reset_index().rename(columns={'index': 'column', 0: 'num_NA'}))
