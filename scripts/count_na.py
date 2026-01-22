import os
import pandas as pd
import numpy as np
import re

# Set data locations
data_dir = "../data/raw"
processed_dir = "../data/processed"

# Load file paths
nhanes3_exam_path = os.path.join(data_dir, "nhanes/nh3_adult_youth_exam/nhanes3_exam_processed.csv")
nhanes3_youth_path = os.path.join(data_dir, "nhanes/nh3_adult_youth_exam/nhanes3_youth_processed.csv")
nhanes3_adult_path = os.path.join(data_dir, "nhanes/nh3_adult_youth_exam/nhanes3_adult_processed.csv")

# Use only SEQN values in processed nh3_ref file
# Infer the correct filename (latest): find file matching processed/nhanes3_*.csv
import glob
ref_file = sorted(glob.glob(os.path.join(processed_dir, "nh3/nh3_ref.csv")))[-1]
nh3_ref = pd.read_csv(ref_file)
ref_seqns = nh3_ref['seqn'].str.split('-').str[1].astype(int).values
print('-'*100)
print(nh3_ref['seqn'])
print(ref_seqns)
print('-'*100)

# Read raw data
exam = pd.read_csv(nhanes3_exam_path)
adult = pd.read_csv(nhanes3_adult_path)

print(exam.isna().sum().sort_values(ascending=True))
print('-'*100)
print(adult.isna().sum().sort_values(ascending=True))
print('-'*100)


mvn_dict = {
    'BMPHT': 88888,
    'BMPWT': 888888,
    'BMPBMI': 8888,
    'BMPARMC': 8888,
    'BMPARML': 8888,
    'BMPLEG': 8888,
    'BMPWAIST': 88888,
    'BMPSITHT': 88888
}

# Concatenate the full NHANES 3 cohort
adult = adult.rename(columns={'SEQN': 'seqn'})
exam = exam.rename(columns={'SEQN': 'seqn'})

sub_nh3 = adult[adult['seqn'].isin(ref_seqns)]
na_counts = sub_nh3.isna().sum().sort_values(ascending=True)

sub_nh3 = exam[exam['seqn'].isin(ref_seqns)]
na_counts = sub_nh3.isna().sum().sort_values(ascending=True)

