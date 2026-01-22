# +
import pandas as pd
import numpy as np

import seaborn as sns
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score

import datetime
import pickle
import os

import warnings
warnings.filterwarnings("ignore")

def map_subregion(df):
    country_map = pd.read_csv('data/continents2.csv')[['name', 'sub-region']]
    
    country_to_subregion_1 = dict(zip(country_map['name'], country_map['sub-region']))
    
    # create a dictionary mapping countries to subregions
    country_to_subregion_2 = {
        'Akrotiri and Dhekelia': 'Western Asia',
        'Borneo': 'Southeast Asia',
        'Bosnia and Herzegovina': 'Southern Europe',
        'Brunei': 'Southeast Asia',
        'Cape Verde': 'Western Africa',
        'Caribbean': 'Latin America and the Caribbean',
        'Channel Islands': 'Northern Europe',
        'East Africa': 'Eastern Africa',
        'Falkland Islands': 'South America',
        'Ivory Coast': 'Western Africa',
        'Madeira': 'Southern Europe',
        'Myanmar (Burma)': 'Southeast Asia',
        'Netherlands Antilles': 'Latin America and the Caribbean',
        'North Korea': 'Eastern Asia',
        'Pacific Islands': 'Oceania',
        'Palestine': 'Western Asia',
        'Republic of Kosovo': 'Southern Europe',
        'Saint Helena': 'Western Africa',
        'Serbia/Montenegro': 'Southern Europe',
        'Sicily': 'Southern Europe',
        'Swaziland': 'Southern Africa',
        'The Guianas': 'Latin America and the Caribbean',
        'UK': 'Northern Europe',
        'USA': 'Northern America',
        'West Indies': 'Latin America and the Caribbean',
        'Sri Lanka': 'Southern Asia',
        'Bulgaria': 'Eastern Europe',
        'Slovenia': 'Southern Europe',
        'Macedonia': 'Southern Europe'
    }
    
    merged_dict = country_to_subregion_1 | country_to_subregion_2

    
    # map the countries to subregions using the dictionary
    df['DMD_BORN_SubRegion'] = df['DMD_Born_UK'].map(merged_dict)


    return df

column_dict = {
    'RIDAGEYR': 'Age',
    'RIAGENDR': 'Female',
    'RIDRETH_White': 'White',
    'RIDRETH_Black': 'Black',
    'RIDRETH_Asian': 'Asian',
    'RIDRETH_South Asian': 'South Asian',
    'RIDRETH_Other': 'Other',
    'BMXHT': 'Height (cm)',
    'BMXWT': 'Weight (kg)',
    'BMXBMI': 'BMI',
    'BMXWAIST': 'Waist Circumference (cm)',
    'BMXHIP': 'Hip Circumference (cm)',
    'BMXSIT': 'Sitting Height (cm)',
    'DMD_Mom_Smoked': 'Mother Smoked During Pregnancy',
    'DMD_Finish_HS': 'Finished High School',
    'DMD_Home_Smoke': 'Smoking in Home',
    'DMD_Veteran': 'Veteran Status',
    'DMD_PM10': 'Exposure to PM10',
    'DMD_PM2.5': 'Exposure to PM2.5',
    'DMD_PM2.5_10': 'Exposure to PM2.5 and PM10',
    'DMD_Work_Fumes': 'Worked with Fumes',
    'DMD_Work_Exhaust': 'Worked with Exhaust',
    'DMD_Work_Smoking_Exp': 'Workplace Smoking Exposure',
    'DMD_Work_Breathing_Probs': 'Breathing Problems at Work',
    'BMXWT10_About average': 'Weight at 10 Y: About Average',
    'BMXWT10_Do not know': 'Weight at 10 Y: Do Not Know',
    'BMXWT10_Plumper': 'Weight at 10 Y: Plumper',
    'BMXWT10_Prefer not to answer': 'Weight at 10 Y: Prefer Not to Answer',
    'BMXWT10_Thinner': 'Weight at 10 Y: Thinner',
    'BMXHT10_About average': 'Height at 10 Y: About Average',
    'BMXHT10_Do not know': 'Height at 10 Y: Do Not Know',
    'BMXHT10_Prefer not to answer': 'Height at 10 Y: Prefer Not to Answer',
    'BMXHT10_Shorter': 'Height at 10 Y: Shorter',
    'BMXHT10_Taller': 'Height at 10 Y: Taller',
    'DMD_INC_Do not know': 'Household Income: Do Not Know',
    'DMD_INC_Four or more': 'Household Income: 52,000 - 100,000',
    'DMD_INC_None': 'Household Income: None',
    'DMD_INC_One': 'Household Income: < 18,000',
    'DMD_INC_Prefer not to answer': 'Household Income: Prefer Not to Answer',
    'DMD_INC_Three': 'Household Income: 31,000 - 51,999',
    'DMD_INC_Two': 'Household Income: 18,000 - 30,999',
    'DMD_BORN_SubRegion_Australia and New Zealand': 'Born: Australia and New Zealand',
    'DMD_BORN_SubRegion_Central Asia': 'Born: Central Asia',
    'DMD_BORN_SubRegion_Eastern Africa': 'Born: Eastern Africa',
    'DMD_BORN_SubRegion_Eastern Asia': 'Born: Eastern Asia',
    'DMD_BORN_SubRegion_Eastern Europe': 'Born: Eastern Europe',
    'DMD_BORN_SubRegion_Latin America and the Caribbean': 'Born: Latin America and the Caribbean',
    'DMD_BORN_SubRegion_Melanesia': 'Born: Melanesia',
    'DMD_BORN_SubRegion_Northern Africa': 'Born: Northern Africa',
    'DMD_BORN_SubRegion_Northern America': 'Born: Northern America', 
    'DMD_BORN_SubRegion_Northern Europe': 'Born: Northern Europe',
    'DMD_BORN_SubRegion_Oceania': 'Born: Oceania',
    'DMD_BORN_SubRegion_Polynesia': 'Born: Polynesia',
    'DMD_BORN_SubRegion_South America': 'Born: South America',
    'DMD_BORN_SubRegion_South-eastern Asia': 'Born: South-Eastern Asia',
    'DMD_BORN_SubRegion_Southeast Asia': 'Born: Southeast Asia',
    'DMD_BORN_SubRegion_Southern Africa': 'Born: Southern Africa',
    'DMD_BORN_SubRegion_Southern Asia': 'Born: Southern Asia',
    'DMD_BORN_SubRegion_Southern Europe': 'Born: Southern Europe',
    'DMD_BORN_SubRegion_Sub-Saharan Africa': 'Born: Sub-Saharan Africa',
    'DMD_BORN_SubRegion_Western Africa': 'Born: Western Africa',
    'DMD_BORN_SubRegion_Western Asia': 'Born: Western Asia',
    'DMD_BORN_SubRegion_Western Europe': 'Born: Western Europe'}

def one_hot_encoding(df, features, target):
    X = df[features]
    y = df[target]
    # Get a list of categorical columns
    categorical_columns = list(X.select_dtypes(include=['object']).columns)
    
    # Convert categorical columns to one-hot encoded columns
    for column in categorical_columns:
        one_hot = pd.get_dummies(X[column], prefix=column)
        X = X.drop(column, axis=1)
        X = X.join(one_hot)
    
    # Map original column names to descriptive names
    descriptive_cols = [column_dict[col] for col in X.columns]
    
    # Scale features using StandardScaler
    scaler = StandardScaler()
    X_scaled = pd.DataFrame(scaler.fit_transform(X))
    #X_scaled.columns=descriptive_cols
    X_scaled.columns = X.columns
    X_scaled.index=X.index
    
    return X_scaled, y

def train_test_split_stratified(df, features, target):    
    X_scaled, y = one_hot_encoding(df, features, target)

    X_train, X_test, y_train, y_test = train_test_split(X_scaled, y, 
                                                        stratify=X_scaled.join(ukb['RIDRETH'])['RIDRETH'],
                                                        test_size=0.2, random_state=42)
    
    return X_train, X_test, y_train, y_test, X_scaled.columns
    
def bootstrap(X_train, y_train, n_bootstraps, penalization):

    # Initialize array to store bootstrap coefficients
    bootstrap_coefs = np.zeros((n_bootstraps, X_train.shape[1]))

    # Perform bootstrap resampling and Lasso regression on training data
    for i in range(n_bootstraps):
        # Generate bootstrap sample
        indices = np.random.choice(X_train.index, size=X_train.shape[0], replace=True)
        X_boot = X_train.loc[indices]
        y_boot = y_train[indices]
        # Fit Lasso regression to bootstrap sample
        if penalization == 'Ridge':
            model = Ridge(alpha=0.01)
            
        elif penalization == 'Lasso':
            model = Lasso(alpha=0.01)
            
        else:
            model = LinearRegression()
            
        model.fit(X_boot, y_boot)
        bootstrap_coefs[i, :] = model.coef_
        
    # Calculate standard errors and confidence intervals
    z_score = 1.96  # 95% confidence interval
    ci_beta = z_score*np.std(bootstrap_coefs, axis=0, ddof=1)

    return ci_beta

def train_regression(df, features, target, penalization=Lasso, n_boot = 5):
    X_train, X_test, y_train, y_test, encoded_features = train_test_split_stratified(df, features, target)

    # Define the hyperparameter grid
    param_grid = {'alpha': [0.01, 0.1, 1, 10, 100]}

    if penalization == 'Ridge':
        # Perform grid search with cross-validation
        model = Ridge()
        grid_search = GridSearchCV(model, param_grid, cv=5, scoring='neg_mean_absolute_error')
        grid_search.fit(X_train, y_train)

        # Train the model with the best hyperparameters
        best_alpha = grid_search.best_params_['alpha']
        model = Ridge(alpha=best_alpha)
        model.fit(X_train, y_train)
        mean_score = -grid_search.best_score_
        
    elif penalization == "Lasso":
        model = Lasso()
        grid_search = GridSearchCV(model, param_grid, cv=5, scoring='neg_mean_absolute_error')
        grid_search.fit(X_train, y_train)

        # Train the model with the best hyperparameters
        best_alpha = grid_search.best_params_['alpha']
        model = Lasso(alpha=best_alpha)
        model.fit(X_train, y_train)
        mean_score = -grid_search.best_score_        
        
    else:
        model = LinearRegression()
        scores = cross_val_score(model, X_train, y_train, cv=5, scoring='neg_mean_absolute_error')
        mean_score = np.mean(scores)
        best_alpha = 0
        model.fit(X_train, y_train)
    
    ci_beta = bootstrap(X_train, y_train, n_bootstraps=n_boot, penalization = penalization)
    
    model_dict = {}
    model_dict[target] = {
                'Best Hyperparameters': best_alpha,
                'Features': encoded_features, 
                'Model': model,
                'CI Betas': ci_beta,
                'X_test': X_test,
                'y_test': y_test
                 }
    
    return model_dict

def update_model_dict(model, key, m_dict):
    if key in models:
        models[key].update(model_dict) # Append to the existing list
    else:
        models[key] = model_dict.copy()
    return models 

def sort_dict_by_key(d):
    return dict(sorted(d.items(), key=lambda x: x[0]))

def sort_dict_by_value(d):
    return dict(sorted(d.items(), key=lambda x: x[1], reverse=True))

def custom_sort_key(key):
    if key.startswith('RI'):
        return 4
    elif key.startswith('BMX'):
        return 3
    elif key.startswith('DMD'):
        return 2
    else:
        return 4
        
def sorted_feature_imp(features, model_coef, ci_coefs):    
    d = dict(zip(features, model_coef))
    d_ci = dict(zip(features, ci_coefs))
    
    filtered_dict = {k: v for k, v in d.items() if abs(v) > 0.01}
    sorted_dict = sort_dict_by_value(filtered_dict)
    sorted_keys = sorted(sorted_dict.keys(), key=lambda k: (custom_sort_key(k), abs(sorted_dict[k])), reverse=True)
    
    sorted_dict = {k: sorted_dict[k] for k in sorted_keys}

    new_dict = {column_dict[k]: v for k, v in sorted_dict.items()}
    
    # Order d_ci based on sorted_keys
    d_ci = {k: d_ci[k] for k in sorted_keys}

    return list(new_dict.keys()), list(new_dict.values()), list(d_ci.values())

# Define a function to calculate the mean absolute error and confidence intervals
def calculate_mae_ci(y_true, y_pred):
    mae = mean_absolute_error(y_true, y_pred)
    se = np.sqrt(np.mean((y_true - y_pred)**2))
    n = len(y_true)
    t = stats.t.ppf(1-0.025, n-1)
    ci = t * se / np.sqrt(n)
    return mae, ci