# Load necessary libraries
library(tidyverse)
library(xgboost)
library(caret)

# Set working directory
setwd("/n/groups/patel/aashna/PFT-ML-2024/")

# Load datasets
train_data <- read_csv('projects/reference_pfts/splits/UKB_train_ref_pft.csv', show_col_types = FALSE)
ukb_test_data <- read_csv('projects/gli_global/predictions/UKB_test_ref_pft_gli_preds.csv', show_col_types = FALSE)
nhanes_test_data <- read_csv('projects/gli_global/predictions/NHANES_test_ref_pft_gli_preds.csv') #splits/NHANES_test_ref_pft_.csv', show_col_types = FALSE)
non_ref_ukb <- read_csv('projects/gli_global/predictions/UKB_all_ref_pft_gli_preds.csv')
non_ref_nhanes <- read_csv('projects/gli_global/predictions/NHANES_all_ref_pft_gli_preds.csv', show_col_types = FALSE)

process_data <- function(data, cohort){
  if (cohort == 'NHANES') {
    data$oldRIDRETH <- data$RIDRETH
    data$RIDRETH_new <- ifelse(data$RIDRETH == "Hispanic", "Other", data$RIDRETH)
    if (!('DMD_SMOKE_EXP' %in% colnames(data))){
      data <- data %>% mutate(DMD_Work_Smoke = OCQ_Work_Smoke) #%>% select(-OCQ_Work_Smoke)
    }
  }
  if (cohort == 'UKB') {
    if (!('DMD_Born_US' %in% colnames(data))){
      data$DMD_Born_US <- data$DMD_Born_UK
    }
    data$RIDRETH_new <- data$RIDRETH
  }  
  
  if (!('DMD_SMOKE_EXP' %in% colnames(data))){
    data <- data %>%
      mutate(DMD_SMOKE_EXP = DMD_Home_Smoke | DMD_Work_Smoke)
  }
  
  data$Smoke_Exposure <- as.integer(data$DMD_SMOKE_EXP)
  data <- data %>% mutate(RIDRETH_new = factor(RIDRETH, levels = c("White", "Black", "Asian", "Hispanic", "Other")))
  data$RIDRETH_encoded <- as.numeric(data$RIDRETH_new)
  
  #drop columns that start with pred_
  data <- data %>% select(-starts_with("pred_"))
  return(data)
}

train_data <- process_data(train_data, 'UKB')
ukb_test_data <- process_data(ukb_test_data, 'UKB')
nhanes_test_data <- process_data(nhanes_test_data, 'NHANES')
non_ref_ukb <- process_data(non_ref_ukb, 'UKB')
non_ref_nhanes <- process_data(non_ref_nhanes, 'NHANES')

# Function to train XGBoost model
xgboost_train <- function(train_matrix, train_labels, file_name) {
  tune_grid <- expand.grid(
    nrounds = c(50, 100, 200),
    max_depth = c(3, 5, 7),
    eta = c(0.1, 0.01),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )
  
  train_control <- trainControl(
    method = "cv",
    number = 3,
    verboseIter = TRUE,
    returnData = FALSE,
    returnResamp = "all",
    allowParallel = TRUE
  )
  
  xgb_train <- train(
    x = train_matrix,
    y = train_labels,
    trControl = train_control,
    tuneGrid = tune_grid,
    method = "xgbTree",
    metric = "RMSE", 
  )
  
  best_model <- xgb_train$finalModel
  saveRDS(best_model, file = file_name)
  
  return(best_model)
}

# Define the models and corresponding feature sets
models <- list(
  bh = c("BMXHT", "RIDAGEYR", "RIAGENDR"),
  bhs = c("RIAGENDR", "RIDAGEYR", "BMXHT", "BMXSIT"),
  bhws = c("RIAGENDR", "RIDAGEYR", "BMXHT", "BMXWAIST", "BMXSIT"),
  bhwsu = c("RIAGENDR", "RIDAGEYR", "BMXHT", "BMXWAIST", "BMXSIT", "DMD_Born_US"),
  bhwsfus = c("RIAGENDR", "RIDAGEYR", "BMXHT", "BMXWAIST", "BMXSIT", "DMD_Born_US", "Smoke_Exposure")
)

# Prepare and evaluate models
for (pft in c('FVC', 'FEV1', 'FEV1/FVC')) {
  for (race in c('w', 'wo')) {
    for (model_name in names(models)) {
      features <- models[[model_name]]
        if (race == 'w'){
          features <- c(features, 'RIDRETH_encoded')
        }

    # Prepare training data
    train_matrix <- as.matrix(train_data %>% select(all_of(features)))
    train_labels <- train_data[[pft]]
    
    print(pft)
    
    if (pft == 'FEV1/FVC') {
      print('not')
      pft_name= 'FEV1.FVC'
    }
    else {
      pft_name = pft
    }
    
    file_name <- paste0('models/', model_name, '_', pft_name, '_', race, '_race_xgb_model.rds')
    if (file.exists(file_name)) {
      print('model exists')
      best_model <- readRDS(file_name)
    } else {
      # Train the model
      best_model <- xgboost_train(train_matrix, train_labels, file_name)
    }
    # Prepare test data
    ukb_test_matrix <- as.matrix(ukb_test_data %>% select(all_of(features)))
    ukb_test_labels <- ukb_test_data[[pft_name]]
    
    nhanes_test_matrix <- as.matrix(nhanes_test_data %>% select(all_of(features)))
    nhanes_test_labels <- nhanes_test_data[[pft_name]]
    
    # Predict on UKB test data
    ukb_preds <- predict(best_model, newdata = ukb_test_matrix)
    ukb_test_data[[paste0("pred_", pft_name, "_", model_name, "_", race, "_race")]] <- ukb_preds
    
    # Predict on NHANES test data
    nhanes_preds <- predict(best_model, newdata = nhanes_test_matrix)
    nhanes_test_data[[paste0("pred_", pft_name, "_", model_name, "_", race, "_race")]] <- nhanes_preds
    
    non_ref_ukb_preds <- predict(best_model, newdata = as.matrix(non_ref_ukb %>% select(all_of(features))))
    non_ref_ukb[[paste0("pred_", pft_name, "_", model_name, "_", race, "_race")]] <- non_ref_ukb_preds
    
    non_ref_nhanes_preds <- predict(best_model, newdata = as.matrix(non_ref_nhanes %>% select(all_of(features))))
    non_ref_nhanes[[paste0("pred_", pft_name, "_", model_name, "_", race, "_race")]] <- non_ref_nhanes_preds
    
    # Calculate RMSE
    ukb_mae <- mean(abs(ukb_preds - ukb_test_labels))
    nhanes_mae <- mean(abs(nhanes_preds - nhanes_test_labels))
    
    print(model_name)
    print(ukb_mae)
    print(nhanes_mae)
    }
  }
}

nhanes_test_data$RIDRETH <- nhanes_test_data$oldRIDRETH
non_ref_nhanes$RIDRETH <- non_ref_nhanes$oldRIDRETH
# Save updated dataframes
write_csv(ukb_test_data, 'projects/reference_pfts/predictions/UKB_test_ref_pft_preds.csv')
write_csv(nhanes_test_data, 'projects/reference_pfts/predictions/NHANES_test_ref_pft_preds.csv')
write_csv(non_ref_nhanes, 'projects/reference_pfts/predictions/NHANES_all_ref_pft_preds.csv')
write_csv(non_ref_ukb, 'projects/reference_pfts/predictions/UKB_all_ref_pft_preds.csv')



