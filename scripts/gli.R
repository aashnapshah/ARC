library(gamlss)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(xgboost)

# Change directory
setwd("/n/groups/patel/aashna/PFT-ML-2024/")

# Load data
file_path = 'projects/reference_pfts/splits/UKB_train_ref_pft.csv'
ukb_train <- read.csv(file_path)

ukb_test <- read_csv('projects/reference_pfts/predictions/UKB_test_ref_pft_with_preds.csv', show_col_types = FALSE)
nhanes_test <- read_csv('projects/reference_pfts/predictions/NHANES_test_ref_pft_with_preds.csv', show_col_types = FALSE)
ukb_test_nonr <- read_csv('projects/reference_pfts/predictions/UKB_all_ref_pft_with_preds.csv', show_col_types = FALSE)
nhanes_test_nonr <- read_csv('projects/reference_pfts/predictions/NHANES_all_ref_pft_with_preds.csv', show_col_types = FALSE)

ukb_test$RIAGENDR <- as.integer(ifelse(ukb_test$RIAGENDR == "True", 1, 0))
nhanes_test$RIAGENDR <- as.integer(ifelse(nhanes_test$RIAGENDR == "True", 1, 0))
ukb_test_nonr$RIAGENDR <- as.integer(ifelse(ukb_test_nonr$RIAGENDR == "True", 1, 0))
nhanes_test_nonr$RIAGENDR <- as.integer(ifelse(nhanes_test_nonr$RIAGENDR == "True", 1, 0))

# Select the columns in ukb_train
common_columns <- intersect(names(ukb_train), names(ukb_test))
common_columns <- intersect(common_columns, names(nhanes_test))
common_columns <- c("BMXHT", "RIAGENDR", "RIDAGEYR", "RIDRETH", "FVC")
ukb_train <- ukb_train[, common_columns, drop = FALSE]
ukb_test_ <- ukb_test[, common_columns, drop = FALSE]
nhanes_test_ <- nhanes_test[, common_columns, drop = FALSE]
ukb_test_nonr_ <- ukb_test_nonr[, common_columns, drop = FALSE]
nhanes_test_nonr_ <- nhanes_test_nonr[, common_columns, drop = FALSE]

# Calculate ethnicity probabilities and weights
ethnicity_probabilities <- table(ukb_train$RIDRETH) / nrow(ukb_train)
ukb_train$RIDRETH_Prob <- ethnicity_probabilities[ukb_train$RIDRETH]

probabilities <- ukb_train %>%
  group_by(RIDRETH, RIAGENDR) %>%
  summarise(prob = n() / nrow(ukb_train)) %>%
  ungroup()

probabilities$weight <- 1 / probabilities$prob

ukb_train <- merge(ukb_train, probabilities, by = c("RIDRETH", "RIAGENDR"), all.x = TRUE)
weights <- ukb_train$weight %>% as.vector() %>% setNames(NULL)



my_gamlss <- function(PFT, f, sample_weights = NULL) {
  set.seed(100)
  use_formula <- as.formula(sprintf("%s ~ %s", PFT, f))
  print(use_formula)
  mod_all <- gamlss(
    formula = use_formula,
    sigma.formula = ~1+cs(RIDAGEYR, 2),
    nu.formula=~1,
    data=ukb_train,
    weights = sample_weights,
    family=BCCGo,
    method = RS()
  )
}

formulas <- c()
formulas["Age_Sex_HT"] <- "1 + RIAGENDR + log(BMXHT) + scs(log(RIDAGEYR), penalty=4) + RIAGENDR:scs(log(RIDAGEYR), penalty=4)"
formulas["Age_Sex_HT_Race"] <- "1 + RIAGENDR + RIDRETH + log(BMXHT) + scs(log(RIDAGEYR), penalty=4) + RIAGENDR:scs(log(RIDAGEYR), penalty=4)"  

# List of PFTs
pfts <- c('FVC') #, 'FEV1', 'FEV1/FVC')

# Process each PFT
for (pft in pfts) {
  new_col_neutral <- paste0("gli_neutral_refitted_", pft)
  new_col_race <- paste0("gli_race_refitted_", pft)
  
  # Train models
  model_race_free <- my_gamlss(pft, formulas$Age_Sex_HT, sample_weights = NULL)
  model_race <- my_gamlss(pft, formulas$Age_Sex_HT_Race, sample_weights = NULL)
  
  # Generate predictions for training data
  ukb_train[[new_col_neutral]] <- predict(model_race_free, newdata = ukb_train_, type = "response")
  ukb_train[[new_col_race]] <- predict(model_race, newdata = ukb_train_, type = "response")
  
  # Generate predictions for test data
  ukb_test_[[new_col_neutral]] <- predict(model_race_free, newdata = ukb_test_, type = "response")
  ukb_test_[[new_col_race]] <- predict(model_race, newdata = ukb_test_, type = "response")
  
  # Adjust and generate predictions for NHANES test data
  nhanes_test_$RIDRETH <- ifelse(nhanes_test_$RIDRETH == "Hispanic", "Other", nhanes_test_$RIDRETH)
  nhanes_test_[[new_col_neutral]] <- predict(model_race_free, newdata = nhanes_test_, type = "response")
  nhanes_test_[[new_col_race]] <- predict(model_race, newdata = nhanes_test_, type = "response")
  
  # Generate predictions for test data
  ukb_test_nonr_[[new_col_neutral]] <- predict(model_race_free, newdata = ukb_test_nonr_, type = "response")
  ukb_test_nonr_[[new_col_race]] <- predict(model_race, newdata = ukb_test_nonr_, type = "response")
  
  # Adjust and generate predictions for NHANES test data
  nhanes_test_nonr_$RIDRETH <- ifelse(nhanes_test_nonr_$RIDRETH == "Hispanic", "Other", nhanes_test_nonr_$RIDRETH)
  nhanes_test_nonr_[[new_col_neutral]] <- predict(model_race_free, newdata = nhanes_test_nonr_, type = "response")
  nhanes_test_nonr_[[new_col_race]] <- predict(model_race, newdata = nhanes_test_nonr_, type = "response")
}

# Save updated dataframes
#write_csv(ukb_test, 'projects/reference_pfts/splits/UKB_test_ref_pft_with_preds_gli.csv')
#write_csv(nhanes_test, 'projects/reference_pfts/splits/NHANES_test_ref_pft_with_preds_gli.csv')

#write_csv(nhanes_test, 'projects/reference_pfts/predictions/NHANES_test_pft_with_preds_gli.csv')
#write_csv(ukb_test, 'projects/reference_pfts/predictions/UKB_test_ref_pft_with_preds_gli.csv')
#write_csv(nhanes_test_nonr, 'projects/reference_pfts/predictions/NHANES_all_pft_with_preds_gli.csv')
#write_csv(ukb_test_nonr, 'projects/reference_pfts/predictions/UKB_all_ref_pft_with_preds_gli.csv')

df_black <- ukb_test_ %>% filter(RIDRETH=="Black")
print(mean(abs(df_black$gli_neutral_refitted_FVC - df_black$FVC)))

# Save updated datasets to CSV
#write.csv(ukb_train, "../reference_pfts/UKB_train_gli_ref_pft.csv", row.names = FALSE)
#write.csv(ukb_test, "processed/UKB_test_gli_ref_pft.csv", row.names = FALSE)
#write.csv(nhanes_test, "processed/NHANES_test_gli_ref_pft.csv", row.names = FALSE)

df_train_old <- read_csv('/n/groups/patel/aashna/PFT-ML-2024/projects/reference_pfts/old/predictions/UKB_train_gli_ref_pft.csv')
print(sum(ukb_train$gli_neutral_refitted_FVC == df_train_old$gli_neutral_refitted_FVC)==nrow(ukb_train))
print(sum(ukb_train$RIAGENDR == df_train_old$RIAGENDR)==nrow(ukb_train))
print(sum(ukb_train$RIDRETH == df_train_old$RIDRETH)==nrow(ukb_train))
print(sum(ukb_train$FVC == df_train_old$FVC)==nrow(ukb_train))
print(sum(ukb_train$BMXHT == df_train_old$BMXHT)==nrow(ukb_train))
print(sum(ukb_train$RIDAGEYR == df_train_old$RIDAGEYR)==nrow(ukb_train))

df_test_old <- read_csv('/n/groups/patel/aashna/PFT-ML-2024/projects/reference_pfts/old/predictions/UKB_test_gli_ref_pft.csv')
df_black <- df_test_old %>% filter(RIDRETH=="Black")
print(mean(abs(df_black$gli_neutral_refitted_FVC - df_black$FVC)))

df_test_nhanes_old <- read_csv('/n/groups/patel/aashna/PFT-ML-2024/projects/reference_pfts/old/predictions/NHANES_test_gli_ref_pft.csv')
df_black <- df_test_nhanes_old %>% filter(RIDRETH=="Black")
print(mean(abs(df_black$gli_neutral_refitted_FVC - df_black$FVC)))

