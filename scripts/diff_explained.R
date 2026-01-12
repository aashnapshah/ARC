# race_adjusted.R

# ---- Setup ----
packages <- c("binom","caret","dplyr","fst","gamlss","gamlss.dist",
              "ggplot2","GGally","ggbeeswarm","ggridges","ggsci",
              "glmnet","gridExtra","kableExtra","knitr","lubridate",
              "MASS","plotly","plyr","pROC","readr","RNHANES","rspiro",
              "scales","srvyr","stringr","survey","splines", "tidyr", "combinat")

package.check <- lapply(packages, function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, dependencies = TRUE)
  require(x, character.only = TRUE)
})

# ---- Helper functions ----
calculate_estimate <- function(model) {
    mse <- mean(residuals(model)^2)
    return(mse)
}

calculate_mse <- function(model) {
    mse <- mean(residuals(model)^2)
    return(mse)
}

calculate_RSS <- function(model) {
    mse <- sum(residuals(model)^2)
    return(mse)
}

# ---- Section 1: Stepwise add predictors by sex ----
cohorts <- c('nhanes', 'ukb')

for (i in seq_along(cohorts)){
  file_path <- sprintf("../../processed/%s/_healthy_encoded_2023-07-18.csv", cohorts[i])
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% dplyr::mutate(FEV1 = FEV1 * 1000, FVC = FVC * 1000)
  for (j in 1:2) {
    if (j == 1) {
      sex <- 'female'
      df_pft_sex <- cohort_df[cohort_df$RIAGENDR == 1, ]
    } else {
      sex <- 'male'
      df_pft_sex <- cohort_df[cohort_df$RIAGENDR == 0, ]
    }
    for (pft in c("FEV1", "FVC")) {
      categories <- c('Anthropometrics', 'Exposures',  'Sociodemographics', 'Ancestry', 'Genetics', 'PRS')
      if (cohorts[i] == 'nhanes3') {
        const_pred <- sprintf("%s ~ RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Other + RIDRETH_White", pft)
        pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
          dplyr::mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
        pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(df_pft_sex))
      } else if (cohorts[i] == 'nhanes4') {
        const_pred <- sprintf("%s ~ RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
        pred_map <- read.table("data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
          dplyr::mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
        pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(df_pft_sex))
      } else {
        const_pred <- sprintf("%s ~ RIDAGEYR + RIDRETH_Black + RIDRETH_Asian + RIDRETH_South.Asian + RIDRETH_Other + RIDRETH_White", pft)
        pred_map <- read.table("data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
          dplyr::mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
        pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(df_pft_sex))
      }

      pft_sex_race_df <- data.frame()
      formula <- as.formula(const_pred)
      model <- lm(formula, data = df_pft_sex)
      ESTIMATE <- calculate_estimate(model)
      df <- as.data.frame(summary(model)$coefficients)
      df <- cbind(Covariate = rownames(df), df)
      df['Model'] <- 'Baseline'
      df['Category'] <- 'Baseline'
      df['R2'] <- summary(model)$r.squared
      df['RSS'] <- ESTIMATE
      pft_sex_race_df <- rbind(pft_sex_race_df, df)

      for (cat in categories) {
        predictors <- as.character(pred_map[pred_map$Category == cat, ]$Predictor.Name.R)
        if (length(predictors) == 0) next
        selected_variables <- c()
        remaining_variables <- setdiff(predictors, selected_variables)
        best_EST <- Inf
        while (length(remaining_variables) > 0) {
          best_variable <- ""
          for (variable in remaining_variables) {
            new_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "), "+", variable)))
            new_EST <- calculate_estimate(new_model)
            if (new_EST <= best_EST) {
              best_variable <- variable
              best_EST <- new_EST
            }
          }
          selected_variables <- c(selected_variables, best_variable)
          model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
          ESTIMATE <- best_EST
          remaining_variables <- setdiff(remaining_variables, best_variable)
          df <- as.data.frame(summary(model)$coefficients)
          df <- cbind(Covariate = rownames(df), df)
          df['Model'] <- sprintf("+ %s", best_variable)
          df['Category'] <- cat
          df['RSS'] <- ESTIMATE
          df['R2'] <- summary(model)$r.squared
          pft_sex_race_df <- rbind(pft_sex_race_df, df)
          row.names(pft_sex_race_df) <- NULL
        }
        print(cat); print(cohorts[i]); print(sex); print(pft)
      }
      outfile <- sprintf('data/%s/%s_race_est_%s_%s_2023-15-06.csv', cohorts[i], cohorts[i], tolower(pft), sex)
      print(outfile)
      #write.csv(pft_sex_race_df, file = outfile, row.names = FALSE)
    }
  }
}

# ---- Section 2: Main execution, no race covariates ----
for (i in seq_along(cohorts)) {
  if (cohorts[i] == 'nhanes') {
    file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-07-18.csv", cohorts[i], cohorts[i])
  } else {
    file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-12-04.csv", cohorts[i], cohorts[i])
  }
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% dplyr::mutate(FEV1 = FEV1 * 1000, FVC = FVC * 1000)
  for (pft in c("FEV1", "FVC")) {
    pft_race_df <- data.frame()
    if (cohorts[i] == 'nhanes3') {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Other + RIDRETH_White", pft)
    } else if (cohorts[i] == 'nhanes4') {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    } else if (cohorts[i] == 'nhanes') {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    } else {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    }
    const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR", pft)
    formula <- as.formula(const_pred)
    model <- lm(formula, data = cohort_df)
    ESTIMATE <- calculate_estimate(model)
    df <- as.data.frame(summary(model)$coefficients)
    df <- cbind(Covariate = rownames(df), df)
    df['Model'] <- 'Baseline'
    df['Category'] <- 'Baseline'
    df['R2'] <- summary(model)$r.squared
    df['MSE'] <- ESTIMATE
    df['RSS'] <- calculate_RSS(model)
    pft_race_df <- rbind(pft_race_df, df)
    categories <- c('Anthropometrics', 'Exposures', 'Sociodemographics', 'PRS')
    for (cat in categories) {
      pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE, sep = ",") %>%
        dplyr::mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name))
      pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(cohort_df))
      predictors <- as.character(pred_map[pred_map$Category == cat, ]$Predictor.Name.R)
      predictors <- predictors[!grepl("_norm(al)?$|composite?$", predictors)]
      if (length(predictors) == 0) next
      selected_variables <- c()
      remaining_variables <- setdiff(predictors, selected_variables)
      best_EST <- Inf
      while (length(remaining_variables) > 0) {
        best_variable <- ""
        for (variable in remaining_variables) {
          new_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "), "+", variable)))
          new_EST <- calculate_estimate(new_model)
          if (new_EST <= best_EST) {
            best_variable <- variable
            best_EST <- new_EST
          }
        }
        selected_variables <- c(selected_variables, best_variable)
        model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
        ESTIMATE <- best_EST
        remaining_variables <- setdiff(remaining_variables, best_variable)
        df <- as.data.frame(summary(model)$coefficients)
        df <- cbind(Covariate = rownames(df), df)
        df['Model'] <- sprintf("+ %s", best_variable)
        df['Category'] <- cat
        df['MSE'] <- ESTIMATE
        df['R2'] <- summary(model)$r.squared
        df['RSS'] <- calculate_RSS(model)
        pft_race_df <- rbind(pft_race_df, df)
        row.names(pft_race_df) <- NULL
      }
    }
    print(cohorts[i])
    print(pft)
    outfile <- sprintf('data/%s/%s_race_est_c_%s_wo_race_2023-12-15.csv', cohorts[i], cohorts[i], tolower(pft))
    print(outfile)
    print(pft_race_df)
    write.csv(pft_race_df, file = outfile, row.names = FALSE)
  }
}

# ---- Section 3: Stepwise with all variables (category All), sorted, with race covariates ----
for (i in seq_along(cohorts)){
  if (cohorts[i] == 'nhanes') {
    file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-07-18.csv", cohorts[i], cohorts[i])
  } else {
    file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-12-04.csv", cohorts[i], cohorts[i])
  }
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% dplyr::mutate(FEV1 = FEV1*1000, FVC = FVC*1000)
  for (pft in c("FEV1", "FVC")) {
    pft_race_df <- data.frame()
    if (cohorts[i] == 'nhanes3') {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Other + RIDRETH_White", pft)
    } else if (cohorts[i] == 'nhanes4') {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    } else if (cohorts[i] == 'nhanes') {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    } else {
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    }
    formula <- as.formula(const_pred)
    model <- lm(formula, data = cohort_df)
    ESTIMATE <- calculate_estimate(model)
    df <- as.data.frame(summary(model)$coefficients)
    df <- cbind(Covariate = rownames(df), df)
    df['Model'] <- 'Baseline'
    df['Category'] <- 'Baseline'
    df['R2'] <- summary(model)$r.squared
    df['MSE'] <- ESTIMATE
    df['RSS'] <- calculate_RSS(model)
    pft_race_df <- rbind(pft_race_df, df)

    categories <- c('Anthropometrics', 'Exposures', 'Sociodemographics', 'PRS')
    pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE, sep = ",") %>%
      dplyr::mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name))
    pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(cohort_df))
    pred_map <- subset(pred_map, Category %in% categories)
    predictors <- as.character(pred_map$Predictor.Name.R)
    predictors <- sort(predictors[!grepl("BMXARML|(al)?$|composite?$", predictors)])
    selected_variables <- c()
    remaining_variables <- setdiff(predictors, selected_variables)
    best_EST <- Inf
    while (length(remaining_variables) > 0) {
      best_variable <- ""
      for (variable in remaining_variables) {
        new_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "), "+", variable)))
        new_EST <- calculate_estimate(new_model)
        best_variable <- variable
        best_EST <- new_EST
      }
      print(remaining_variables)
      selected_variables <- c(selected_variables, best_variable)
      model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
      ESTIMATE <- best_EST
      remaining_variables <- setdiff(remaining_variables, best_variable)
      df <- as.data.frame(summary(model)$coefficients)
      df <- cbind(Covariate = rownames(df), df)
      df['Model'] <- sprintf("+ %s", best_variable)
      df['Category'] <- "All"
      df['R2'] <- summary(model)$r.squared
      df['MSE'] <- ESTIMATE
      df['RSS'] <- calculate_RSS(model)
      pft_race_df <- rbind(pft_race_df, df)
      row.names(pft_race_df) <- NULL
    }
    print(cohorts[i])
    print(pft)
    print(pft_race_df)
    outfile <- sprintf('data/%s/%s_race_est_nc_sorted_w_race_%s_2023-12-15.csv', cohorts[i], cohorts[i], tolower(pft))
    print(outfile)
    print(pft_race_df)
    write.csv(pft_race_df, file = outfile, row.names = FALSE)
  }
}

# ---- Section 4: Per-race models for UKB ----
cohorts <- c('ukb')
for (i in seq_along(cohorts)){
  if (cohorts[i] == 'nhanes') {
    file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-07-18.csv", cohorts[i], cohorts[i])
  } else {
    file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-12-04.csv", cohorts[i], cohorts[i])
  }
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% dplyr::mutate(FEV1 = FEV1*1000, FVC = FVC*1000)
  for (pft in c("FEV1", "FVC")) {
    print('here')
    white_cohort <- cohort_df[cohort_df$RIDRETH_White == 1, ]
    black_cohort <- cohort_df[cohort_df$RIDRETH_Black == 1, ]
    asian_cohort <- cohort_df[cohort_df$RIDRETH_Asian == 1, ]
    #hispanic_cohort <- cohort_df[cohort_df$RIDRETH_Hispanic == 1, ]
    other_cohort <- cohort_df[cohort_df$RIDRETH_Other == 1, ]
    race_cohorts <- list(
      white = white_cohort,
      black = black_cohort,
      asian = asian_cohort,
      other = other_cohort
      # , hispanic = hispanic_cohort
    )
    pft_race_df <- data.frame()
    for (race in names(race_cohorts)) {
      cohort <- race_cohorts[[race]]
      const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR * RIDAGEYR", pft)
      formula <- as.formula(const_pred)
      model <- lm(formula, data = cohort)
      ESTIMATE <- calculate_estimate(model)
      df <- as.data.frame(summary(model)$coefficients)
      df <- cbind(Covariate = rownames(df), df)
      df['Model'] <- 'Baseline'
      df['Category'] <- 'Baseline'
      df['R2'] <- summary(model)$r.squared
      df['MSE'] <- ESTIMATE
      df['RSS'] <- calculate_RSS(model)
      df['Race'] <- race
      pft_race_df <- rbind(pft_race_df, df)
      categories <- c('Anthropometrics', 'PRS', 'Exposures', 'Sociodemographics')
      selected_variables <- c()
      for (cat in categories) {
        pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE, sep = ",") %>%
          dplyr::mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name))
        pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(cohort_df))
        pred_map <- subset(pred_map, Category %in% cat)
        predictors <- as.character(pred_map$Predictor.Name.R)
        predictors <- sort(predictors[!grepl("_norm(al)?$|composite?$", predictors)])
        # For each variable in category, add and refit model (no stepwise, just accumulate)
        for (variable in predictors) {
          selected_variables <- c(selected_variables, variable)
          updated_formula <- sprintf("%s + %s", const_pred, paste(selected_variables, collapse = " + "))
          model <- update(model, formula = as.formula(updated_formula))
          df <- as.data.frame(summary(model)$coefficients)
          df <- cbind(Covariate = rownames(df), df)
          df['Model'] <- variable
          df['Category'] <- cat
          df['R2'] <- summary(model)$r.squared
          df['MSE'] <- calculate_estimate(model)
          df['RSS'] <- calculate_RSS(model)
          df['Race'] <- race
          pft_race_df <- rbind(pft_race_df, df)
          row.names(pft_race_df) <- NULL
          print(paste("Cohort:", cohorts[i], "PFT:", pft, "Race:", race, "Category:", cat))
        }
      }
    }
    print(pft_race_df)
    outfile <- sprintf('data/%s/%s_race_variance_explained_%s_2023-12-15.csv', cohorts[i], cohorts[i], tolower(pft))
    write.csv(pft_race_df, file = outfile, row.names = FALSE)
    print(outfile)
    print(pft_race_df)
  }
}
