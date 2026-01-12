---
title: "race_adjusted"
author: "Aashna Shah"
date: "2023-06-14"
output: html_document
---

```{r setup, include=FALSE}
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
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```{r}
# Function to calculate the RSS
calculate_estimate <- function(model) {
    mse <- mean(residuals(model)^2)
  return(mse)
}
cohorts <- c('nhanes', 'ukb')

for (i in seq_along(cohorts)){
  file_path <- sprintf("../../processed/%s/_healthy_encoded_2023-07-18.csv", cohorts[i])
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% mutate(FEV1 = FEV1*1000, FVC = FVC*1000)
  for (j in 1:2) {
  if (j == 1) {
    sex <- 'female'
    df_pft_sex <- cohort_df[cohort_df$RIAGENDR==1,]
    }
  
  else {
    sex <- 'male'
    df_pft_sex <- cohort_df[cohort_df$RIAGENDR==0,]
  }
  # NOT SEX
  for (pft in c("FEV1", "FVC")){
    categories <- c('Anthropometrics', 'Exposures',  'Sociodemographics', 'Ancestry', 'Genetics', 'PRS')

    if (cohorts[i] == 'nhanes3') {
      const_pred <- sprintf("%s ~ RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Other + RIDRETH_White", pft)
      pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
                    mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
      pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(df_pft_sex))
    }
    else if (cohorts[i] == 'nhanes4') {
      const_pred <- sprintf("%s ~ RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
      pred_map <- read.table("data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
                    mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
      pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(df_pft_sex))
    }
    else {
      const_pred <- sprintf("%s ~ RIDAGEYR + RIDRETH_Black + RIDRETH_Asian + RIDRETH_South.Asian + RIDRETH_Other + RIDRETH_White", pft)
      pred_map <- read.table("data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
                    mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
      pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(df_pft_sex))
    }
    pft_sex_race_df = data.frame()
    # Initial model with only the intercept
    formula <- as.formula(const_pred)
    model <- lm(formula, data = df_pft_sex)
    
    # Calculate the initial RSS
    ESTIMATE <- calculate_estimate(model)
    
    df = as.data.frame(summary(model)$coefficients)
    df <- cbind(Covariate = rownames(df), df)
    df['Model'] <- 'Baseline'
    df['Category'] <- 'Baseline'
    df['R2'] <- summary(model)$r.squared
    df['RSS'] <- ESTIMATE
    pft_sex_race_df <- rbind(pft_sex_race_df,df)
          
    for (cat in categories) {
          predictors <- as.character(pred_map[pred_map$Category==cat, ]$Predictor.Name.R)
          if (length(predictors)==0) {next}
          # Create an empty list to store the selected variables
          selected_variables <- c()
          
          # Number of remaining variables
          remaining_variables <- setdiff(predictors, selected_variables)
          
          best_EST <- Inf
          while (length(remaining_variables) > 0) {
            # Initialize variables for tracking the best addition
            best_variable <- ""
            
            # Loop through the remaining variables and find the one that reduces the RSS the most
            for (variable in remaining_variables) {
              new_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "), "+", variable)))
              new_EST <- calculate_estimate(new_model)
                
              if (new_EST <= best_EST) {
                #print(new_EST)
                #print(best_EST)
                best_variable <- variable
                best_EST <- new_EST
              }
            }
            
          # Add the best variable to the selected variables
          selected_variables <- c(selected_variables, best_variable)
          
          # Update the model and RSS with the added variable
          model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
          ESTIMATE <- best_EST

          # Remove the added variable from the remaining variables
          remaining_variables <- setdiff(remaining_variables, best_variable)
          df = as.data.frame(summary(model)$coefficients)
          df <- cbind(Covariate = rownames(df), df)
          df['Model'] <- sprintf("+ %s", best_variable)
          df['Category'] <- cat
          df['RSS'] <- ESTIMATE
          df['R2'] <- summary(model)$r.squared

          pft_sex_race_df <- rbind(pft_sex_race_df,df)
          row.names(pft_sex_race_df) <- NULL
    
          }
    # Create the final model with the selected variables
    #final_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
    
    # Check the variable order and intercept
    #selected_variables <- names(final_model$coefficients)[-1]
    #final_intercept <- coef(final_model)[[1]]
    
    # Print the results
      print(cat)
      print(cohorts[i])
      print(sex)
      print(pft)

    }
    outfile = sprintf('data/%s/%s_race_est_%s_%s_2023-15-06.csv', cohorts[i], cohorts[i], tolower(pft), sex)
    print(outfile)

    #write.csv(pft_sex_race_df, file = outfile, row.names=FALSE)
  }
  }
}

```

```{r}
calculate_mse <- function(model) {
    mse <- mean(residuals(model)^2)
  return(mse)
}

calculate_RSS <- function(model) {
    mse <- sum(residuals(model)^2)
  return(mse)
}

cohorts <- c('nhanes', 'ukb')

for (i in seq_along(cohorts)){
  if (cohorts[i] == 'nhanes') {
      file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-07-18.csv", cohorts[i], cohorts[i])
  }
  else {
      file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-12-04.csv", cohorts[i], cohorts[i])
  }
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% mutate(FEV1 = FEV1*1000, FVC = FVC*1000)

  for (pft in c("FEV1", "FVC")){
    pft_race_df = data.frame()
    # Initial model with only the intercept
    if (cohorts[i] == 'nhanes3') {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Other + RIDRETH_White"
                                , pft)
        }
        else if (cohorts[i] == 'nhanes4') {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White"
                                , pft)
        }
        else if (cohorts[i] == 'nhanes') {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White"
                                , pft)
        }
        else {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White"
                                , pft)
        } 
    const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR ", pft)
    formula <- as.formula(const_pred)
    model <- lm(formula, data = cohort_df)
    
    # Calculate the initial RSS
    ESTIMATE <- calculate_estimate(model)
    
    df = as.data.frame(summary(model)$coefficients)
    df <- cbind(Covariate = rownames(df), df)
    df['Model'] <- 'Baseline'
    df['Category'] <- 'Baseline'
    df['R2'] <- summary(model)$r.squared
    df['MSE'] <- ESTIMATE
    df['RSS'] <- calculate_RSS(model) 
    pft_race_df <- rbind(pft_race_df,df)
          
    categories <- c('Anthropometrics' , 'Exposures', 'Sociodemographics', 'PRS') 
    for (cat in categories) {
          pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
                          mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
          pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(cohort_df))
          predictors <- as.character(pred_map[pred_map$Category==cat, ]$Predictor.Name.R)
          predictors <- predictors[!grepl("_norm(al)?$|composite?$", predictors)]

          if (length(predictors)==0) {next}
          # Create an empty list to store the selected variables
          selected_variables <- c()
          
          # Number of remaining variables
          remaining_variables <- setdiff(predictors, selected_variables)
          
          best_EST <- Inf
          while (length(remaining_variables) > 0) {
            # Initialize variables for tracking the best addition
            best_variable <- ""
            
            # Loop through the remaining variables and find the one that reduces the RSS the most
            for (variable in remaining_variables) {
              new_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "), "+", variable)))
              new_EST <- calculate_estimate(new_model)
                
              if (new_EST <= best_EST) {
                best_variable <- variable
                best_EST <- new_EST
              }
            }
            
          # Add the best variable to the selected variables
          selected_variables <- c(selected_variables, best_variable)
          
          # Update the model and RSS with the added variable
          model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
          ESTIMATE <- best_EST
      
          # Remove the added variable from the remaining variables
          remaining_variables <- setdiff(remaining_variables, best_variable)
          df = as.data.frame(summary(model)$coefficients)
          df <- cbind(Covariate = rownames(df), df)
          df['Model'] <- sprintf("+ %s", best_variable)
          df['Category'] <- cat
          df['MSE'] <- ESTIMATE
          df['R2'] <- summary(model)$r.squared
          df['RSS'] <- calculate_RSS(model) 

          pft_race_df <- rbind(pft_race_df,df)
          row.names(pft_race_df) <- NULL
          }
    }
    
  print(cohorts[i])
  print(pft)
  outfile = sprintf('data/%s/%s_race_est_c_%s_wo_race_2023-12-15.csv', cohorts[i], cohorts[i], tolower(pft))
  print(outfile)
  print(pft_race_df)
  write.csv(pft_race_df, file = outfile, row.names=FALSE)
  }
}

```

```{r}
paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))
```

```{r}
calculate_estimate <- function(model) {
    mse <- mean(residuals(model)^2)
  return(mse)
}

for (i in seq_along(cohorts)){
  if (cohorts[i] == 'nhanes') {
      file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-07-18.csv", cohorts[i], cohorts[i])
  }
  else {
      file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-12-04.csv", cohorts[i], cohorts[i])
  }
  
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% mutate(FEV1 = FEV1*1000, FVC = FVC*1000)

  for (pft in c("FEV1", "FVC")){
    pft_race_df = data.frame()
    # Initial model with only the intercept
    if (cohorts[i] == 'nhanes3') {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Other + RIDRETH_White", pft)
        }
        else if (cohorts[i] == 'nhanes4') {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
        }
    else if (cohorts[i] == 'nhanes') {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Hispanic + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
    }
    
        else {
          const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR + RIDRETH_Black + RIDRETH_Asian + RIDRETH_Other + RIDRETH_White", pft)
        }
    
    #const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR*RIDAGEYR", pft)
    formula <- as.formula(const_pred)
    model <- lm(formula, data = cohort_df)
    
    # Calculate the initial RSS
    ESTIMATE <- calculate_estimate(model)
    
    df = as.data.frame(summary(model)$coefficients)
    df <- cbind(Covariate = rownames(df), df)
    df['Model'] <- 'Baseline'
    df['Category'] <- 'Baseline'
    df['R2'] <- summary(model)$r.squared
    df['MSE'] <- ESTIMATE
    df['RSS'] <- calculate_RSS(model) 

    pft_race_df <- rbind(pft_race_df,df)
          
    categories <- c('Anthropometrics' , 'Exposures', 'Sociodemographics', 'PRS') #, #'Ancestry', #'Anthropometrics at 10', 'Genetics')
    
    #for (cat in categories) {
    pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE,  sep = ",") %>%
                    mutate(Predictor.Name.R = gsub(" ", ".", Predictor.Name)) 
    pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(cohort_df))
    pred_map <- subset(pred_map, Category %in% categories)

    predictors <- as.character(pred_map$Predictor.Name.R)
          #predictors <- as.character(pred_map[pred_map$Category==cat, ]$Predictor.Name.R)
    predictors <- sort(predictors[!grepl("BMXARML|(al)?$|composite?$", predictors)])

    # Create an empty list to store the selected variables
    selected_variables <- c()
    
    # Number of remaining variables
    remaining_variables <- setdiff(predictors, selected_variables)
    
    best_EST <- Inf
    while (length(remaining_variables) > 0) {
      # Initialize variables for tracking the best addition
      best_variable <- ""
      
      # Loop through the remaining variables and find the one that reduces the RSS the most
      for (variable in remaining_variables) {
        new_model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "), "+", variable)))
        new_EST <- calculate_estimate(new_model)
          
        #if (new_EST <= best_EST) {
        #  best_variable <- variable
        #  best_EST <- new_EST
        best_variable <- variable
        best_EST <- new_EST
        #}
      }
    print(remaining_variables)
    # Add the best variable to the selected variables
    selected_variables <- c(selected_variables, best_variable)
    
    # Update the model and RSS with the added variable
    model <- update(model, formula = as.formula(paste(sprintf("%s +", const_pred), paste(selected_variables, collapse = " + "))))
    ESTIMATE <- best_EST

    # Remove the added variable from the remaining variables
    remaining_variables <- setdiff(remaining_variables, best_variable)
    df = as.data.frame(summary(model)$coefficients)
    df <- cbind(Covariate = rownames(df), df)
    df['Model'] <- sprintf("+ %s", best_variable)
    df['Category'] <- "All"
    df['R2'] <- summary(model)$r.squared
    df['MSE'] <- ESTIMATE
    df['RSS'] <- calculate_RSS(model) 
    
    pft_race_df <- rbind(pft_race_df,df)
    row.names(pft_race_df) <- NULL
    }
  print(cohorts[i])
  print(pft)
  print(pft_race_df)
  outfile = sprintf('data/%s/%s_race_est_nc_sorted_w_race_%s_2023-12-15.csv', cohorts[i], cohorts[i], tolower(pft))
  print(outfile)
  print(pft_race_df)
  write.csv(pft_race_df, file = outfile, row.names=FALSE)
  }
}
```

```{r}
calculate_estimate <- function(model) {
    mse <- mean(residuals(model)^2)
  return(mse)
}
cohorts <- c('ukb')

for (i in seq_along(cohorts)){
  if (cohorts[i] == 'nhanes') {
      file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-07-18.csv", cohorts[i], cohorts[i])
  }
  else {
      file_path <- sprintf("../../processed/%s/%s_healthy_encoded_2023-12-04.csv", cohorts[i], cohorts[i])
  }
  
  cohort_df <- read.table(file_path, row.names = NULL, header = TRUE, sep = ",")
  cohort_df <- cohort_df %>% mutate(FEV1 = FEV1*1000, FVC = FVC*1000)
  
  print('here')
  white_cohort <- cohort_df[cohort_df$RIDRETH_White == 1, ]
  black_cohort <- cohort_df[cohort_df$RIDRETH_Black == 1, ]
  asian_cohort <- cohort_df[cohort_df$RIDRETH_Asian == 1, ]
  #hispanic_cohort <- cohort_df[cohort_df$RIDRETH_Hispanic == 1, ]
  other_cohort <- cohort_df[cohort_df$RIDRETH_Other == 1, ]
  # Add other cohorts as needed

  race_cohorts <- list(white = white_cohort, black = black_cohort, asian = asian_cohort, other = other_cohort) #, hispanic = hispanic_cohort )
  pft_race_df <- data.frame()  # Initialize an empty data frame
  
  for (race in names(race_cohorts)) {
    cohort <- race_cohorts[[race]]
    
    const_pred <- sprintf("%s ~ RIAGENDR + RIDAGEYR + RIAGENDR * RIDAGEYR", pft)
    formula <- as.formula(const_pred)
    
    model <- lm(formula, data = cohort)
    
    # Calculate the initial RSS
    ESTIMATE <- calculate_estimate(model)
    
    df = as.data.frame(summary(model)$coefficients)
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
      pred_map <- read.table("../../data/predictor_category_map.csv", row.names = NULL, header = TRUE, sep = ",") %>% mutate(Predictor.Name.R = gsub(" ", ".",Predictor.Name))
      pred_map <- subset(pred_map, Predictor.Name.R %in% colnames(cohort_df))
      pred_map <- subset(pred_map, Category %in% cat)  # Adjust this line
      
      predictors <- as.character(pred_map$Predictor.Name.R)
      predictors <- sort(predictors[!grepl("_norm(al)?$|composite?$", predictors)])
      
      # Initialize a data frame to store results for the current category
      category_df <- data.frame()
      
      # Loop through the predictors for the current category
      for (variable in predictors) {
        selected_variables <- c(selected_variables, variable)

        # Update the model and RSS with the added variable
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

        # Print statements moved inside the loop to see results for each variable
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
```
