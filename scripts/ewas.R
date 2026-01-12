# DWAS.R
# Author: Aashna Shah
# Date: 2023-04-11

# ---- Package Setup ----
packages <- c(
  "binom","caret","dplyr","fst","gamlss","gamlss.dist",
  "ggplot2","GGally","ggbeeswarm","ggridges","ggsci","ggpubr", 
  "glmnet","gridExtra","kableExtra","knitr","lubridate",
  "MASS","plotly","plyr","pROC","readr","RNHANES","rspiro",
  "scales","srvyr","stringr","survey","splines", "tidyr"
)

package.check <- lapply(packages, function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, dependencies = TRUE)
  require(x, character.only = TRUE)
})

# ---- Data Load and Preprocessing ----

nhanes_pft <- read.table("../../processed/nhanes/nhanes_healthy_scaled_2023-07-18.csv", row.names = NULL, header = TRUE, sep = ",")
nhanes_pft$RIDRETH <- factor(ifelse(nhanes_pft$RIDRETH_Other == 1, '_Other',
                             ifelse(nhanes_pft$RIDRETH_Black == 1, '_Black',
                                    ifelse(nhanes_pft$RIDRETH_Hispanic == 1, '_Hispanic',
                                    ifelse(nhanes_pft$RIDRETH_Asian == 1, '_Asian',
                                    ifelse(nhanes_pft$RIDRETH_White == 1, '_White', NA)))))) 

nhanes_pft <- nhanes_pft %>% mutate(
  FEV1 = FEV1 * 1000,
  FVC = FVC * 1000,
  FEV1.FVC = FEV1.FVC * 1000
) %>%
  select(-starts_with("RIDRETH_"), -starts_with("DMD_INC_"))

nhanes_pft$RIDRETH <- relevel(nhanes_pft$RIDRETH, ref = '_White')

nhanes3_pft <- read.table("../../processed/nhanes3/nhanes3_healthy_scaled_2023-07-18.csv", row.names = NULL, header = TRUE, sep = ",")
nhanes3_pft$RIDRETH <- factor(ifelse(nhanes3_pft$RIDRETH_Other == 1, '_Other',
                             ifelse(nhanes3_pft$RIDRETH_Black == 1, '_Black',
                                    ifelse(nhanes3_pft$RIDRETH_Hispanic == 1, '_Hispanic',
                                    ifelse(nhanes3_pft$RIDRETH_White == 1, '_White', NA)))))
nhanes3_pft <- nhanes3_pft %>% mutate(
  FEV1 = FEV1 * 1000,
  FVC = FVC * 1000,
  FEV1.FVC = FEV1.FVC * 1000
) %>%
  select(-starts_with("RIDRETH_"), -starts_with("DMD_INC_"))
nhanes3_pft$RIDRETH <- relevel(nhanes3_pft$RIDRETH, ref = '_White')

nhanes4_pft <- read.table("../../processed/nhanes4/nhanes4_healthy_scaled_2023-07-18.csv", row.names = NULL, header = TRUE, sep = ",")
nhanes4_pft$RIDRETH <- factor(ifelse(nhanes4_pft$RIDRETH_Other == 1, '_Other',
                             ifelse(nhanes4_pft$RIDRETH_Black == 1, '_Black',
                                    ifelse(nhanes4_pft$RIDRETH_Hispanic == 1, '_Hispanic',
                                    ifelse(nhanes4_pft$RIDRETH_Asian == 1, '_Asian',
                                    ifelse(nhanes4_pft$RIDRETH_White == 1, '_White', NA))))))
nhanes4_pft <- nhanes4_pft %>% mutate(
  FEV1 = FEV1 * 1000,
  FVC = FVC * 1000,
  FEV1.FVC = FEV1.FVC * 1000
) %>%
  select(-starts_with("RIDRETH_"), -starts_with("DMD_INC_"))
nhanes4_pft$RIDRETH <- relevel(nhanes4_pft$RIDRETH, ref = '_White')

ukb_pft <- read.table("../../processed/ukb/ukb_healthy_scaled_2023-12-04.csv", row.names = NULL, header = TRUE, sep = ",")
ukb_pft$RIDRETH <- factor(ifelse(ukb_pft$RIDRETH_White == 1, '_White',
                                 ifelse(ukb_pft$RIDRETH_Black == 1, '_Black',
                                        ifelse(ukb_pft$RIDRETH_Asian == 1, '_Asian',
                                               ifelse(ukb_pft$RIDRETH_Other == 1, '_Other', NA)))))
ukb_pft <- ukb_pft %>% mutate(
  FEV1 = FEV1 * 1000,
  FVC = FVC * 1000,
  FEV1.FVC = FEV1.FVC * 1000
) %>%
  select(-starts_with("RIDRETH_"), -starts_with("DMD_INC_"))
ukb_pft$RIDRETH <- relevel(ukb_pft$RIDRETH, ref = '_White')
# ukb_pft$DMD_INC <- relevel(ukb_pft$DMD_INC, ref = '_Four.or.more')

# ---- EWAS loop for NHANES and UKB ----

cohorts <- c('nhanes', 'ukb')
cohorts_df <- list(nhanes_pft, ukb_pft)

for (i in seq_along(cohorts)) {
  cohort_df <- cohorts_df[[i]]

  for (pft in c("FEV1", "FVC", "FEV1.FVC")) {
    pft_df <- data.frame()
    for (col_name in colnames(cohort_df)) {
      if (col_name %in% c('RIDAGEYR', 'RIAGENDR', 'X', 
                          'FVC', 'FEV1', 'FEV1.FVC', 
                          'SDMVPSU', 'SDMVSTRA', 'WTMEC6YR')) {
        next
      }
      formula <- as.formula(sprintf('%s ~ RIDAGEYR + RIAGENDR + %s', pft, col_name))
      model <- lm(formula, data = cohort_df)
      df <- as.data.frame(summary(model)$coefficients)
      df <- cbind(Covariate = rownames(df), df)
      df['Model'] <- col_name

      if (nrow(df) < 4) next
      pft_df <- rbind(pft_df, df)
      row.names(pft_df) <- NULL
    }
    pft_name <- gsub("\\.", "", pft)
    output_file <- sprintf("data/%s/%s_ewas_results_%s_refW_2023-12-04.csv", cohorts[i], cohorts[i], pft_name)
    print(pft_df)
    write.csv(pft_df, file = output_file, row.names = FALSE)
  }
}

# ---- EWAS loop for NHANES3/4 (survey weights) ----

cohorts <- c('nhanes3', 'nhanes4')
for (i in seq_along(cohorts)) {
  file_name <- sprintf("../../processed/%s_healthy_scaled_2023-06-08.csv", cohorts[i])
  cohort_df <- read.table(file_name, row.names = NULL, header = TRUE, sep = ",")

  cohort_df_survey <- svydesign(
    id = ~SDMVPSU,
    strata = ~SDMVSTRA,
    weights = ~WTMEC6YR,
    nest = TRUE,
    data = cohort_df
  )

  for (pft in c("FEV1", "FVC", "FEV1.FVC")) {
    pft_df <- data.frame()
    for (col_name in colnames(cohort_df)) {
      if (col_name %in% c('RIDAGEYR', 'RIAGENDR', 'X', 'FVC', 'FEV1', 'FEV1.FVC',
                          'SDMVPSU', 'SDMVSTRA', 'WTMEC6YR')) {
        next
      }
      formula <- as.formula(sprintf('%s ~ cs((RIDAGEYR), df = 7) + RIAGENDR + cs((RIDAGEYR), df = 7):RIAGENDR + %s', pft, col_name))
      model <- svyglm(formula, design = cohort_df_survey)

      df <- as.data.frame(summary(model)$coefficients)
      df <- cbind(Covariate = rownames(df), df)
      df['Model'] <- col_name

      if (nrow(df) < 4) next

      pft_df <- rbind(pft_df, df)
      row.names(pft_df) <- NULL
    }
    pft_name <- gsub("\\.", "", pft)
    output_file <- sprintf("data/%s/%s_ewas_survey_results_%s_2023-06-09.csv", cohorts[i], cohorts[i], pft_name)
    write.csv(pft_df, file = output_file, row.names = FALSE)
  }
}
