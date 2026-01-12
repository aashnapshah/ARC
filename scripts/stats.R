# descriptive_tables.R
# Author: Aashna Shah
# Date: 2023-06-09

library(survey)
library(dplyr)
library(gtsummary)
library(knitr)
library(gt)

# Load mapping dictionary
df <- read.table("../data/column_description_dictionary.csv", header = TRUE, sep = ",")
mapping_dict <- df %>% dplyr::pull(Description) %>% setNames(df$X)
additional_dict <- c(
  FEV1 = 'FEV1',
  FVC = 'FVC',
  `FEV1.FVC` = 'FEV1/FVC',
  BMXWT10 = 'Weight at 10',
  BMXHT10 = 'Height at 10',
  DMD_INC = 'Household Income',
  DMD_BORN_SubRegion = 'Subregion Born In',
  DMD_PM2.5 = 'Exposure to PM 2.5',
  DMD_PM2.5_10 = 'Exposure to PM 10'
)
mapping_dict <- c(mapping_dict, additional_dict)

### NHANES (Main) ###
nhanes_pft <- read.table("../processed/nhanes/nhanes_healthy_2023-07-18.csv", row.names = NULL, header = TRUE, sep = ",")
nhanes_pft <- nhanes_pft %>%
  mutate(SEX = ifelse(RIAGENDR == 'True', "Female", "Male"),
         DMD_HH_Size = as.numeric(DMD_HH_Size),
         RIDRETH = factor(RIDRETH, levels = c('White', 'Black', 'Hispanic', 'Asian', 'Other')))
logical_cols <- sapply(nhanes_pft, \(x) all(x %in% c('True', 'False')))
nhanes_pft[logical_cols] <- lapply(nhanes_pft[logical_cols], as.logical)

nhanes_columns <- c("RIAGENDR", "RIDAGEYR", "FEV1", "FVC", "FEV1.FVC", 
  "BMXHT", "BMXSIT", "BMXBMI", # "BMXARMC", 
  # "BMXARML", "BMXLEG", "BMXWAIST", 
  "DMD_IncPovRatio", 
  "DMD_HH_Size", "DMD_Born_US", "DMD_Finish_HS", 
  # "DMD_Married",
  "DMD_Health_Ins", # "DMD_Priv_Ins", 
  "DMD_English_Lang", "DMD_Home_Smoke"
  #, "DMD_Veteran", "DMD_Military_HC"
)
nhanes_table <- nhanes_pft %>%
  tbl_summary(
    by = RIDRETH, 
    type = list(DMD_HH_Size ~ "continuous"),
    include = all_of(nhanes_columns),
    statistic = list(
      all_continuous() ~ "{mean} [{sd}]",
      all_categorical() ~ "{p}%"
    ),
    digits = all_continuous() ~ 1
  )
nhanes_table$table_body <- nhanes_table$table_body %>%
  mutate(label = ifelse(label %in% names(mapping_dict), mapping_dict[label], label))

print(nhanes_table)

# P-values by race comparisons
comparison_races <- c('Black', 'Asian', 'Other', 'Hispanic')
get_p_value_race <- function(data, variable, reference_group, comparison_group) {
  ref_group_data <- data[data$RIDRETH == reference_group, ][[variable]]
  comp_group_data <- data[data$RIDRETH == comparison_group, ][[variable]]
  t_test_res <- tryCatch(
    t.test(ref_group_data, comp_group_data, var.equal = FALSE),
    error = function(e) NA
  )
  if (is.list(t_test_res)) p_value <- t_test_res$p.value else p_value <- NA
  return(p_value)
}
p_values_list <- list()
for (col in nhanes_columns) {
  p_values_race <- sapply(comparison_races, function(race) {
    get_p_value_race(nhanes_pft, col, 'White', race)
  })
  p_values_list[[col]] <- p_values_race
}
p_values_df <- as.data.frame(p_values_list)
rownames(p_values_df) <- comparison_races
print(p_values_df)

### NHANES III ###
nhanes3_pft <- read.table("../processed/nhanes3/nhanes3_healthy_2023-07-18.csv", row.names = NULL, header = TRUE, sep = ",")
nhanes3_pft <- nhanes3_pft %>%
  mutate(SEX = ifelse(RIAGENDR == 'True', "Female", "Male"),
         DMD_HH_Size = as.numeric(DMD_HH_Size),
         RIDRETH = factor(RIDRETH, levels = c('White', 'Black', 'Hispanic', 'Other')))
logical_cols <- sapply(nhanes3_pft, \(x) all(x %in% c('True', 'False')))
nhanes3_pft[logical_cols] <- lapply(nhanes3_pft[logical_cols], as.logical)

nhanes3_pft_survey <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = nhanes3_pft
)
nhanes3_columns <- c("RIDAGEYR", "FEV1", "FVC", "FEV1.FVC", 
  "BMXHT", "BMXSIT", "BMXBMI", "BMXARMC", 
  "BMXARML", "BMXLEG", "BMXWAIST", "DMD_IncPovRatio", 
  "DMD_HH_Size", "DMD_Born_US", "DMD_Finish_HS", 
  "DMD_Married", "DMD_Health_Ins", "DMD_Priv_Ins", 
  "DMD_English_Lang", "DMD_Home_Smoke", "DMD_Veteran", "DMD_Military_HC"
)
nhanes3_table <- nhanes3_pft %>% tbl_strata(
  strata = RIDRETH,
  ~.x %>% tbl_summary(
    by = SEX,
    type = list(DMD_HH_Size ~ "continuous"),
    include = all_of(nhanes3_columns),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{p}%"
    ),
    digits = all_continuous() ~ 1
  )
)
nhanes3_table$table_body <- nhanes3_table$table_body %>%
  mutate(label = ifelse(label %in% names(mapping_dict), mapping_dict[label], label))
nhanes3_table_survey <- nhanes3_pft_survey %>% tbl_strata(
  strata = RIDRETH,
  ~.x %>% tbl_svysummary(
    by = SEX,
    type = list(DMD_HH_Size ~ "continuous"),
    include = all_of(nhanes3_columns),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{p}%"
    ),
    digits = all_continuous() ~ 1
  )
)
nhanes3_table_survey$table_body <- nhanes3_table_survey$table_body %>%
  mutate(label = ifelse(label %in% names(mapping_dict), mapping_dict[label], label))

print(nhanes3_table_survey)
print(nhanes3_table)

# Save GT tables as images
nhanes3_table %>%
  bold_labels() %>% as_gt() %>%
  gt::tab_header("Table 1: NHANES III Summary Statistics") %>%
  gt::tab_options(
    table.font.size = "small",
    table.font.style = "sans-serif",
    data_row.padding = gt::px(1)
  ) %>%
  gt::gtsave(file = "tables/nhanes3_descriptive_stats.png")
nhanes3_table_survey %>%
  bold_labels() %>% as_gt() %>%
  gt::tab_header("Table 2: NHANES III Survey Summary Statistics") %>%
  gt::tab_options(
    table.font.size = "small",
    table.font.style = "sans-serif",
    data_row.padding = gt::px(1)
  ) %>%
  gt::gtsave(file = "tables/nhanes3_survey_descriptive_stats.png")

### NHANES IV ###
nhanes4_pft <- read.table("../processed/nhanes4/nhanes4_healthy_2023-07-18.csv", row.names = NULL, header = TRUE, sep = ",")
nhanes4_pft <- nhanes4_pft %>%
  mutate(
    DMD_HH_Size = as.numeric(DMD_HH_Size),
    SEX = ifelse(RIAGENDR == 'True', "Female", "Male"),
    RIDRETH = factor(RIDRETH, levels = c('White', 'Black', 'Hispanic', 'Asian', 'Other'))
  )
logical_cols <- sapply(nhanes4_pft, \(x) all(x %in% c('True', 'False')))
nhanes4_pft[logical_cols] <- lapply(nhanes4_pft[logical_cols], as.logical)

nhanes4_pft_survey <- svydesign(
  id = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMEC6YR,
  nest = TRUE,
  data = nhanes4_pft
)

nhanes4_table <- nhanes4_pft %>% tbl_strata(
  strata = RIDRETH,
  ~.x %>% tbl_summary(
    by = SEX,
    type = list(DMD_HH_Size ~ "continuous"),
    include = nhanes_columns,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{p}%"
    ),
    digits = all_continuous() ~ 1
  )
)
nhanes4_table$table_body <- nhanes4_table$table_body %>%
  mutate(label = ifelse(label %in% names(mapping_dict), mapping_dict[label], label))

nhanes4_table_survey <- nhanes4_pft_survey %>% tbl_strata(
  strata = RIDRETH,
  ~.x %>% tbl_svysummary(
    by = SEX,
    type = list(DMD_HH_Size ~ "continuous"),
    include = nhanes_columns,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{p}%"
    ),
    digits = all_continuous() ~ 1
  )
)
nhanes4_table_survey$table_body <- nhanes4_table_survey$table_body %>%
  mutate(label = ifelse(label %in% names(mapping_dict), mapping_dict[label], label))

# Save GT tables for NHANES IV
nhanes4_table %>%
  bold_labels() %>% as_gt() %>%
  gt::tab_header("Table 2: NHANES IV Summary Statistics") %>%
  gt::tab_options(
    table.font.size = "small",
    table.font.style = "sans-serif",
    data_row.padding = gt::px(1)
  ) %>%
  gt::gtsave(file = "tables/nhanes4_descriptive_stats.png")
nhanes4_table_survey %>%
  bold_labels() %>% as_gt() %>%
  gt::tab_header("Table 4: NHANES IV Survey Summary Statistics") %>%
  gt::tab_options(
    table.font.size = "small",
    table.font.style = "sans-serif",
    data_row.padding = gt::px(1)
  ) %>%
  gt::gtsave(file = "tables/nhanes4_survey_descriptive_stats.png")

### UK Biobank ###
ukb_pft <- read.table("../processed/ukb/ukb_healthy_2023-12-04.csv", row.names = NULL, header = TRUE, sep = ",") %>%
  select(starts_with('RI') | starts_with('F') | starts_with('BMX') | starts_with('DMD') | starts_with('ratio') | starts_with('fev1_score') | starts_with('composite'))
logical_cols <- sapply(ukb_pft, \(x) all(x %in% c('True', 'False')))
ukb_pft[logical_cols] <- lapply(ukb_pft[logical_cols], as.logical)
ukb_pft <- ukb_pft %>%
  mutate(
    SEX = ifelse(RIAGENDR == TRUE, "Female", "Male"),
    RIDRETH = factor(RIDRETH, levels = c('White', 'Black', 'Asian', 'Other')),
    Income = as.numeric(DMD_INC)
  )
# Simple summary table, a subset of columns
ukb_table_simple <- ukb_pft %>%
  select(RIDRETH, RIDAGEYR, SEX, FEV1, FVC, `FEV1.FVC`) %>%
  tbl_strata(
    strata = RIDRETH,
    ~.x %>% tbl_summary(
      by = SEX,
      type = list(where(is.numeric) ~ "continuous"),
      statistic = list(all_continuous() ~ "{mean} ({sd})"),
      digits = all_continuous() ~ 1
    )
  )
print(ukb_table_simple)

# More detailed UKB Table with p-values
cols <- c("RIDAGEYR", "RIAGENDR", "FEV1", "FVC", "FEV1.FVC", 
          "BMXHT", "BMXSIT", "BMXBMI", "BMXWAIST",
          # "BMXHIP", "BMXWT10", "BMXHT10", 
          "DMD_Finish_HS", # "DMD_PM10"
          "DMD_Born_UK", "DMD_Work_Smoke", "DMD_Home_Smoke",
          # "DMD_Work_Breathing_Probs",
          "ratio_score", "ratio_score_norm", "fev1_score",
          "fev1_score_norm", "composite_normal", "composite",
          "Income")

ukb_table <- ukb_pft %>%
  tbl_summary(
    by = RIDRETH,
    percent = "column",
    type = list(Income ~ "continuous"),
    include = cols,
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{p}%"
    ),
    digits = all_continuous() ~ 1
  )

# List of races for comparison in UKB data
comparison_races_ukb <- c('Black', 'Asian', 'Other')
get_p_value_race <- function(data, variable, reference_group, comparison_group) {
  ref_group_data <- data[data$RIDRETH == reference_group, ][[variable]]
  comp_group_data <- data[data$RIDRETH == comparison_group, ][[variable]]
  t_test_res <- tryCatch(
    t.test(ref_group_data, comp_group_data, var.equal = FALSE),
    error = function(e) NA
  )
  if (is.list(t_test_res)) p_value <- t_test_res$p.value else p_value <- NA
  return(p_value)
}
p_values_list <- list()
for (col in cols) {
  p_values_race <- sapply(comparison_races_ukb, function(race) {
    get_p_value_race(ukb_pft, col, 'White', race)
  })
  p_values_list[[col]] <- p_values_race
}
p_values_df <- as.data.frame(p_values_list)
rownames(p_values_df) <- comparison_races_ukb
print(p_values_df)

# Save UK Biobank table as GT
ukb_table %>%
  bold_labels() %>% as_gt() %>%
  gt::tab_header("Table 3: UKBioBank Summary Statistics") %>%
  gt::tab_options(
    table.font.size = "small",
    table.font.style = "sans-serif",
    data_row.padding = gt::px(1)
  ) %>%
  gt::gtsave(file = "tables/ukb_descriptive_stats.png")

# End of script
