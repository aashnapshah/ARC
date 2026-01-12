# Load dependencies
source('/n/groups/patel/aashna/PFT-ML-2024/helpers/params.R')
source('/n/groups/patel/aashna/PFT-ML-2024/helpers/functions.R')

# Utility to replace NAs with a default (FALSE by default)
na_to_F <- function(vec, alt = FALSE) replace(vec, is.na(vec), alt)

# FIELD CODEBOOK - maps readable variable names to UKB field IDs
codebook <- c(
  ID              = "eid",    Sex              = "31",     Age           = "21003",
  WaistCircum     = "48",     HipCircum        = "49",     Height        = "50",
  Sit_Box         = "51",     Assessment_Date  = "53",     Assessment_Center = "54",
  CurrentSmoker   = "1239",   HH_Smoking_Exp   = "1259",   HH_Smoking_Exp_Hr_Per_Wk = "1269",
  Outside_Smoking_Exp_Hr_Per_Wk = "1279",
  Size_Age10      = "1687",   Height_Age10     = "1697",   Maternal_Smoking = "1787",
  Wheezing_Whistling_Past_Yr  = "2316",
  Smoked_100      = "2644",   Spiro_Quality_Rank = "3059",
  Spiro_Accept_1  = "3061",   Spiro_Accept_2   = "3061.0.1", Spiro_Accept_3 = "3061.0.2",
  FVC_1           = "3062",   FVC_2           = "3062.0.1", FVC_3 = "3062.0.2",
  FEV1_1          = "3063",   FEV1_2          = "3063.0.1", FEV1_3 = "3063.0.2",
  BoxHeight       = "3077",   Yr_Arrived_UK    = "3659",   Education = "6138",
  SR_Cancer       = "20001",  SR_NonCancer     = "20002",  SR_Cancer_Age = "20009",
  Sitting_Height  = "20015",  Birth_Weight     = "20022",  Country_Birth = "20115",
  Spiro_Reprod    = "20152",  Ethnicity        = "21000",  BMI          = "21001",
  Weight          = "21002",  Cough_Most_Days  = "22502",  Phlegm_Most_Days = "22504",
  Genetic_Ancestry = "22006", Best_FEV1        = "20150",  Best_FVC     = "20151",
  Work_Fumes      = "22610",  Work_Smoking_Exp = "22611",  Work_Diesel  = "22615",
  Work_Breathing_Probs = "22616", Work_Stop_Probs_Improved = "22618",
  PM10 = "24005", PM2.5 = "24006", PM2.5_10 = "24008"
)
codebook <- setNames(names(codebook), codebook)

# Read UKB data
ukb_raw <- readRDS("data/joined_ukb.rds")

# Make columns readable
ukb_cols <- gsub('.0.0', '', gsub('f.', '', colnames(ukb_raw), fixed = TRUE), fixed = TRUE)
coded <- ukb_cols %in% names(codebook)
colnames(ukb_raw)[coded] <- codebook[ukb_cols[coded]]

# Format date columns
ukb_raw$Assessment_Date <- as.Date(ukb_raw$Assessment_Date, "%Y-%m-%d")

for (col in grep("_days", colnames(ukb_raw), value = TRUE)) {
  ukb_raw[[col]] <- ifelse(ukb_raw[[col]] == "0370-01-01", NA, ukb_raw[[col]])
  ukb_raw[[col]] <- as.Date(ukb_raw[[col]], "%Y-%m-%d")
  diff_col <- sub("_days", "_diff", col)
  ukb_raw[[diff_col]] <- as.numeric(ukb_raw[[col]] - ukb_raw[["Assessment_Date"]])
}

# Derive fields and flags
ukb_processed <- ukb_raw %>%
  dplyr::mutate(
    ANCESTRY = pop,
    RIDRETH = dplyr::case_when(
      Ethnicity %in% c("British","Irish","Any other white background") ~ "White",
      grepl("White and|mixed|Mixed", Ethnicity)                        ~ "Multiracial",
      Ethnicity %in% c("Indian","Pakistani","Bangladeshi","Asian or Asian British","Any other Asian background") ~ "South Asian",
      Ethnicity %in% c("Caribbean","African","Black or Black British","Any other Black background") ~ "Black",
      grepl("Other ethnic group|Do not know|Prefer not to answer", Ethnicity) ~ "Other",
      TRUE ~ Ethnicity
    ),
    IS_FEMALE = Sex == "Female",
    RIDAGEYR  = Age,
    BMXHT     = Height,
    BMXWT     = Weight,
    BMXBMI    = BMI,
    BMXWAIST  = WaistCircum,
    BMXHIP    = HipCircum,
    BMXSIT    = Sitting_Height,
    BMXWT10   = Size_Age10,
    BMXHT10   = Height_Age10,

    # Spirometry acceptability and best value extracts
    Accept_1 = na_to_F(Spiro_Accept_1 == 0 | Spiro_Accept_1 == 32),
    Accept_2 = na_to_F(Spiro_Accept_2 == 0 | Spiro_Accept_2 == 32),
    Accept_3 = na_to_F(Spiro_Accept_3 == 0 | Spiro_Accept_3 == 32),

    FEV1_1a = ifelse(Accept_1 & FEV1_1 < 10, FEV1_1, NA),
    FEV1_2a = ifelse(Accept_2 & FEV1_2 < 10, FEV1_2, NA),
    FEV1_3a = ifelse(Accept_3 & FEV1_3 < 10, FEV1_3, NA),
    FVC_1a  = ifelse(Accept_1 & FVC_1 < 10, FVC_1, NA),
    FVC_2a  = ifelse(Accept_2 & FVC_2 < 10, FVC_2, NA),
    FVC_3a  = ifelse(Accept_3 & FVC_3 < 10, FVC_3, NA),

    FEV1 = pmax(FEV1_1a, FEV1_2a, FEV1_3a, na.rm = TRUE),
    FVC  = pmax(FVC_1a, FVC_2a, FVC_3a, na.rm = TRUE),

    FEV1_diff = apply(data.frame(FEV1, FEV1_1, FEV1_2, FEV1_3), 1, function(row) sort(abs(row[1] - row[-1]))[2]),
    FVC_diff  = apply(data.frame(FVC,  FVC_1,  FVC_2,  FVC_3),  1, function(row) sort(abs(row[1] - row[-1]))[2]),

    # Repeatability flags
    FEV1_Rep = na_to_F(FEV1_diff <= ifelse(FVC > 1, 0.1501, 0.100)),
    FVC_Rep  = na_to_F(FVC_diff  <= ifelse(FVC > 1, 0.1501, 0.100)),

    FEV1 = ifelse(FEV1_Rep, FEV1, NA),
    FVC  = ifelse(FVC_Rep,  FVC,  NA),

    # Exclusion criteria
    EXCL_SMOKED_100 = na_to_F(Smoked_100 == 1),
    EXCL_SMOKER     = na_to_F(CurrentSmoker %in% c(1, 2)),
    DMD_Mom_Smoked  = na_to_F(Maternal_Smoking == 1),
    DMD_Born_UK     = ifelse(is.finite(Country_Birth), Country_Birth == 354, TRUE),
    DMD_Finish_HS   = na_to_F(Education %in% c(1,2,5,6)),
    DMD_Home_Smoke  = na_to_F(HH_Smoking_Exp),
    DMD_Veteran     = FALSE,
    DMD_PM10        = PM10,
    DMD_PM2.5       = PM2.5,
    DMD_PM2.5_10    = PM2.5_10,
    DMD_Work_Fumes          = Work_Fumes,
    DMD_Work_Smoking_Exp    = Work_Smoking_Exp,
    DMD_Work_Breathing_Probs = Work_Breathing_Probs,
    DMD_Work_Stop_Probs_Improved = Work_Stop_Probs_Improved,
    MCQ_ASTHMA           = na_to_F(asthma_relative == "before"),
    MCQ_ASTHMA_INCID     = na_to_F(asthma_relative == "after"),
    MCQ_BRONCHITIS_EMPH  = na_to_F(chronic_bronch_relative == "before"),
    MCQ_BRONCHITIS_EMPH_INCID = na_to_F(chronic_bronch_relative == "after"),
    MCQ_OTHER_COPD       = na_to_F(otherCOPD_relative == "before"),
    MCQ_OTHER_COPD_INCID = na_to_F(otherCOPD_relative == "after"),
    MCQ_COPD             = MCQ_BRONCHITIS_EMPH | MCQ_OTHER_COPD,
    MCQ_COPD_INCID       = MCQ_BRONCHITIS_EMPH_INCID | MCQ_OTHER_COPD_INCID,
    MCQ_BRONCHIECTASIS   = na_to_F(bronchiectasis_relative == "before"),
    MCQ_BRONCHIECTASIS_INCID = na_to_F(bronchiectasis_relative == "after"),
    MCQ_ASBESTOSIS       = na_to_F(asbestosis_relative == "before"),
    MCQ_ASBESTOSIS_INCID = na_to_F(asbestosis_relative == "after"),
    MCQ_PULM_FIB         = na_to_F(pulm_fibrosis_relative == "before"),
    MCQ_PULM_FIB_INCID   = na_to_F(pulm_fibrosis_relative == "after"),
    MCQ_LUNG_CANCER      = na_to_F(lung_cancer_relative == "before"),
    MCQ_LUNG_CANCER_INCID = na_to_F(lung_cancer_relative == "after"),
    SYMP_COUGH           = na_to_F(Cough_Most_Days == 1),
    SYMP_PHLEGM          = na_to_F(Phlegm_Most_Days == 1),
    SYMP_WHEEZING        = na_to_F(Wheezing_Whistling_Past_Yr == 1)
  )

# High-level inclusion/exclusion filtering and flag assignment
ukb <- ukb_processed %>%
  dplyr::mutate(
    SEQN           = sprintf("U-%s", 1:nrow(ukb_processed)),
    CYCLE_YEAR     = "2006-2010",
    FEV1FVC        = FEV1 / FVC,
    MISSING_AGE    = RIDAGEYR < 3,
    MISSING_DEMO   = is.na(RIDRETH) | is.na(IS_FEMALE),
    MISSING_BMX    = is.na(BMXHT) | is.na(BMXHT) | is.na(BMXBMI) | is.na(BMXWAIST),
    MISSING_SPIRO  = is.na(FEV1) | is.na(FVC),
    MISSING        = MISSING_AGE | MISSING_DEMO | MISSING_BMX | MISSING_SPIRO,
    EXCL_SMOKING   = EXCL_SMOKED_100 | EXCL_SMOKER,
    EXCL_MEDICAL   = MCQ_ASTHMA | MCQ_COPD | MCQ_BRONCHIECTASIS | MCQ_ASBESTOSIS | MCQ_PULM_FIB | MCQ_LUNG_CANCER | FEV1FVC < 0.7,
    EXCL_SYMPTOMS  = SYMP_COUGH | SYMP_WHEEZING | SYMP_PHLEGM,
    EXCL_AGE       = !na_to_F(RIDAGEYR >= 6 & RIDAGEYR <= 79),
    EXCLUDED       = EXCL_AGE | EXCL_SYMPTOMS | EXCL_MEDICAL | EXCL_SMOKING,
    WTMEC6YR       = 1
  ) %>%
  dplyr::filter(!is.na(RIDRETH))

write.table(
  ukb,
  file = sprintf("processed/ukb_yixuan_%s.csv", Sys.Date()),
  quote = FALSE, row.names = FALSE, sep = ","
)

rm(ukb_processed)
