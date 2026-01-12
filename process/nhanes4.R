# Pulmonary Function Reference Values Without Race Stratification in the United States

setwd("/Users/aashnashah/Dropbox/Research/ml_spirometry/")
source('helpers/functions.R')
source('helpers/params.R')

# Years of NHANES to include
cycle_years_input <- 2007:2012
# CDC age bins
cdc_breaks <- c(0, 18, 45, 65, 100)
cdc_labels <- c("0-17","18-44","45-64","65+")
census_breaks <- c(0, 18, 20, 25, 30, 35, 40, 45,
                   50, 55, 60, 65, 70, 75, 80, 85)
census_labels <- c("0-17","18-19","20-24","25-29",
                   "30-34","35-39","40-44",
                   "45-49","50-54","55-59","60-64",
                   "65-69","70-74","75-79","80-84")

keep_var <- list(
  DEMO = c("RIDRETH1","RIDRETH3","RIAGENDR","RIDAGEYR","DMDMARTL",
           "DMDBORN2","DMDBORN4","DMDCITZN","DMDYRSUS","DMDEDUC2",
           "INDHHIN2","INDFMIN2","INDFMPIR","DMDHHSIZ","DMDFMSIZ",
           "DMQMILIT","DMQMILIZ","SIALANG","SIAINTRP",
           "RIDEXPRG","WTINT2YR","WTMEC2YR","SDMVPSU","SDMVSTRA"),
  BMX = c("BMXWT","BMXHT","BMXBMI","BMXLEG","BMXARML",
          "BMXARMC","BMXWAIST","BMDAVSAD"),
  # ECQ = c("ECD010","ECQ020","ECD070A","ECD070B","ECQ080","ECQ090","ECQ095"),
  HIQ = c("HIQ011","HIQ031A","HIQ031B","HIQ031F","HIQ031D","HIQ210"),
  HUQ = c("HUQ010","HUQ030","HUQ050","HUQ060","HUQ071","HUD080"),
  # INQ = c("INQ030","INQ060"),
  MCQ = c("MCQ010","MCQ025","MCQ035","MCQ040","MCQ050","MCQ051",
          "MCQ160B","MCQ160G","MCQ180G","MCQ160K","MCQ170K","MCQ230A",
          "MCQ230B","MCQ230C","MCQ230D","MCQ240N","MCQ300B"),
  OCQ = c("OCD231","OCD241","OCD391","OCD392",
          "OCQ275","OCQ290Q",
          "OCQ510","OCQ520","OCQ530","OCQ540","OCQ550",
          "OCQ560","OCQ570","OCQ580"),
  PAQ = c("PAQ560","PAQ605","PAQ610","PAD615","PAQ620","PAQ625",
          "PAD630","PAQ635","PAQ640","PAD645","PAQ650",
          "PAQ655","PAD660","PAQ665","PAQ670","PAD675","PAD680"),
  RDQ = c("RDQ031","RDQ050","RDQ070","RDQ090","RDQ100",
          "RDQ134", "RDQ135","RDQ140"),
  SMQ = c("SMQ020","SMQ040"),
  SMQFAM = c("SMD410","SMD415","SMD430"),
  SMQRTU = c("SMQ690A","SMQ690B","SMQ690C"),
  SPX = c(
    "SPXNSTAT","SPXNCMT","SPXNFVC","SPXNEV","SPXNFEV5",
    "SPXNFEV7","SPXNFEV1","SPXNFEV3","SPXNFEV6","SPXNPEF",
    "SPXNF257","SPXNFET","SPXNQFVC","SPXNQFV1","SPXNACC","SPXNQEFF",
    "SPDBRONC",
    "SPXBSTAT","SPXBCMT","SPXBFVC","SPXBEV","SPXBFEV5",
    "SPXBFEV7","SPXBFEV1","SPXBFEV3","SPXBFEV6","SPXBPEF",
    "SPXBF257","SPXBFET","SPXBQFVC","SPXBQFV1","SPXBACC","SPXBQEFF"
  )
)

# columns to convert to numeric after loading
num_cols <- c("SEQN","RIDAGEYR","SDMVPSU","SDMVSTRA", keep_var$BMX,
              "SPXNFVC","SPXNEV","SPXNFEV5",
              "SPXNFEV7","SPXNFEV1","SPXNFEV3","SPXNFEV6","SPXNPEF",
              "SPXNF257","SPXNFET",
              "SPXBFVC","SPXBEV","SPXBFEV5",
              "SPXBFEV7","SPXBFEV1","SPXBFEV3","SPXBFEV6","SPXBPEF",
              "SPXBF257","SPXBFET")

# Maps years to letter codes
code_map <- setNames(LETTERS, seq(1999, 2050, 2))
# Converts vector of years into paired cycles
cycle_years <- paste(cycle_years_input[c(T,F)],
                     cycle_years_input[c(F,T)],
                     sep="-") %>%
  setNames(code_map[as.character(cycle_years_input)[c(T,F)]])
# Add letter codes to data file names
data_file_list <- sapply(names(keep_var), function(x)
  sprintf("%s_%s", x, names(cycle_years))) %>% as.vector
# Download and merge NHANES data using `RNHANES` package
raw_data <- nhanes_load_data(data_file_list, year = cycle_years,
                             destination = "data/nhanes4", cache = TRUE,
                             demographics = FALSE, recode = TRUE)


# Kept variables and their descriptions
joined_data <- list()
# Loop through all years and all variables
for (y in 1:length(cycle_years)) {
  for (var in names(keep_var)) {
    letter <- names(cycle_years)[y]
    var_names <- c("SEQN", keep_var[[var]])
    varletter <- sprintf("%s_%s", var, letter)
    file_name <- sprintf("%s_%s", var, letter)
    # Collect individual data file from combined file
    df <- raw_data[[file_name]]
    if ("DMDBORN4" %in% colnames(df))
      df <- df %>% mutate(DMDBORN2 = DMDBORN4)
    # Define columns for the missing variables and fill them with NA
    var_absent <- setdiff(var_names, colnames(df))
    df[, var_absent] <- NA
    df <- df %>%
      select(dplyr::all_of(var_names)) %>%
      mutate(SEQN = as.character(SEQN))
    if (any(table(df$SEQN) != 1))
      stop("Multiple lines per respondent")
    # If it's DEMO, save to joined_data.
    # Otherwise, left-join with the previous/existing table
    if (var != names(keep_var)[1]) {
      df <- joined_data[[cycle_years[y]]] %>%
        merge(df, by = "SEQN", all.x = TRUE) %>%
        mutate_all(as.character)
    }
    joined_data[[cycle_years[y]]] <- df
  }
}

# Check that all columns line up with each other before rbind
cols_match <- sapply(joined_data, function(x) {
  vars <- unlist(keep_var)
  all(colnames(x) == c("SEQN", vars))
})

if (all(cols_match)) {
  all_data <- lapply(joined_data, function(dat) {
    # Process age, gender, and race variables for stratified reweighting
    dat <- dat %>%
      mutate(WTINT2YR = as.numeric(WTINT2YR)) %>% #interview weights
      mutate(WTMEC2YR = as.numeric(replace(WTMEC2YR, WTMEC2YR=="Not MEC Examined", 0))) %>% #exam weights
      mutate(RIDRETH=ifelse(is.na(RIDRETH3), RIDRETH1, RIDRETH3)) %>%
      mutate(IS_FEMALE = RIAGENDR == "Female") %>%
      mutate(IS_BLACK = RIDRETH1 == "Non-Hispanic Black") %>%
      mutate(RIDAGEYR = replace(RIDAGEYR, nchar(RIDAGEYR) > 2, "90") %>% as.integer) %>%
      mutate(cdc_age = cut(RIDAGEYR, right = FALSE,
                           breaks = cdc_breaks, labels = cdc_labels)) %>%
      mutate(PREGNANT = na_to_F(grepl("Yes", RIDEXPRG))) %>%
      mutate(EXCL_SMOKER = na_to_F(grepl("day", SMQ040, ignore.case = TRUE))) %>% #current smoker
      mutate(EXCL_SMOKED_100 = na_to_F(SMQ020=="Yes")) %>% # smoked at least 100 cigarettes ever
      mutate(EXCL_SMOKED_RECENT = na_to_F(SMQ690A=="Cigarettes") | # smoked X in the past 5 days
               na_to_F(SMQ690B=="Pipes") |
               na_to_F(SMQ690C=="Cigars"))
    
    dat <- dat %>%
      mutate(
        MISSING_AGE = RIDAGEYR < 3 | RIDAGEYR > 80,
        MISSING_DEMO = is.na(RIDRETH) | is.na(IS_FEMALE),
        MISSING_BMX = is.na(BMXHT) | is.na(BMXHT) | is.na(BMXBMI) | is.na(BMXWAIST),
        MISSING_SPIRO = !na_to_F(SPXNSTAT == "Complete Exam"),
        MISSING = MISSING_AGE | MISSING_DEMO | MISSING_BMX | MISSING_SPIRO,
        SPXNFEV1 = ifelse(MISSING, NA, SPXNFEV1),
        SPXNFVC = ifelse(MISSING, NA, SPXNFVC)
      ) %>%
      mutate(
        MCQ_ASTHMA = na_to_F(MCQ035 == "Yes"), #current asthma
        MCQ_BRONCHITIS = na_to_F(MCQ170K == "Yes"), #current chronic bronchitis
        MCQ_EMPHYSEMA = na_to_F(MCQ180G == "Yes"), #ever emphysema
        MCQ_LUNG_CANCER = !is.na(MCQ240N), #ever lung cancer
        # MCQ_HEART_FAILURE = na_to_F(MCQ160B == "Yes"), #current HF
        SYMP_COUGH = na_to_F(RDQ031=="Yes"), #coughing most days 3 mo.
        SYMP_WHEEZING = na_to_F(RDQ070=="Yes"), #wheezing/whistling in chest past year
        SYMP_PHLEGM = na_to_F(RDQ050=="Yes"), #phlegm most days 3 mo.; not found to affect PFTs
        #NIGHT_COUGH = na_to_F(RDQ140=="Yes"), #long-term dry cough at night; not found to affect PFTs
        EXCL_SMOKING = EXCL_SMOKER | EXCL_SMOKED_100 | EXCL_SMOKED_RECENT,
        EXCL_MEDICAL = MCQ_ASTHMA | MCQ_BRONCHITIS | MCQ_EMPHYSEMA |MCQ_LUNG_CANCER, #| MCQ_HEART_FAILURE,
        EXCL_SYMPTOMS = SYMP_WHEEZING | SYMP_PHLEGM | SYMP_COUGH,
        EXCL_AGE = !na_to_F(RIDAGEYR >= 6 & RIDAGEYR <= 79),
        EXCLUDED = EXCL_AGE | EXCL_SYMPTOMS | EXCL_MEDICAL | EXCL_SMOKING
      )
    
    dat <- suppressWarnings(dat %>%
                              mutate(DMD_Health_Ins = na_to_F(HIQ011=="Yes")) %>%
                              mutate(DMD_Priv_Ins = na_to_F(HIQ031A=="Covered by private insurance")) %>%
                              mutate(Medicare = na_to_F(HIQ031B=="Covered by Medicare")) %>%
                              mutate(Medicaid = na_to_F(HIQ031D=="Covered by Medicaid")) %>%
                              mutate(No_Ins_Yr = na_to_F(HIQ210=="Yes")) %>%
                              mutate(DMD_Home_Smoke = na_to_F(SMD410=="Yes")) %>%
                              mutate(Home_Smokers = as.integer(na_to_F(
                                SMD415 %>% replace(SMD415=="3 or more", 3), alt=0))) %>%
                              mutate(Home_Cigs = as.integer(na_to_F(
                                SMD430 %>% replace(SMD430=="40 or more", 40), alt=0))) %>%
                              mutate(DMD_FmHx_Asthma = na_to_F(MCQ300B=="Yes")) %>%
                              # mutate(DMD_Hospitalized = na_to_F(HUQ071=="Yes")) %>%
                              mutate(DMD_Hosp_Num = na_to_F(
                                HUD080 %>%
                                  replace(HUD080=="6 times or more", 6) %>%
                                  replace(HUD080=="Don't know", NA),
                                alt=0) %>% as.integer) %>%
                              mutate(Spanish_Lang = na_to_F(SIALANG=="Spanish")) %>%
                              mutate(Interp_Lang = na_to_F(SIAINTRP=="Yes")) %>%
                              mutate(DMD_English_Lang = na_to_F(SIALANG=="English") & !Interp_Lang) %>%
                              mutate(DMD_Veteran = na_to_F(DMQMILIZ == "Yes" | DMQMILIT == "Yes")) %>%
                              mutate(DMD_Military_HC = na_to_F(HIQ031F == "Covered by military health care")) %>%
                              mutate(DMD_Born_US = na_to_F(grepl("Born in 50 US", DMDBORN2))) %>%
                              mutate(Born_Spanish = na_to_F(DMDBORN2 %in% c("Born in Mexico","Born in Other Spanish Speaking Country"))) %>%
                              mutate(Born_Other = na_to_F(DMDBORN2 == "Born in Other Non-Spanish Speaking Country")) %>%
                              mutate(DMD_Citizen = na_to_F(grepl("naturalization", DMDCITZN))) %>%
                              mutate(DMD_IncPovRatio = as.numeric(INDFMPIR)) %>% #na_to_F(as.numeric(INDFMPIR)) %>%
                              mutate(DMD_Married = na_to_F(DMDMARTL=="Married")) %>%
                              # mutate(DMD_Start_College = na_to_F(grepl("college", DMDEDUC2, ignore.case = TRUE))) %>%
                              mutate(DMD_Finish_HS = na_to_F(grepl("college", DMDEDUC2, ignore.case = TRUE)) |
                                       na_to_F(grepl("GED", DMDEDUC2))) %>%
                              # mutate(DMD_Start_HS = na_to_F(!grepl("Less than 9th grade", DMDEDUC2))) %>%
                              mutate(DMD_Partnered = na_to_F(DMDMARTL %in% c("Married", "Living with partner"))) %>%
                              mutate(DMD_Years_in_US = DMDYRSUS %>%
                                       replace(is.na(DMDYRSUS), RIDAGEYR) %>%
                                       replace(DMDYRSUS %in% c("Don't know","Don't Know","Refused"), 20) %>%
                                       replace(grepl("Less than 1 y", DMDYRSUS), 0.5) %>%
                                       replace(grepl("less than 5 y", DMDYRSUS), 3) %>%
                                       replace(grepl("less than 10 y", DMDYRSUS), 7.5) %>%
                                       replace(grepl("less than 15 y", DMDYRSUS), 12.5) %>%
                                       replace(grepl("less than 20 y", DMDYRSUS), 17.5) %>%
                                       replace(grepl("less than 30 y", DMDYRSUS), 25) %>%
                                       replace(grepl("less than 40 y", DMDYRSUS), 35) %>%
                                       replace(grepl("less than 50 y", DMDYRSUS), 45) %>%
                                       replace(grepl("50 years or more", DMDYRSUS), 60) %>%
                                       as.numeric()
                              ) %>%
                              mutate(DMD_HH_Size = DMDHHSIZ %>%
                                       replace(DMDHHSIZ=="7 or more people in the Household", 7) %>%
                                       as.numeric) %>%
                              mutate(V_Work_Days = na_to_F(as.numeric(PAQ610), alt=0)) %>%
                              mutate(V_Work_Min = na_to_F(as.numeric(PAD615), alt=0)) %>%
                              mutate(V_Work_Min_Week = V_Work_Days * V_Work_Min) %>%
                              mutate(M_Work_Days = na_to_F(as.numeric(PAQ625), alt=0)) %>%
                              mutate(M_Work_Min = na_to_F(as.numeric(PAD630), alt=0)) %>%
                              mutate(M_Work_Min_Week = M_Work_Days * M_Work_Min) %>%
                              mutate(WB_Days = na_to_F(as.numeric(PAQ640), alt=0)) %>%
                              mutate(WB_Min = na_to_F(as.numeric(PAD645), alt=0)) %>%
                              mutate(WB_Min_Week = WB_Days * WB_Min) %>%
                              mutate(V_Rec_Days = na_to_F(as.numeric(PAQ655), alt=0)) %>%
                              mutate(V_Rec_Min = na_to_F(as.numeric(PAD660), alt=0)) %>%
                              mutate(V_Rec_Min_Week = V_Rec_Days * V_Rec_Min) %>%
                              mutate(M_Rec_Days = na_to_F(as.numeric(PAQ670), alt=0)) %>%
                              mutate(M_Rec_Min = na_to_F(as.numeric(PAD675), alt=0)) %>%
                              mutate(M_Rec_Min_Week = M_Rec_Days * M_Rec_Min) %>%
                              mutate(PAQ_MV_Min_Week = V_Work_Min_Week + M_Work_Min_Week +
                                       WB_Min_Week + V_Rec_Min_Week + M_Rec_Min_Week) %>%
                              mutate(PAQ_Sed_Min = na_to_F(as.numeric(PAD680), alt=0)) %>%
                              #mutate(OCQ_SHS_Smell = na_to_F(OCQ290Q > 0)) %>%
                              mutate(OCQ_Work_Smoke = na_to_F(OCQ275 == "Yes")) %>%
                              #mutate(OCQ_Mineral_Dust = na_to_F(OCQ510 == "Yes")) %>%
                              #mutate(OCQ_Organic_Dust = na_to_F(OCQ530 == "Yes")) %>%
                              #mutate(OCQ_Exhaust_Fumes = na_to_F(OCQ550 == "Yes")) %>%
                              #mutate(OCQ_Other_Fumes = na_to_F(OCQ570 == "Yes")) %>%
                              mutate(OCQ_Mineral_Dust_YR = na_to_F(OCQ520, alt=0) %>% as.integer) %>%
                              mutate(OCQ_Organic_Dust_YR = na_to_F(OCQ540, alt=0) %>% as.integer) %>%
                              mutate(OCQ_Exhaust_Fumes_YR = na_to_F(OCQ560, alt=0) %>% as.integer) %>%
                              mutate(OCQ_Other_Fumes_YR = na_to_F(OCQ580, alt=0) %>% as.integer) %>%
                              mutate(OCQ_Work_Air = (OCQ_Mineral_Dust_YR>0) | (OCQ_Organic_Dust_YR>0) | (OCQ_Exhaust_Fumes_YR>0) | (OCQ_Other_Fumes_YR>0)
                              ))
      #)
      #     dat <- dat %>% mutate(
      #       OCQ_EXPOSURE = OCQ_Mineral_Dust | OCQ_Organic_Dust |
      #         OCQ_Exhaust_Fumes | OCQ_Other_Fumes
      #     ) %>% mutate(
      #       OCQ_EXPOSURE_YR = pmax(OCQ_Mineral_Dust_YR, OCQ_Organic_Dust_YR,
      #                              OCQ_Exhaust_Fumes_YR, OCQ_Other_Fumes_YR)
      #     )
    
    weights <- dat %>%
      group_by(MISSING, RIDRETH1, IS_FEMALE, cdc_age) %>%
      dplyr::summarize(W = sum(WTMEC2YR), .groups = "drop") %>% as.data.frame()
    scaling_table <- weights[!weights$MISSING, -1]
    total_weights <- weights[!weights$MISSING, "W"] + weights[weights$MISSING, "W"]
    scaling_table$W <- total_weights/weights[!weights$MISSING, "W"]
    
    dat$WTMEC2YR_original <- dat$WTMEC2YR
    for (i in 1:nrow(scaling_table)) {
      row <- scaling_table[i,]
      to_change <- dat$RIDRETH1 == unlist(row["RIDRETH1"]) &
        dat$IS_FEMALE == unlist(row["IS_FEMALE"]) &
        dat$cdc_age == unlist(row["cdc_age"])
      dat[to_change, "WTMEC2YR"] <- dat[to_change, "WTMEC2YR"] * as.numeric(row["W"])
    }
    dat$WTMEC2YR[dat$MISSING] <- 0
    dat$WTMEC6YR <- dat$WTMEC2YR / (length(cycle_years_input)/2)
    dat$WTMEC6YR_original <- dat$WTMEC2YR_original / (length(cycle_years_input)/2)
    return(dat)
    
  }) %>% bind_rows(.id = "CYCLE_YEAR")
} else {
  stop("Column names do not match")
}

# Convert numeric cols to numeric type
convert_numeric <- function(dat) {
  for (col in num_cols) {
    if (col %in% colnames(dat))
      dat[,col] <- as.numeric(dat[,col])
  }
  return(dat)
}

nhanes4 <- convert_numeric(all_data) %>%
  mutate(
    WEIGHTS=WTMEC6YR_original/mean(WTMEC6YR_original),
    # Set to NA if MISSING or doesn't meet ATS standards
    FEV1_NQC = na_to_F(grepl("Exceeds ATS", SPXNQFV1)) | na_to_F(grepl("Meets ATS", SPXNQFV1)),
    FVC_NQC = na_to_F(grepl("Exceeds ATS", SPXNQFVC)) | na_to_F(grepl("Meets ATS", SPXNQFVC)),
    FEV1_BQC = na_to_F(grepl("Exceeds ATS", SPXBQFV1)) | na_to_F(grepl("Meets ATS", SPXBQFV1)),
    FVC_BQC = na_to_F(grepl("Exceeds ATS", SPXBQFVC)) | na_to_F(grepl("Meets ATS", SPXBQFVC)),
    SPXNFEV1 = ifelse(MISSING_SPIRO | SPXNFEV1>10000 | !FEV1_NQC, NA, SPXNFEV1),
    SPXNFVC = ifelse(MISSING_SPIRO | SPXNFVC>10000 | !FVC_NQC, NA, SPXNFVC),
    SPXBFEV1 = ifelse(MISSING_SPIRO | SPXBFEV1>10000 | !FEV1_BQC, NA, SPXBFEV1),
    SPXBFVC = ifelse(MISSING_SPIRO | SPXBFVC>10000 | !FVC_BQC, NA, SPXBFVC),
    # Flag as missing spiro if missing complete exam, FEV1, or FVC (both nl and bronchodilated)
    MISSING_SPIRO = MISSING_SPIRO |
      (is.na(SPXNFEV1) & is.na(SPXBFEV1)) |
      (is.na(SPXNFVC) & is.na(SPXBFVC)),
    MISSING = MISSING | MISSING_SPIRO,
    SPXBFEV1 = SPXBFEV1/1000,
    SPXBFVC = SPXBFVC/1000,
    SPXBFEV1FVC = SPXBFEV1/SPXBFVC,
    SPXNFEV1 = SPXNFEV1/1000,
    SPXNFVC = SPXNFVC/1000,
    SPXNFEV1FVC = SPXNFEV1/SPXNFVC,
    BRONCH_TEST = !is.na(SPXBFEV1),
    BRONCH_RESP = na_to_F((SPXBFEV1 > 1.12*SPXNFEV1 & SPXBFEV1 > SPXNFEV1 + 0.2) |
                            (SPXBFVC > 1.12*SPXNFVC & SPXBFVC > SPXNFVC + 0.2)),
    FEV1 = pmax(SPXNFEV1, SPXBFEV1, na.rm=TRUE),
    FVC = pmax(SPXNFVC, SPXBFVC, na.rm=TRUE),
    FEV1FVC = pmax(SPXNFEV1FVC, SPXBFEV1FVC, na.rm=TRUE),
    EXCL_MEDICAL = EXCL_MEDICAL | na_to_F(FEV1FVC < 0.7),
    SEQN = sprintf("B-%s", SEQN)
  ) %>% 
  mutate(
    RIDRETH = RIDRETH %>%
      replace(grepl("White", RIDRETH), "White") %>%
      replace(grepl("Black", RIDRETH), "Black") %>%
      replace(grepl("Chinese", RIDRETH), "Asian") %>%
      replace(grepl("Asian", RIDRETH), "Asian") %>%
      replace(grepl("Mexican", RIDRETH), "Hispanic") %>%
      replace(RIDRETH == "Other Hispanic", "Hispanic") %>%
      replace(RIDRETH == "Other", "Other") %>%
      replace(RIDRETH == "Other Race - Including Multi-Racial", "Other") %>%
      replace(RIDRETH == "Multiracial", "Other"),
    ETHNICITY = RIDRETH %>% 
      replace(RIDRETH=="Asian", "Other") %>% 
      replace(RIDRETH %in% c("Asian", "Other"), "Asian/Other")) %>% arrange(SEQN)

saveRDS(nhanes4, file = sprintf("processed/nhanes4_%s.rds", Sys.Date()))
setwd("/Users/aashnashah/Dropbox/Research/PFT-ML/")
write.table(nhanes4, file = sprintf("processed/nhanes4_%s.csv", Sys.Date()), quote = F, row.names = F, sep = "\t")

rm(raw_data, joined_data, all_data, df, keep_var)
