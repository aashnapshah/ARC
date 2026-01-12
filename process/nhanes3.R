# Pulmonary Function Reference Values Without Race Stratification in the United States

setwd("/Users/aashnashah/Dropbox/Research/ml_spirometry/")
source('helpers/params.R')
source('helpers/functions.R')

nhanes3_exam <- read.csv(file="data/nh3_adult_youth_exam/nhanes3_exam_processed.csv", header=TRUE) %>%
#  select(SEQN, WTPFHX6, MXPLANG, contains("DMAR"), contains("HS"), contains("BMP"),
#         contains("SDP"), contains("MYPB")) %>%
  mutate(
    MLANG = na_to_F(MXPLANG==1),
    BMXHT = ifelse(BMPHT==88888, NA, BMPHT),
    BMXWT = ifelse(BMPWT==888888, NA, BMPWT),
    BMXBMI = ifelse(BMPBMI==8888, NA, BMPBMI),
    BMXARMC = ifelse(BMPARMC==8888, NA, BMPARMC),
    BMXARML = ifelse(BMPARML==8888, NA, BMPARML),
    BMXLEG = ifelse(BMPLEG==8888, NA, BMPLEG),
    BMXWAIST = ifelse(BMPWAIST==88888, NA, BMPWAIST),
    BMXSITHT = ifelse(BMPSITHT==88888, NA, BMPSITHT),
    MXPLANG = na_to_F(MXPLANG==1)
    # PREGNANT = na_to_F(MAPF12R==1) | na_to_F(MYPC17==1) #children
  )

nhanes3_youth <- read.csv(file="data/nh3_adult_youth_exam/nhanes3_youth_processed.csv", header=TRUE) %>%
  select(SEQN, DMPPIR, HFA6XCR, HFA8R, HFB11, HFB13, HFA12, HYLANG, HFF1, HFF2R,
         contains("HYE"), contains("HYG")) %>%
  merge(nhanes3_exam, by="SEQN", all.x=TRUE, all.y=FALSE, sort = TRUE) %>%
  mutate(
    RIDRETH = factor(DMARETHN,
                     labels = c("Non-Hispanic White","Non-Hispanic Black","Mexican American","Other")) %>%
      as.character(),
    IS_BLACK = DMARETHN == 2,
    IS_FEMALE = factor(HSSEX, labels = c("FALSE","TRUE")) %>% as.logical(),
    RIDAGEYR = HSAGEIR * ifelse(na_to_F(HSAGEU==1), 1/12, 1),
    EXCL_SMOKED_100 = na_to_F(MYPB3==1),
    EXCL_SMOKED_RECENT = na_to_F(MYPB11==1) | na_to_F(MYPB27A==1) | na_to_F(MYPB27B==1),
    EXCL_SMOKER = na_to_F(MYPB5==1),
    # EXCL_EVER_SMOKER = na_to_F(MYPB1==1),
    EXCL_CIGARS_PIPES = FALSE,
    DMD_HH_Size = na_to_F(pmin(10, HSHSIZER), 0),
    DMD_Born_US = na_to_F(HFA6XCR==1),
    DMD_Finish_HS = na_to_F(HFA8R > 12),
    DMD_IncPovRatio = pmin(5, na_to_F(DMPPIR, alt=5)),
    DMD_Married = na_to_F(HFA12==1 | HFA12==2),
    DMD_Health_Ins = na_to_F(HFB13 == 1),
    DMD_Priv_Ins = na_to_F(HFB11==1),
    DMD_English_Lang = MLANG & na_to_F(HYLANG==1),
    DMD_Home_Smoke = na_to_F(HFF1==1),
    DMD_Veteran = FALSE,
    DMD_Military_HC = FALSE,
    MCQ_ASTHMA = na_to_F(HYE1G==1),
    MCQ_COPD = na_to_F(HYE1H==1), #MCQ_BRONCHITIS
    MCQ_HEART_FAILURE = FALSE,
    MCQ_LUNG_CANCER = FALSE,
    SYMP_WHEEZING = na_to_F(HYG8==1) | na_to_F(HYG12==1), #overall, with cough
    SYMP_COUGH = na_to_F(HYG2==1) | na_to_F(HYG6==1) | na_to_F(HYG7==555),
    SYMP_PHLEGM = na_to_F(HYG4==1)#,
    # SYMP_SOB = FALSE
  ) %>% select(
    SEQN, RIDRETH, IS_FEMALE, RIDAGEYR, WTPFHX6, starts_with("BMX"), starts_with("SDP"),
    starts_with("DMD"), starts_with("EXCL"), starts_with("MCQ"), starts_with("SYMP"),
  )

nhanes3_adult <- read.csv(file="data/nh3_adult_youth_exam/nhanes3_adult_processed.csv", header=TRUE) %>%
  select(SEQN, DMPPIR, HFA6XCR, HFA8R, HFA12, HFB11, HFB13,
         HALANG, HFA13, HFB9, HFF1, HFF2R,
         contains("HAR"), contains("HAC"), contains("HAL")) %>%
  merge(nhanes3_exam, by="SEQN", all.x=TRUE, all.y=FALSE, sort = TRUE) %>%
  mutate(
    RIDRETH = factor(DMARETHN,
                     labels = c("Non-Hispanic White","Non-Hispanic Black","Mexican American","Other")) %>%
      as.character(),
    IS_BLACK = DMARETHN == 2,
    IS_FEMALE = factor(HSSEX, labels = c("FALSE","TRUE")) %>% as.logical(),
    RIDAGEYR = HSAGEIR * ifelse(na_to_F(HSAGEU==1), 1/12, 1),
    EXCL_SMOKED_100 = na_to_F(HAR1==1) | na_to_F(MYPB3==1),
    EXCL_SMOKED_RECENT = na_to_F(MYPB11==1) | na_to_F(MYPB27A==1) | na_to_F(MYPB27B==1),
    EXCL_SMOKER = na_to_F(HAR3==1) | na_to_F(MYPB5==1),
    # EXCL_EVER_SMOKER = na_to_F(HAR1==1) | na_to_F(MYPB1==1),
    EXCL_CIGARS_PIPES = na_to_F(HAR24==1) | na_to_F(HAR27==1),
    DMD_HH_Size = na_to_F(pmin(10, HSHSIZER), 0),
    DMD_Born_US = na_to_F(HFA6XCR==1) | na_to_F(HFA6XCR==2),
    DMD_Finish_HS = na_to_F(HFA8R > 12),
    DMD_IncPovRatio = pmin(5, na_to_F(DMPPIR, alt=5)),
    DMD_Married = na_to_F(HFA12==1 | HFA12==2),
    DMD_Health_Ins = na_to_F(HFB13==1),
    DMD_Priv_Ins = na_to_F(HFB11==1),
    DMD_English_Lang = MLANG & na_to_F(HALANG==1),
    DMD_Home_Smoke = na_to_F(HFF1==1),
    DMD_Veteran = na_to_F(HFA13==1),
    DMD_Military_HC = na_to_F(HFB9==1),
    MCQ_ASTHMA = na_to_F(HAC1E==1),
    # MCQ_BRONCHITIS = na_to_F(HAC1F==1),
    # MCQ_EMPHYSEMA = na_to_F(HAC1G==1),
    MCQ_COPD = na_to_F(HAC1F==1) | na_to_F(HAC1G==1),
    MCQ_HEART_FAILURE = na_to_F(HAC1C==1),
    MCQ_LUNG_CANCER = na_to_F(HAC3OS == 15),
    SYMP_WHEEZING = na_to_F(HAL6==1 | HAL10==1), #overall, with cough
    SYMP_COUGH = na_to_F(HAL1==1),
    SYMP_PHLEGM = na_to_F(HAL3==1)#,
    # SYMP_SOB = na_to_F(HAL5==1)
  ) %>% select(
    SEQN, RIDRETH, IS_FEMALE, RIDAGEYR, WTPFHX6, starts_with("BMX"), starts_with("SDP"),
    starts_with("DMD"), starts_with("EXCL"), starts_with("MCQ"), starts_with("SYMP"),
  )

#https://wwwn.cdc.gov/nchs/data/nhanes3/9a/Readme.txt
nh3spiro_raw <- read.csv("data/nh3_spirometry/nh3spiro_processed.csv", header=TRUE)

missing_TQF <- between(nh3spiro_raw$SPPTQF, 10, 99) | nh3spiro_raw$SPPTQF > 109
#large extrapolated volume, cough, late PEF, or no plateau/<6s
missing_CQF <- sapply(nh3spiro_raw$SPPCQF, function(cqf) {
  any(intToBits(cqf)[1:4] > 0)
})

nh3spiro_raw$Rep_FEV1 <- sapply(nh3spiro_raw$SPPCQF, function(cqf) {
  intToBits(cqf)[6] == 0
})
nh3spiro_raw$Rep_FVC <- sapply(nh3spiro_raw$SPPCQF, function(cqf) {
  intToBits(cqf)[7] == 0
})

nh3spiro <- suppressWarnings(
  nh3spiro_raw %>%
    filter(!missing_CQF) %>%
    group_by(SEQN) %>%
    dplyr::summarize(
      FEV1=ifelse(length(SPPFEV1[Rep_FEV1]) >= 3, max(SPPFEV1[Rep_FEV1]), NA),
      FVC=ifelse(length(SPPFVC[Rep_FVC] >= 3), max(SPPFVC[Rep_FVC]), NA),
      .groups="drop") %>%
    filter(is.finite(FEV1), is.finite(FVC),
           FEV1<10000, FVC<10000)
)

nhanes3 <- rbind(nhanes3_adult, nhanes3_youth) %>%
  arrange(SEQN) %>%
  merge(nh3spiro, by="SEQN", all.x=TRUE, all.y=FALSE, sort = TRUE) %>%
  mutate(
    SEQN = sprintf("A-%s", SEQN),
    MISSING_AGE = RIDAGEYR < 3 | RIDAGEYR >= 90,
    MISSING_DEMO = is.na(RIDRETH) | is.na(IS_FEMALE),
    MISSING_BMX = is.na(BMXHT) | is.na(BMXHT) | is.na(BMXBMI) | is.na(BMXWAIST),
    MISSING_SPIRO = is.na(FEV1) | is.na(FVC),
    MISSING = MISSING_AGE | MISSING_DEMO | MISSING_BMX | MISSING_SPIRO,
    FEV1 = FEV1/1000,
    FVC = FVC/1000,
    FEV1FVC = pmin(FEV1/FVC, 1),
    CYCLE_YEAR = ifelse(SDPPHASE==1, "1988-1991", "1991-1994"),
    SDMVPSU = SDPPSU6,
    SDMVSTRA = SDPSTRA6,
    # WEIGHTS = WTPFHX6/mean(WTPFHX6, na.rm=TRUE),
    WTMEC6YR = WTPFHX6,
    SDMVPSU = SDPPSU6,
    SDMVSTRA = SDPSTRA6,
    EXCL_SMOKING = EXCL_SMOKED_100 | EXCL_SMOKED_RECENT |
      EXCL_SMOKER | EXCL_CIGARS_PIPES, # | EXCL_EVER_SMOKER,
    EXCL_MEDICAL = MCQ_ASTHMA | MCQ_COPD | MCQ_LUNG_CANCER | FEV1FVC < 0.7, # | MCQ_HEART_FAILURE,
    EXCL_SYMPTOMS = SYMP_WHEEZING | SYMP_COUGH | SYMP_PHLEGM, #| SYMP_SOB,
    EXCL_AGE = !na_to_F(RIDAGEYR >= 6 & RIDAGEYR <= 79),
    EXCLUDED = EXCL_AGE | EXCL_SYMPTOMS | EXCL_MEDICAL | EXCL_SMOKING
  ) %>% select(
    SEQN, CYCLE_YEAR, RIDRETH, IS_FEMALE, RIDAGEYR,
    WTMEC6YR, starts_with("SDMV"), FEV1, FVC, FEV1FVC,
    starts_with("MISSING"), starts_with("EXCL"),
    starts_with("MCQ"), starts_with("SYMP"),
    starts_with("BMX"), starts_with("DMD")
  ) %>% filter(!is.na(WTMEC6YR))

write.table(nhanes3, file = sprintf("processed/nhanes3_%s.csv", Sys.Date()), quote = F, row.names = F, sep = ",")

rm(nhanes3_adult, nhanes3_exam, nhanes3_youth, nh3spiro_raw)

