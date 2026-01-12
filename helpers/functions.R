# Multiethnic Spirometric Reference Values Without Race Stratification in the United States

select <- dplyr::select
summarize <- dplyr::summarize

na_cols <- function(df) {
  apply(df, 2, function(col) {
    mean(is.na(col))
  }) %>% sort %>% return()
}

save_table <- function(table, title) {
  write.table(table, sep = "\t", row.names = FALSE, quote = FALSE,
              sprintf("tables/raw/%s.tsv", title))
}

save_latex <- function(table, title) {
  save_kable(table, sprintf("tables/latex/%s.pdf", title))
}

# Convert missing to "No" or FALSE or anything else
na_to_F <- function(vec, alt = FALSE) {
  return(replace(vec, is.na(vec), alt))
}
# Rounded, weighted versions of aggregation functions
rw_quantile <- function(x, w, q, na.rm = TRUE) { return(wtd.quantile(as.numeric(x), q=q, na.rm = na.rm, weight = w)) }
rw_IQR <- function(x, w, digits=5) { round(rw_quantile(x, w, q=c(0.25, 0.75)) %>% diff, digits) }
rw_median <- function(x, w, digits=5) { round(rw_quantile(x, w, q=c(0.5)), digits) }
rw_percent <- function(x, w, digits=5) { round(100*weighted.mean(x, w, na.rm=TRUE), digits)}
rw_percentyes <- function(x, w, digits=5) { round(100*weighted.mean(x=="Yes", w, na.rm=TRUE), digits)}
r_sum <- function(x) { x %>% sum() %>% round() %>% return() }

# Succint footnote markers for latex tables
format_1 <- function(x) {format(x, digits=1, nsmall=1)}
comma_to_num <- function(x) {
  return(as.numeric(trimws(gsub(",","",x))))
}

get_funnel <- function(df) {
  c(nrow(df),
    sum(!df$MISSING_SPIRO),
    sum(!df$MISSING_SPIRO & !df$MISSING_DEMO & !df$MISSING_BMX),
    sum(!df$MISSING_SPIRO & !df$MISSING_DEMO & !df$MISSING_BMX &
          !df$EXCL_AGE),
    sum(!df$MISSING_SPIRO & !df$MISSING_DEMO & !df$MISSING_BMX &
          !df$EXCL_AGE & !df$EXCL_SMOKING),
    sum(!df$MISSING_SPIRO & !df$MISSING_DEMO & !df$MISSING_BMX &
          !df$EXCL_AGE & !df$EXCL_SMOKING & !df$EXCL_MEDICAL),
    sum(!df$MISSING_SPIRO & !df$MISSING_DEMO & !df$MISSING_BMX &
          !df$EXCL_AGE & !df$EXCL_SMOKING & !df$EXCL_MEDICAL & !df$EXCL_SYMPTOMS)
  )
}

put_num <- function(x, sigs=c(3,3,3)) {
  out <- quantile(x, probs = c(0.5, 0.25, 0.75), na.rm=TRUE)
  sprintf("%s (%s–%s)",
          out[1] %>% signif(sigs[1]),
          out[2] %>% signif(sigs[2]),
          out[3] %>% signif(sigs[3])
  ) %>% return()
}

put_perc <- function(x, sigs=3) {
  sprintf("%s (%s)", sum(x), 100 * signif(mean(x),sigs))
}

my_lm <- function(PFT, f, weights = NULL) {
  set.seed(100)
  use_formula <- as.formula(sprintf("%s ~ %s", PFT, f))
  lm(formula = use_formula, data=healthy[train,], weights = weights) %>% return()
}

my_gamlss <- function(PFT, f, weights = NULL) {
  set.seed(100)
  use_formula <- as.formula(sprintf("%s ~ %s", PFT, f))
  mod_all <- gamlss(
    formula = use_formula,
    sigma.formula = ~1+scs(RIDAGEYR, 2),
    nu.formula=~1,
    data=na.omit(healthy[train,]),
    weights = weights,
    family=BCCGo,
    method = RS()
  )
  return(mod_all)
  # form_text <- mod_all$mu.formula %>% as.character() %>% tail(1)
  # form_text <- gsub(", ", ",", form_text)
  # lower <- paste0("~", gsub(" \\+ [a-zA-Z0-9_\\(\\),]+:+[a-zA-Z0-9_\\(\\),]+","", form_text))
  # mod_final <- suppressWarnings(stepGAIC(
  #   object=mod_all,
  #   scope=list(
  #     lower=as.formula(lower),
  #     upper=mod_all$mu.formula
  #   ),
  #   k=log(mod_all$N),
  #   parameter="mu",
  #   direction="backward"
  # ))
  # return(mod_final)
}

my_rf <- function(PFT, f) {
  set.seed(100)
  use_formula <- as.formula(sprintf("%s ~ %s", PFT, f))
  randomForest(formula=use_formula, data=healthy[train,], ntree=1000, importance=FALSE)
}

wtd.zscore <- function(x1, y1, w, subset=TRUE) {
  numerator <- (y1[subset]/x1[subset,"mu"])^(x1[subset,"nu"])-1
  denominator <- x1[subset,"nu"] * x1[subset,"sigma"]
  return(weighted.mean(numerator/denominator, w[subset]))
}

wtd.rmse <- function(x1, y1, w, subset=TRUE) {
  return(sqrt(weighted.mean((y1[subset]-x1[subset])^2, w[subset], na.rm = TRUE)))
}

wtd.rmse_diff <- function(x1, y1, w, subset=TRUE) {
  subset_w <- healthy$RIDRETH == "White" & subset
  white <- sqrt(weighted.mean((y1[subset_w]-x1[subset_w])^2, w[subset_w], na.rm = TRUE))
  subset_b <- healthy$RIDRETH == "Black" & subset
  black <- sqrt(weighted.mean((y1[subset_b]-x1[subset_b])^2, w[subset_b], na.rm = TRUE))
  return(black-white)
}

wtd.percerror <- function(x1, y1, w, subset=TRUE) {
  return(weighted.mean(abs(y1[subset]-x1[subset])/y1[subset], w[subset], na.rm = TRUE))
}

wtd.percerror_diff <- function(x1, y1, w, subset=TRUE) {
  subset_w <- healthy$RIDRETH == "White" & subset
  white <- weighted.mean(abs(y1[subset_w]-x1[subset_w])/y1[subset_w], w[subset_w], na.rm = TRUE)
  subset_b <- healthy$RIDRETH == "Black" & subset
  black <- weighted.mean(abs(y1[subset_b]-x1[subset_b])/y1[subset_b], w[subset_b], na.rm = TRUE)
  return(black-white)
}

wtd.repeatability <- function(x1, y1, w, subset=TRUE) {
  return(weighted.mean(abs(y1[subset]-x1[subset]) > 0.150, w[subset], na.rm = TRUE))
}

wtd.repeatability_diff <- function(x1, y1, w, subset=TRUE) {
  subset_w <- (healthy$RIDRETH == "White") & subset
  white <- weighted.mean(abs(y1[subset_w]-x1[subset_w]) > 0.150, w[subset_w], na.rm = TRUE)
  subset_b <- (healthy$RIDRETH == "Black") & subset
  black <- weighted.mean(abs(y1[subset_b]-x1[subset_b]) > 0.150, w[subset_b], na.rm = TRUE)
  return(black-white)
}

wtd.p20 <- function(x1, y1, w, subset=TRUE) {
  return(weighted.mean(between(x1[subset]/y1[subset], 0.8, 1.2), w[subset], na.rm = TRUE))
}

wtd.p20_diff <- function(x1, y1, w, subset=TRUE) {
  subset_w <- (healthy$RIDRETH == "White") & subset
  white <- weighted.mean(between(x1[subset_w]/y1[subset_w], 0.8, 1.2), w[subset_w], na.rm = TRUE)
  subset_b <- (healthy$RIDRETH == "Black") & subset
  black <- weighted.mean(between(x1[subset_b]/y1[subset_b], 0.8, 1.2), w[subset_b], na.rm = TRUE)
  return(black-white)
}


get_gamlss_params <- function(mod, newdata) {
  sapply(c("mu","sigma","nu"), function(param) {
    predict(object=mod, newdata=newdata, parameter=param, type="response")
  }, USE.NAMES = TRUE) %>% return()
}

get_quantile <- function(params, q) {
  return(params[,"mu"]*exp(log(qnorm(q)*params[,"sigma"]*params[,"nu"]+1)/params[,"nu"]))
}

get_z_score <- function(params, values) {
  return(((values/params[,"mu"])^(params[,"nu"])-1)/(params[,"nu"] * params[,"sigma"]))
}

# Formats values for specified formulae into plottable data frames
repeatability <- function(x, y) {
  return(abs(x - y) < 0.200)
}

paren_format <- function(df, digits=2, sigs=3, sep="–") {
  # Generates formatted table: rows = race/ethnicity, 1 col = val (ci_l - ci_u)
  df %>%
    apply(1, function(row) {
      # sigs2 <- sigs
      # # if (max(row) < 10) sigs <- sigs2 <- 3
      # if (signif(round(row[2], sigs2-1), sigs2-1) == signif(round(row[3], sigs2-1), sigs2-1)) {
      #   sigs2 <- sigs + 1
      # }
      # if (max(row) > 100) {
      #   sigs <- ceiling(log10(max(row)))
      #   sigs2 <- sigs + 1
      # }
      sprintf("%s (%s%s%s)",
              signif(round(row[1], sigs), sigs),
              signif(round(row[2], sigs), sigs),
              # signif(round(row[2], sigs2-1), sigs2-1),
              sep,
              signif(round(row[3], sigs), sigs))
      # signif(round(row[3], sigs2-1), sigs2-1))
    }) %>%
    as.data.frame() %>%
    return()
}

get_waitlist <- function(FVC_perc, age, bmi, sick=FALSE) {
  betas <- 0.015097*age - 0.051781*bmi - 0.019675*FVC_perc +
    #pHTN, diabetes, fxn status (some, total), CMV, 6MWT, incr PCO2
    sick*(2.376700 + 0.158821 + 0.182250 + 0.115024 + 1.2138 + 0.330752 + 0.076370) +
    #PA systolic, PCO2
    ifelse(sick, 90, 20)*0.015889 + ifelse(sick, 80, 40)*0.005448
  sapply(betas, function(beta) {
    sum(waitlist_survival^exp(beta))
  })
}

get_post_tx <- function(FVC_perc, age, bmi, sick=FALSE) {
  alphas <- 0.003510*age - 0.002751*FVC_perc +
    #PCWP, CMV, pHTN
    sick*(0.033046 + 0.312846 + 0.623207) +
    #SCr, fxn status (total),
    ifelse(sick, 5, 1)*0.061986 + ifelse(sick, 0, -0.488525)
  sapply(alphas, function(alpha) {
    sum(post_tx_survival^exp(alpha))
  })
}

get_LAS <- function(FVC_perc, age, bmi, sick=FALSE) {
  100*(get_post_tx(FVC_perc, age, bmi, sick) - 2*get_waitlist(FVC_perc, age, bmi, sick) + 730)/1095
}

get_mean_tx <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean(f=sprintf("~%s(100*%s/%s_PRD_NH3, RIDAGEYR, BMXBMI)", fx, PFT, PFT),
             des, paren, end_op), #NHANES original
    get_mean(f=sprintf("~%s(100*%s/%s_PRD_NH3_bl, RIDAGEYR, BMXBMI)", fx, PFT, PFT),
             des, paren, end_op), #NHANES blended
    get_mean(f=sprintf("~%s(100*%s/%s_PRD_GLI, RIDAGEYR, BMXBMI)", fx, PFT, PFT),
             des, paren, end_op), #GLI original
    get_mean(f=sprintf("~%s(100*%s/%s_PRD_GLI_bl, RIDAGEYR, BMXBMI)", fx, PFT, PFT),
             des, paren, end_op), #GLI blended
    get_mean(f=sprintf("~%s(100*%s/%s_Age_Sex_HT_Race_GAMLSS, RIDAGEYR, BMXBMI)", fx, PFT, PFT),
             des, paren, end_op), #New race-based
    get_mean(f=sprintf("~%s(100*%s/%s_Age_Sex_HT_WT_BMI_GAMLSS, RIDAGEYR, BMXBMI)", fx, PFT, PFT),
             des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_mean_tx_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean_bs(f=sprintf(
      "~%s(100*des$variables$%s/des$variables$%s_PRD_NH3, RIDAGEYR, BMXBMI)",
      fx, PFT, PFT), des, paren, end_op), #NHANES original
    get_mean_bs(f=sprintf(
      "~%s(100*des$variables$%s/des$variables$%s_PRD_NH3_bl, RIDAGEYR, BMXBMI)",
      fx, PFT, PFT), des, paren, end_op), #NHANES blended
    get_mean_bs(f=sprintf(
      "~%s(100*des$variables$%s/des$variables$%s_PRD_GLI, RIDAGEYR, BMXBMI)",
      fx, PFT, PFT), des, paren, end_op), #GLI original
    get_mean_bs(f=sprintf(
      "~%s(100*des$variables$%s/des$variables$%s_PRD_GLI_bl, RIDAGEYR, BMXBMI)",
      fx, PFT, PFT), des, paren, end_op), #GLI blended
    get_mean_bs(f=sprintf(
      "~%s(100*des$variables$%s/des$variables$%s_Age_Sex_HT_Race_GAMLSS, RIDAGEYR, BMXBMI)",
      fx, PFT, PFT), des, paren, end_op), #New race-based
    get_mean_bs(f=sprintf(
      "~%s(100*des$variables$%s/des$variables$%s_Age_Sex_HT_WT_BMI_GAMLSS, RIDAGEYR, BMXBMI)",
      fx, PFT, PFT), des, paren, end_op),  #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_mean_waitlist <- function(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_waitlist") {
  get_mean_tx(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_waitlist")
}
get_mean_waitlist_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_waitlist") {
  get_mean_tx_bs(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_waitlist")
}

get_mean_post_tx <- function(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_post_tx") {
  get_mean_tx(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_post_tx")
}
get_mean_post_tx_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_post_tx") {
  get_mean_tx_bs(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_post_tx")
}

get_mean_LAS <- function(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_LAS") {
  get_mean_tx(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_LAS")
}
get_mean_LAS_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_LAS") {
  get_mean_tx_bs(des, PFT, rPFT, paren=NULL, end_op=NULL, fx="get_LAS")
}


get_mean_zs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean(f=sprintf("~I(%s_zs_NH3)", PFT),
             des, paren, end_op), #NHANES original
    get_mean(f=sprintf("~I(%s_zs_NH3_bl)", PFT),
             des, paren, end_op), #NHANES blended
    get_mean(f=sprintf("~I(%s_zs_GLI)", PFT),
             des, paren, end_op), #GLI original
    get_mean(f=sprintf("~I(%s_zs_GLI_bl)", PFT),
             des, paren, end_op), #GLI blended
    get_mean(f=sprintf("~I(zs_%s_Age_Sex_HT_Race_GAMLSS)", PFT),
             des, paren, end_op), #New race-based
    get_mean(f=sprintf("~I(zs_%s_Age_Sex_HT_WT_BMI_GAMLSS)", PFT),
             des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_mean_zs_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean_bs(f=sprintf("des$variables$%s_zs_NH3",
                          PFT), des, paren, end_op), #NHANES original
    get_mean_bs(f=sprintf("des$variables$%s_zs_NH3_bl",
                          PFT), des, paren, end_op), #NHANES original #NHANES blended
    get_mean_bs(f=sprintf("des$variables$%s_zs_GLI",
                          PFT), des, paren, end_op), #NHANES original #GLI original
    get_mean_bs(f=sprintf("des$variables$%s_zs_GLI_bl",
                          PFT), des, paren, end_op), #NHANES original #GLI blended
    get_mean_bs(f=sprintf("des$variables$zs_%s_Age_Sex_HT_Race_GAMLSS",
                          PFT), des, paren, end_op), #NHANES original des, paren, end_op), #New race-based
    get_mean_bs(f=sprintf("des$variables$zs_%s_Age_Sex_HT_WT_BMI_GAMLSS",
                          PFT), des, paren, end_op) #NHANES originaldes, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}


below_LLN <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop(f=sprintf("~I(%s < %s_LLN_NH3)", rPFT, PFT), des, paren), #NHANES original
    get_prop(f=sprintf("~I(%s < %s_LLN_NH3_bl)", rPFT, PFT), des, paren), #NHANES blended
    get_prop(f=sprintf("~I(%s < %s_LLN_GLI)", rPFT, PFT), des, paren), #GLI original
    get_prop(f=sprintf("~I(%s < %s_LLN_GLI_bl)", rPFT, PFT), des, paren), #GLI blended
    get_prop(f=sprintf("~I(%s < LLN_%s_Age_Sex_HT_Race_GAMLSS)", rPFT, PFT), des, paren), #New race-based
    get_prop(f=sprintf("~I(%s < LLN_%s_Age_Sex_HT_WT_BMI_GAMLSS)", rPFT, PFT), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}


below_LLN_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop_bs(f=sprintf("des$variables$%s < des$variables$%s_LLN_NH3", PFT, PFT), des, paren), #NHANES original
    get_prop_bs(f=sprintf("des$variables$%s < des$variables$%s_LLN_NH3_bl", PFT, PFT), des, paren), #NHANES blended
    get_prop_bs(f=sprintf("des$variables$%s < des$variables$%s_LLN_GLI", PFT, PFT), des, paren), #GLI original
    get_prop_bs(f=sprintf("des$variables$%s < des$variables$%s_LLN_GLI_bl", PFT, PFT), des, paren), #GLI blended
    get_prop_bs(f=sprintf("des$variables$%s < des$variables$LLN_%s_Age_Sex_HT_Race_GAMLSS", PFT, PFT), des, paren), #New race-based
    get_prop_bs(f=sprintf("des$variables$%s < des$variables$LLN_%s_Age_Sex_HT_WT_BMI_GAMLSS", PFT, PFT), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

within_200 <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop(f=sprintf("~repeatability(%s, %s_PRD_NH3)", rPFT, PFT), des, paren), #NHANES original
    get_prop(f=sprintf("~repeatability(%s, %s_PRD_NH3_bl)", rPFT, PFT), des, paren), #NHANES blended
    get_prop(f=sprintf("~repeatability(%s, %s_PRD_GLI)", rPFT, PFT), des, paren), #GLI original
    get_prop(f=sprintf("~repeatability(%s, %s_PRD_GLI_bl)", rPFT, PFT), des, paren), #GLI blended
    get_prop(f=sprintf("~repeatability(%s, Pred_%s_Age_Sex_HT_Race_GAMLSS)", rPFT, PFT), des, paren), #New race-based
    get_prop(f=sprintf("~repeatability(%s, Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS)", rPFT, PFT), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

within_200_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop_bs(f=sprintf("repeatability(des$variables$%s, des$variables$%s_PRD_NH3)", PFT, PFT), des, paren), #NHANES original
    get_prop_bs(f=sprintf("repeatability(des$variables$%s, des$variables$%s_PRD_NH3_bl)", PFT, PFT), des, paren), #NHANES blended
    get_prop_bs(f=sprintf("repeatability(des$variables$%s, des$variables$%s_PRD_GLI)", PFT, PFT), des, paren), #GLI original
    get_prop_bs(f=sprintf("repeatability(des$variables$%s, des$variables$%s_PRD_GLI_bl)", PFT, PFT), des, paren), #GLI blended
    get_prop_bs(f=sprintf("repeatability(des$variables$%s, des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS)", PFT, PFT), des, paren), #New race-based
    get_prop_bs(f=sprintf("repeatability(des$variables$%s, des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS)", PFT, PFT), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_perc_error <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop(f=sprintf("~I(abs(%s_PRD_NH3-%s)/%s)", rPFT, PFT, rPFT),
             des, paren, end_op), #NHANES original
    get_prop(f=sprintf("~I(abs(%s_PRD_NH3_bl - %s)/%s)", rPFT, PFT, rPFT),
             des, paren, end_op), #NHANES blended
    get_prop(f=sprintf("~I(abs(%s_PRD_GLI - %s)/%s)", rPFT, PFT, rPFT),
             des, paren, end_op), #GLI original
    get_prop(f=sprintf("~I(abs(%s_PRD_GLI_bl - %s)/%s)", rPFT, PFT, rPFT),
             des, paren, end_op), #GLI blended
    get_prop(f=sprintf("~I(abs(Pred_%s_Age_Sex_HT_Race_GAMLSS - %s)/%s)", rPFT, PFT, rPFT),
             des, paren, end_op), #New race-based
    get_prop(f=sprintf("~I(abs(Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS - %s)/%s)", rPFT, PFT, rPFT),
             des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_perc_error_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop_bs(f=sprintf("abs(des$variables$%s_PRD_NH3 - des$variables$%s)/des$variables$%s",
                          PFT, PFT, PFT), des, paren, end_op), #NHANES original
    get_prop_bs(f=sprintf("abs(des$variables$%s_PRD_NH3_bl - des$variables$%s)/des$variables$%s",
                          PFT, PFT, PFT), des, paren, end_op), #NHANES blended
    get_prop_bs(f=sprintf("abs(des$variables$%s_PRD_GLI - des$variables$%s)/des$variables$%s",
                          PFT, PFT, PFT), des, paren, end_op), #GLI original
    get_prop_bs(f=sprintf("abs(des$variables$%s_PRD_GLI_bl - des$variables$%s)/des$variables$%s",
                          PFT, PFT, PFT), des, paren, end_op), #GLI blended
    get_prop_bs(f=sprintf("abs(des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS - des$variables$%s)/des$variables$%s",
                          PFT, PFT, PFT), des, paren, end_op), #New race-based
    get_prop_bs(f=sprintf("abs(des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS - des$variables$%s)/des$variables$%s",
                          PFT, PFT, PFT), des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_MAE <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean(f=sprintf("~I(abs(%s - %s_PRD_NH3))", rPFT, PFT),
             des, paren, end_op), #NHANES original
    get_mean(f=sprintf("~I(abs(%s - %s_PRD_NH3_bl))", rPFT, PFT),
             des, paren, end_op), #NHANES blended
    get_mean(f=sprintf("~I(abs(%s - %s_PRD_GLI))", rPFT, PFT),
             des, paren, end_op), #GLI original
    get_mean(f=sprintf("~I(abs(%s - %s_PRD_GLI_bl))", rPFT, PFT),
             des, paren, end_op), #GLI blended
    get_mean(f=sprintf("~I(abs(%s - Pred_%s_Age_Sex_HT_Race_GAMLSS))", rPFT, PFT),
             des, paren, end_op), #New race-based
    get_mean(f=sprintf("~I(abs(%s - Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS))", rPFT, PFT),
             des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}


get_MAE_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean_bs(f=sprintf("abs(des$variables$%s - des$variables$%s_PRD_NH3)",
                          PFT, PFT), des, paren, end_op), #NHANES original
    get_mean_bs(f=sprintf("abs(des$variables$%s - des$variables$%s_PRD_NH3_bl)",
                          PFT, PFT), des, paren, end_op), #NHANES original #NHANES blended
    get_mean_bs(f=sprintf("abs(des$variables$%s - des$variables$%s_PRD_GLI)",
                          PFT, PFT), des, paren, end_op), #NHANES original #GLI original
    get_mean_bs(f=sprintf("abs(des$variables$%s - des$variables$%s_PRD_GLI_bl)",
                          PFT, PFT), des, paren, end_op), #NHANES original #GLI blended
    get_mean_bs(f=sprintf("abs(des$variables$%s - des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS)",
                          PFT, PFT), des, paren, end_op), #NHANES original des, paren, end_op), #New race-based
    get_mean_bs(f=sprintf("abs(des$variables$%s - des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS)",
                          PFT, PFT), des, paren, end_op) #NHANES originaldes, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}


get_MSE <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean(f=sprintf("~I((%s - %s_PRD_NH3)^2)", rPFT, PFT),
             des, paren, end_op), #NHANES original
    get_mean(f=sprintf("~I((%s - %s_PRD_NH3_bl)^2)", rPFT, PFT),
             des, paren, end_op), #NHANES blended
    get_mean(f=sprintf("~I((%s - %s_PRD_GLI)^2)", rPFT, PFT),
             des, paren, end_op), #GLI original
    get_mean(f=sprintf("~I((%s - %s_PRD_GLI_bl)^2)", rPFT, PFT),
             des, paren, end_op), #GLI blended
    get_mean(f=sprintf("~I((%s - Pred_%s_Age_Sex_HT_Race_GAMLSS)^2)", rPFT, PFT),
             des, paren, end_op), #New race-based
    get_mean(f=sprintf("~I((%s - Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS)^2)", rPFT, PFT),
             des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_MSE_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean_bs(f=sprintf("(des$variables$%s - des$variables$%s_PRD_NH3)^2", PFT, PFT),
                des, paren, end_op), #NHANES original
    get_mean_bs(f=sprintf("(des$variables$%s - des$variables$%s_PRD_NH3_bl)^2", PFT, PFT),
                des, paren, end_op), #NHANES blended
    get_mean_bs(f=sprintf("(des$variables$%s - des$variables$%s_PRD_GLI)^2", PFT, PFT),
                des, paren, end_op), #GLI original
    get_mean_bs(f=sprintf("(des$variables$%s - des$variables$%s_PRD_GLI_bl)^2", PFT, PFT),
                des, paren, end_op), #GLI blended
    get_mean_bs(f=sprintf("(des$variables$%s - des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS)^2", PFT, PFT),
                des, paren, end_op), #New race-based
    get_mean_bs(f=sprintf("(des$variables$%s - des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS)^2", PFT, PFT),
                des, paren, end_op) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_p20 <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  s <- 0.2
  results <- list(
    get_prop(f=sprintf("~between(%s/%s_PRD_NH3, 0, 1-%s)", PFT, rPFT, s, s), des, paren), #NHANES original
    get_prop(f=sprintf("~between(%s/%s_PRD_NH3_bl, 0, 1-%s)", PFT, rPFT, s, s), des, paren), #NHANES blended
    get_prop(f=sprintf("~between(%s/%s_PRD_GLI, 0, 1-%s)", PFT, rPFT, s, s), des, paren), #GLI original
    get_prop(f=sprintf("~between(%s/%s_PRD_GLI_bl, 0, 1-%s)", PFT, rPFT, s, s), des, paren), #GLI blended
    get_prop(f=sprintf("~between(%s/Pred_%s_Age_Sex_HT_Race_GAMLSS, 0, 1-%s)",
                       PFT, rPFT, s, s), des, paren), #New race-based
    get_prop(f=sprintf("~between(%s/Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS, 0, 1-%s)",
                       PFT, rPFT, s, s), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}


get_p20_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  s <- 0.2
  results <- list(
    get_prop_bs(f=sprintf("between(des$variables$%s/des$variables$%s_PRD_NH3, 0, 1-%s)", PFT, PFT, s, s), des, paren), #NHANES original
    get_prop_bs(f=sprintf("between(des$variables$%s/des$variables$%s_PRD_NH3_bl, 0, 1-%s)", PFT, PFT, s, s), des, paren), #NHANES blended
    get_prop_bs(f=sprintf("between(des$variables$%s/des$variables$%s_PRD_GLI, 0, 1-%s)", PFT, PFT, s, s), des, paren), #GLI original
    get_prop_bs(f=sprintf("between(des$variables$%s/des$variables$%s_PRD_GLI_bl, 0, 1-%s)", PFT, PFT, s, s), des, paren), #GLI blended
    get_prop_bs(f=sprintf("between(des$variables$%s/des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS, 0, 1-%s)",
                          PFT, PFT, s, s), des, paren), #New race-based
    get_prop_bs(f=sprintf("between(des$variables$%s/des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS, 0, 1-%s)",
                          PFT, PFT, s, s), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_perc_pred <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean(f=sprintf("~I(%s/%s_PRD_NH3)", rPFT, PFT), des, paren), #NHANES original
    get_mean(f=sprintf("~I(%s/%s_PRD_NH3_bl)", rPFT, PFT), des, paren), #NHANES blended
    get_mean(f=sprintf("~I(%s/%s_PRD_GLI)", rPFT, PFT), des, paren), #GLI original
    get_mean(f=sprintf("~I(%s/%s_PRD_GLI_bl)", rPFT, PFT), des, paren), #GLI blended
    get_mean(f=sprintf("~I(%s/Pred_%s_Age_Sex_HT_Race_GAMLSS)", rPFT, PFT), des, paren), #New race-based
    get_mean(f=sprintf("~I(%s/Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS)", rPFT, PFT), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_perc_pred_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean_bs(f=sprintf("des$variables$%s/des$variables$%s_PRD_NH3", PFT, PFT), des, paren), #NHANES original
    get_mean_bs(f=sprintf("des$variables$%s/des$variables$%s_PRD_NH3_bl", PFT, PFT), des, paren), #NHANES blended
    get_mean_bs(f=sprintf("des$variables$%s/des$variables$%s_PRD_GLI", PFT, PFT), des, paren), #GLI original
    get_mean_bs(f=sprintf("des$variables$%s/des$variables$%s_PRD_GLI_bl", PFT, PFT), des, paren), #GLI blended
    get_mean_bs(f=sprintf("des$variables$%s/des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS", PFT, PFT), des, paren), #New race-based
    get_mean_bs(f=sprintf("des$variables$%s/des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS", PFT, PFT), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_GOLD_1p <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  sep <- 0.8 # GOLD 1+
  rPFT <- "SPXBFEV1"
  long_fev1fvc <- "SPXBFEV1FVC"
  results <- list(
    get_prop(f=sprintf("~I(na_to_F(%s/%s_PRD_NH3 < %s & %s < 0.7))",
                       rPFT, PFT, sep, long_fev1fvc), des, paren), #NHANES original
    get_prop(f=sprintf("~I(na_to_F(%s/%s_PRD_NH3_bl < %s & %s < 0.7))",
                       rPFT, PFT, sep, long_fev1fvc), des, paren), #NHANES blended
    get_prop(f=sprintf("~I(na_to_F(%s/%s_PRD_GLI < %s & %s < 0.7))",
                       rPFT, PFT, sep, long_fev1fvc), des, paren), #GLI original
    get_prop(f=sprintf("~I(na_to_F(%s/%s_PRD_GLI_bl < %s & %s < 0.7))",
                       rPFT, PFT, sep, long_fev1fvc), des, paren), #GLI blended
    get_prop(f=sprintf("~I(na_to_F(%s/Pred_%s_Age_Sex_HT_Race_GAMLSS < %s & %s < 0.7))",
                       rPFT, PFT, sep, long_fev1fvc), des, paren), #New race-based
    get_prop(f=sprintf("~I(na_to_F(%s/Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS < %s & %s < 0.7))",
                       rPFT, PFT, sep, long_fev1fvc), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_GOLD_1p_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  sep <- 1 # GOLD 1+
  long_rPFT <- "des$variables$SPXBFEV1/1000"
  long_fev1fvc <- "des$variables$SPXBFEV1/des$variables$SPXBFVC"
  results <- list(
    get_prop_bs(f=sprintf("%s/des$variables$%s_PRD_NH3 < %s & %s < 0.7",
                          long_rPFT, PFT, sep, long_fev1fvc), des, paren), #NHANES original
    get_prop_bs(f=sprintf("%s/des$variables$%s_PRD_NH3_bl < %s & %s < 0.7",
                          long_rPFT, PFT, sep, long_fev1fvc), des, paren), #NHANES blended
    get_prop_bs(f=sprintf("%s/des$variables$%s_PRD_GLI < %s & %s < 0.7",
                          long_rPFT, PFT, sep, long_fev1fvc), des, paren), #GLI original
    get_prop_bs(f=sprintf("%s/des$variables$%s_PRD_GLI_bl < %s & %s < 0.7",
                          long_rPFT, PFT, sep, long_fev1fvc), des, paren), #GLI blended
    get_prop_bs(f=sprintf("%s/des$variables$Pred_%s_Age_Sex_HT_Race_GAMLSS < %s & %s < 0.7",
                          long_rPFT, PFT, sep, long_fev1fvc), des, paren), #New race-based
    get_prop_bs(f=sprintf("%s/des$variables$Pred_%s_Age_Sex_HT_WT_BMI_GAMLSS < %s & %s < 0.7",
                          long_rPFT, PFT, sep, long_fev1fvc), des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_VA_disability <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop(f="~I(VA_Disability_NH3 > 0)", des, paren), #NHANES original
    get_prop(f="~I(VA_Disability_NH3_bl > 0)", des, paren), #NHANES blended
    get_prop(f="~I(VA_Disability_GLI > 0)", des, paren), #GLI original
    get_prop(f="~I(VA_Disability_GLI_bl > 0)", des, paren), #GLI blended
    get_prop(f="~I(VA_Disability_GAMLSS > 0)", des, paren), #New race-based
    get_prop(f="~I(VA_Disability_GAMLSS_bl > 0)", des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_VA_disability_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_prop_bs(f="des$variables$VA_Disability_NH3 > 0", des, paren), #NHANES original
    get_prop_bs(f="des$variables$VA_Disability_NH3_bl > 0", des, paren), #NHANES blended
    get_prop_bs(f="des$variables$VA_Disability_GLI > 0", des, paren), #GLI original
    get_prop_bs(f="des$variables$VA_Disability_GLI_bl > 0", des, paren), #GLI blended
    get_prop_bs(f="des$variables$VA_Disability_GAMLSS > 0", des, paren), #New race-based
    get_prop_bs(f="des$variables$VA_Disability_GAMLSS_bl > 0", des, paren) #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_VA_payment <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean(f="~I(12*va_payment(VA_Disability_NH3, DMD_Married))", des, paren), #NHANES original
    get_mean(f="~I(12*va_payment(VA_Disability_NH3_bl, DMD_Married))", des, paren), #NHANES blended
    get_mean(f="~I(12*va_payment(VA_Disability_GLI, DMD_Married))", des, paren), #GLI original
    get_mean(f="~I(12*va_payment(VA_Disability_GLI_bl, DMD_Married))", des, paren), #GLI blended
    get_mean(f="~I(12*va_payment(VA_Disability_GAMLSS, DMD_Married))", des, paren), #New race-based
    get_mean(f="~I(12*va_payment(VA_Disability_GAMLSS_bl, DMD_Married))", des, paren)#, #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

get_VA_payment_bs <- function(des, PFT, rPFT, paren=NULL, end_op=NULL) {
  results <- list(
    get_mean_bs(f="va_payment(des$variables$VA_Disability_NH3, des$variables$DMD_Married)",
                des, paren), #NHANES original
    get_mean_bs(f="va_payment(des$variables$VA_Disability_NH3_bl, des$variables$DMD_Married)",
                des, paren), #NHANES blended
    get_mean_bs(f="va_payment(des$variables$VA_Disability_GLI, des$variables$DMD_Married)",
                des, paren), #GLI original
    get_mean_bs(f="va_payment(des$variables$VA_Disability_GLI_bl, des$variables$DMD_Married)",
                des, paren), #GLI blended
    get_mean_bs(f="va_payment(des$variables$VA_Disability_GAMLSS, des$variables$DMD_Married)",
                des, paren), #New race-based
    get_mean_bs(f="va_payment(des$variables$VA_Disability_GAMLSS_bl, des$variables$DMD_Married)",
                des, paren)#, #New race-free
  )
  if (!is.null(paren)) {
    return(data.frame(results))
  }
  return(bind_rows(results))
}

# Computes value and CIs for an arbitrary specified formula by race/ethnicity
get_mean_bs <- function(f, des, paren = FALSE, end_op = NULL, names = FALSE) {
  y <- eval(parse(text=f))
  get_bs <- function(y, q) {
    set.seed(100)
    sapply(1:1000, function(i) {
      mean(sample(y, size=length(y), replace=TRUE))
    }) %>% quantile(probs=q)
  }
  overall <- c(mean(y), get_bs(y, c(0.025, 0.975)))
  
  by_race_eth <- des$variables %>% data.frame(y=y) %>% group_by(RIDRETH) %>%
    dplyr::summarize(
      mean=mean(y),
      ci_l=get_bs(y, 0.025),
      ci_u=get_bs(y, 0.975),
      .groups="drop"
    ) %>% data.frame
  output <- rbind(
    c(overall),
    by_race_eth %>% select(-RIDRETH)
  )
  rownames(output) <- c("Overall", by_race_eth$RIDRETH)
  if (!is.null(end_op)) {
    output <- end_op(output)
  }
  if (!is.null(paren)) {
    output <- paren_format(output)
  }
  if (!names) {
    output <- output %>% setNames(1:ncol(output))
  }
  return(output)
}

# Computes value and CIs for an arbitrary specified formula by race/ethnicity
get_prop_bs <- function(f, des, paren = FALSE, end_op = NULL, names = FALSE) {
  y <- eval(parse(text=f))
  overall <- qbinom(c(0.5, 0.025, 0.975), length(y), mean(y))/length(y)
  by_race_eth <- des$variables %>% data.frame(y=y) %>% group_by(RIDRETH) %>%
    dplyr::summarize(
      mean=mean(y),
      ci_l=qbinom(0.025, n(), mean(y))/n(),
      ci_u=qbinom(0.975, n(), mean(y))/n(),
      .groups="drop"
    ) %>% data.frame
  output <- rbind(
    c(overall),
    by_race_eth %>% select(-RIDRETH)
  )
  rownames(output) <- c("Overall", by_race_eth$RIDRETH)
  if (!is.null(end_op)) {
    output <- end_op(output)
  }
  if (!is.null(paren)) {
    output <- paren_format(output)
  }
  if (!names) {
    output <- output %>% setNames(1:ncol(output))
  }
  return(output)
}


# Computes value and CIs for an arbitrary specified formula by race/ethnicity
get_mean <- function(f, des, paren = FALSE, end_op = NULL, names = FALSE) {
  f <- as.formula(f)
  overall <- svymean(x = f, design = des, method = "beta", na.rm=TRUE)
  # Generates output table: rows = race/ethnicity, cols = val, ci_l, ci_u
  output <- rbind(
    Overall=c(overall, sqrt(attributes(overall)$var)),
    svyby(formula = f,
          by = ~RIDRETH, design = des,
          FUN = svymean, na.rm = TRUE
    )[,-1]
  ) %>% setNames(c("mean","se"))
  output$ci_l <- output$mean + qnorm(0.025)*output$se
  output$ci_u <- output$mean + qnorm(0.975)*output$se
  output$se <- NULL
  if (!is.null(end_op)) {
    output <- end_op(output)
  }
  if (!is.null(paren)) {
    output <- paren_format(output)
  }
  if (!names) {
    output <- output %>% setNames(1:ncol(output))
  }
  return(output)
}

# Computes value and CIs for an arbitrary specified formula by race/ethnicity
get_prop <- function(f, des, paren = FALSE, end_op = NULL, names = FALSE) {
  f <- as.formula(f)
  overall <- svyciprop(formula = f, design = des,
                       method = "beta", vartype = "ci", na.rm=TRUE)
  # Generates output table: rows = race/ethnicity, cols = val, ci_l, ci_u
  output <- rbind(
    Overall=c(overall, attributes(overall)$ci),
    svyby(formula = f,
          by = ~RIDRETH, design = des,
          FUN = svyciprop, na.rm = TRUE,
          method = "beta", vartype = "ci"
    )[,-1]
  )
  if (!is.null(end_op)) {
    output <- end_op(output)
  }
  if (!is.null(paren)) {
    output <- paren_format(output)
  }
  if (!names) {
    output <- output %>% setNames(1:ncol(output))
  }
  return(output)
}


va_rating <- function(x) {
  cut(x, breaks = c(-1, 0.4, 0.55, 0.7, 0.8, 100),
      labels = c(1, 0.6, 0.3, 0.1, 0.0)) %>%
    as.character %>% as.numeric %>% return
}

va_payment <- function(rating, if_spouse=FALSE) {
  ifelse(if_spouse,
         va_data[as.character(rating),]$Spouse,
         va_data[as.character(rating),]$No_Spouse
  )
}


ss_qualify <- function(row) {
  # row <- df[df$SEQN == SEQN, c("RIDAGEYR","IS_FEMALE","BMXHT","FEV1","FVC")]
  output_names <- c("By Either","By FEV1","By FVC")
  if (row["RIDAGEYR"] < 18 | as.logical(row["MISSING"])) {
    return(c(NA, NA, NA) %>% setNames(output_names))
  }
  thresholds <- ss[ss$Height_higher > row["BMXHT"] &
                     ss$Height_lower < row["BMXHT"] &
                     ss$Age == ifelse(row["RIDAGEYR"] < 20, "18_to_19", "20_and_older") &
                     ss$Sex == ifelse(row["IS_FEMALE"], "Female", "Male"), "Value"]
  return(c(any(c(row["FEV1"], row["FVC"]) < thresholds),
           row["FEV1"] < thresholds[1],
           row["FVC"] < thresholds[2]) %>%
           setNames(output_names))
}

# Read coefficients for NHANES III reference (Hankinson et al.) and convert to list of lists
nhanes_coefs <- read.csv(file = "data/NHANES_III_coefficients.csv", header = TRUE) %>%
  dlply(1, function(a) {
    dlply(a[,-1], 1, function(b) {
      dlply(b[,-1], 1, function(c) {
        dlply(c[,-1], 1, function(d) {
          d[-1]
        })
      })
    })
  })

# Return individual NHANES III reference estimate based on `nhanes_coef`
# get_nhanes_est_ind <- function(PFT, ethnicity, sex, age, height, output) {
#   if (PFT %in% c("FEV1FEV6", "FEV1FVC")) {
#     age_range <- rep("any", length(age))
#   } else {
#     age_range <- ifelse(sex=="Female",
#                         ifelse(age < 18, "under_18", "18_or_older"),
#                         ifelse(age < 20, "under_20", "20_or_older")
#     )
#   }
#   data.frame(ethnicity=ethnicity, sex=sex, age_range=age_range) %>%
#     apply(1, function(row) {
#       nhanes_coefs[[PFT]][[row["ethnicity"]]][[row["sex"]]][[row["age_range"]]] %>% unlist
#     }) %>% t %>% as.data.frame -> coefs
#   if (output == "LLN") {
#     out <- coefs$Intercept_LLN + coefs$Age*age + coefs$Age_squared*age^2 + coefs$Ht_LLN_cm*height^2
#   } else {
#     out <- coefs$Intercept + coefs$Age*age + coefs$Age_squared*age^2 + coefs$Ht_PRD_cm*height^2
#   }
#   if (PFT == "FEV1FVC") {
#     out <- out/100
#   }
#   return(out)
# }

# Wrapper for multiple NHANES III reference estimates from data frame
get_baseline_est <- function(PFT, df, type = "old") {
  nh3_eth <- df$RIDRETH %>%
    replace(df$RIDRETH == "Hispanic", "Mexican_American") %>%
    replace(grepl("Mexican", df$RIDRETH), "Mexican_American") %>%
    replace(grepl("Asian", df$RIDRETH), "Caucasian") %>%
    replace(grepl("Chinese", df$RIDRETH), "Caucasian") %>%
    replace(grepl("Black",df$RIDRETH), "African_American") %>%
    replace(grepl("White", df$RIDRETH), "Caucasian") %>%
    replace(grepl("Other", df$RIDRETH), "Caucasian")
  nh3_eth_factor <- nh3_eth %>%
    factor(levels = c("Caucasian","African_American","Mexican_American"),
           labels=1:3)
  
  gli_eth <- df$RIDRETH %>%
    replace(df$RIDRETH == "Hispanic", "Caucasian") %>%
    replace(grepl("Mexican", df$RIDRETH), "Caucasian") %>%
    replace(grepl("Asian", df$RIDRETH), "Other") %>%
    replace(grepl("Chinese", df$RIDRETH), "Caucasian") %>%
    replace(grepl("Black",df$RIDRETH), "African_American") %>%
    replace(grepl("White", df$RIDRETH), "Caucasian") %>%
    replace(grepl("Other", df$RIDRETH), "Other")
  gli_eth_factor <- gli_eth %>%
    factor(levels = c("Caucasian","African_American","NE Asian","SE Asian","Other"),
           labels=1:5)
  
  sex <- ifelse(df$IS_FEMALE, "Female", "Male")
  gender <- ifelse(df$IS_FEMALE, 2, 1)
  age <- df$RIDAGEYR
  height <- df$BMXHT
  rs_height <- df$BMXHT/100
  if (PFT == "FEV1") {
    zs_GLI = zscore_GLI(age, rs_height, gender, gli_eth_factor, FEV1=df$FEV1)
    zs_GLI_bl = zscore_GLI(age, rs_height, gender, rep(5, nrow(df)), FEV1=df$FEV1)
    zs_GLI_ca = zscore_GLI(age, rs_height, gender, rep(1, nrow(df)), FEV1=df$FEV1)
  }
  if (PFT == "FVC") {
    zs_GLI = zscore_GLI(age, rs_height, gender, gli_eth_factor, FVC=df$FVC)
    zs_GLI_bl = zscore_GLI(age, rs_height, gender, rep(5, nrow(df)), FVC=df$FVC)
    zs_GLI_ca = zscore_GLI(age, rs_height, gender, rep(1, nrow(df)), FVC=df$FVC)
  }
  if (PFT == "FEV1FVC") {
    zs_GLI = zscore_GLI(age, rs_height, gender, gli_eth_factor, FEV1FVC=df$FEV1FVC)
    zs_GLI_bl = zscore_GLI(age, rs_height, gender, rep(5, nrow(df)), FEV1FVC=df$FEV1FVC)
    zs_GLI_ca = zscore_GLI(age, rs_height, gender, rep(1, nrow(df)), FEV1FVC=df$FEV1FVC)
  }
  PRD_NH3 <- pred_NHANES3(age, rs_height, gender, nh3_eth_factor, param=PFT)
  PRD_NH3_bl <- rowSums(data.frame(
    pred_NHANES3(age, rs_height, gender, rep(1, nrow(df)), param=PFT),
    pred_NHANES3(age, rs_height, gender, rep(2, nrow(df)), param=PFT),
    pred_NHANES3(age, rs_height, gender, rep(3, nrow(df)), param=PFT)))/3
  PRD_NH3_ca <- pred_NHANES3(age, rs_height, gender, 1, param=PFT)
  LLN_NH3 <- LLN_NHANES3(age, rs_height, gender, nh3_eth_factor, param=PFT)
  LLN_NH3_bl <- rowSums(data.frame(
    LLN_NHANES3(age, rs_height, gender, rep(1, nrow(df)), param=PFT),
    LLN_NHANES3(age, rs_height, gender, rep(2, nrow(df)), param=PFT),
    LLN_NHANES3(age, rs_height, gender, rep(3, nrow(df)), param=PFT)))/3
  LLN_NH3_ca <- LLN_NHANES3(age, rs_height, gender, 1, param=PFT)
  NH3_sigma <- (PRD_NH3-LLN_NH3)/qnorm(0.95)
  NH3_sigma_bl <- (PRD_NH3_bl-LLN_NH3_bl)/qnorm(0.95)
  NH3_sigma_ca <- (PRD_NH3_ca-LLN_NH3_ca)/qnorm(0.95)
  col_id <- c("PRD_NH3","PRD_NH3_bl","PRD_NH3_ca",
              "LLN_NH3","LLN_NH3_bl","LLN_NH3_ca",
              "zs_NH3","zs_NH3_bl","zs_NH3_ca",
              "PRD_GLI","PRD_GLI_bl","PRD_GLI_ca",
              "LLN_GLI","LLN_GLI_bl","LLN_GLI_ca",
              "zs_GLI","zs_GLI_bl","zs_GLI_ca")
  output <- suppressWarnings(data.frame(
    PRD_NH3 = PRD_NH3,
    PRD_NH3_bl = PRD_NH3_bl,
    PRD_NH3_ca = PRD_NH3_ca,
    LLN_NH3 = LLN_NH3,
    LLN_NH3_bl = LLN_NH3_bl,
    LLN_NH3_ca = LLN_NH3_ca,
    zs_NH3 = (df[,PFT]-PRD_NH3)/NH3_sigma,
    zs_NH3_bl = (df[,PFT]-PRD_NH3_bl)/NH3_sigma_bl,
    zs_NH3_ca = (df[,PFT]-PRD_NH3_ca)/NH3_sigma_ca,
    PRD_GLI = pred_GLI(age, rs_height, gender, gli_eth_factor, param=PFT),
    PRD_GLI_bl = pred_GLI(age, rs_height, gender, rep(5, nrow(df)), param=PFT),
    PRD_GLI_ca = pred_GLI(age, rs_height, gender, rep(1, nrow(df)), param=PFT),
    LLN_GLI = LLN_GLI(age, rs_height, gender, gli_eth_factor, param=PFT),
    LLN_GLI_bl = LLN_GLI(age, rs_height, gender, rep(5, nrow(df)), param=PFT),
    LLN_GLI_ca = LLN_GLI(age, rs_height, gender, rep(1, nrow(df)), param=PFT),
    zs_GLI = zs_GLI,
    zs_GLI_bl = zs_GLI_bl,
    zs_GLI_ca = zs_GLI_ca
    # PRD = get_nhanes_est_ind(PFT, nh3_eth, sex, age, height, output="PRD"),
    # PRD_bl = get_nhanes_est_ind(PFT, rep("Blended", nrow(df)), sex, age, height, "PRD"),
    # LLN = get_nhanes_est_ind(PFT, nh3_eth, sex, age, height, "LLN"),
    # LLN_bl = get_nhanes_est_ind(PFT, rep("Blended", nrow(df)), sex, age, height, "LLN")
  )) %>% setNames(sprintf("%s_%s", PFT, col_id))
  if (PFT != "FEV1FVC") {
    which_asian <- grepl("Asian", df$RIDRETH) | grepl("Chinese", df$RIDRETH)
    which_NH3 <- grepl("NH3", colnames(output)) & !grepl("NH3_bl", colnames(output))
    output[which_asian, which_NH3] <- 0.88 * output[which_asian, which_NH3]
  }
  return(output)
}

get_var_groups <- function(key) {
  sapply(1:nrow(key), function(x) {
    rep(x, length(eval(parse(text=gsub("-", ":", key$Position[x])))))
  }) %>% unlist %>% return()
}