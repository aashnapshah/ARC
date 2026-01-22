# Can run this R script with:
# $ Rscript ewas_p.R

library("dplyr")
library('tidyverse')
library('reticulate')
library('ggplot2')
library('sysfonts')
library('showtext')

plot_p_values <- function(df, pdf_path=NULL, width=10, height=7) {
  alpha <- 0.05
  num_tests <- nrow(df)
  bonferroni_threshold <- alpha / num_tests
  df$p.value <- as.numeric(df$p.value)
  df$log_p.value <- -log10(df$p.value)
  df$log_p.value[is.infinite(df$log_p.value)] <- max(df$log_p.value[is.finite(df$log_p.value)], na.rm=TRUE) + 2

  # Capitalize for nice facet order
  df$Category <- tools::toTitleCase(df$Category)
  df$Variable <- factor(df$Variable, levels = unique(df$Variable[order(df$Category, df$log_P_value)]))
  cat_levels <- unique(df$Category)
  if ("Race" %in% cat_levels) cat_levels <- c("Race", setdiff(cat_levels, "Race"))
  df$Category <- factor(df$Category, levels=cat_levels)
  df$dataset <- factor(df$dataset)

  # Plot
  p <- ggplot(df, aes(x = Variable, y = log_P_value, color = Category)) +
    geom_point(size=2, alpha=0.8) +
    facet_grid(dataset ~ ., scales = 'free_y') +
    geom_hline(yintercept = -log10(bonferroni_threshold), linetype="dashed", color="red") +
    labs(x = "", y = "-log10(P-value)", color="Category") +
    scale_color_manual(values = custom_colors, aesthetics = c("color", "fill")) +
    coord_cartesian(ylim = c(0, max(df$log_P_value, na.rm=TRUE) + 1)) +
    theme_minimal(base_family=font, base_size=font_size+2) +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=font_size, color='black'),
          axis.text.y = element_text(hjust=1, size=font_size, color='black'),
          legend.position = "bottom",
          legend.box.spacing = unit(0, "pt"),
          legend.margin=margin(0,0,0,0),
          legend.title = element_blank(),
          legend.text = element_text(size = font_size, family = font, color = 'black'),
          panel.border = element_rect(color="gray", fill=NA, size=1.0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color="gray", size=0.05),
          panel.spacing = unit(1, "lines"),
          strip.text = element_text(size=font_size))
  # Save as PDF if requested
  if (!is.null(pdf_path)) {
    ggsave(pdf_path, plot=p, width=width, height=height, device=cairo_pdf)
  }
  return(p)
}

plot_effect_size <- function(df, pdf_path=NULL, width=10, height=7) {
  alpha <- 0.05
  num_comparisons <- nrow(df)
  bonferroni_threshold <- alpha / num_comparisons
  p_col <- grep("^P($|\\.|r|>[|]t[|]$)", names(df), value=TRUE, ignore.case=TRUE)[1]
  df$P.value <- as.numeric(df[[p_col]])
  filtered_df <- df[df$P.value < bonferroni_threshold, ]
  filtered_df$Category <- tools::toTitleCase(filtered_df$Category)
  filtered_df$Variable <- factor(filtered_df$Variable, levels = unique(filtered_df$Variable[order(filtered_df$Category, filtered_df$Estimate)]))
  cat_levels <- unique(filtered_df$Category)
  if ("Race" %in% cat_levels) cat_levels <- c("Race", setdiff(cat_levels, "Race"))
  filtered_df$Category <- factor(filtered_df$Category, levels=cat_levels)
  filtered_df$dataset <- factor(filtered_df$dataset)
  p <- ggplot(filtered_df, aes(x=Estimate, y=Variable, color=dataset)) +
    geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`),
                  width=0.5, alpha = 0.3, position = position_dodge(width = 0.0)) +
    geom_point(position = position_dodge(width=0.0), size=1.5, shape=16, fill="white") +
    geom_vline(xintercept=0, linetype="dashed", color="red", size=0.7, alpha=0.7) +
    facet_wrap(~ Category, scales="free", ncol=2) +
    labs(x="FEV1 Effect Size (mL)", y="Variable", color="") +
    scale_color_manual(values = c("UKB"="#C94CC2", "NHANES"="#118AB2")) +
    theme_minimal(base_family=font, base_size=font_size+2) +
    theme(axis.text.x = element_text(angle=0, hjust=1, size=font_size+3, color='black'),
          axis.text.y = element_text(hjust=1, size=font_size+3, color='black'),
          legend.position = "top",
          legend.text = element_text(size = font_size+10, family = font, color = 'black'),
          legend.box.spacing = unit(0, "pt"),
          legend.margin = margin(0,0,0,0),
          panel.border = element_rect(color = "gray", fill = NA, size = 1.0),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "gray", size = 0.05),
          panel.spacing = unit(1, "lines"),
          strip.text = element_text(size = font_size+7))
  if (!is.null(pdf_path)) {
    ggsave(pdf_path, plot=p, width=width, height=height, device=cairo_pdf)
  }
  return(p)
}

font_add_google("Work Sans", "Work Sans")
font = 'Work Sans'
font_size = 18
showtext_auto()

dir = '../results/ewas/tables/'
cohorts = c('nh', 'ukb')
save_dir = '../results/ewas/plots/'

custom_colors <- c(
  'anthropometrics' = '#3357FF',
  'race' = '#FF5733',
  'demographics' = '#A6D854',
  'exposures' = '#C94CC2',
  'other' = '#118AB2'
)

capitalize_words <- function(x) {
  if (is.null(x)) return(x)
  sapply(x, function(s) {
    s <- gsub("_", " ", as.character(s))
    if (nchar(s) == 0) return(s)
    if (toupper(s) == s) return(s)
    tools::toTitleCase(tolower(s))
  }, USE.NAMES = FALSE)
}

standardize_dataset_label <- function(cohort) {
  lc <- tolower(cohort)
  if (startsWith(lc, 'nh') || lc %in% c('nhanes')) {
    return('NHANES')
  }
  if (lc %in% c('uk', 'ukb', 'ukbiobank')) {
    return('UKB')
  }
  return(toupper(cohort))
}

read_and_standardize <- function(file_path) {
  df <- read.csv(file_path)
  # Standardize column names if needed for compatibility
  colnames(df) <- gsub("Coef\\.|Estimate", "Estimate", colnames(df))
  colnames(df) <- gsub("Std\\.Err\\.|Std. Error", "Std. Error", colnames(df))
  # Add required columns if missing for backward compatibility
  if (!("Category" %in% names(df)) & ("category" %in% names(df))) df$Category <- df$category
  # Prefer 'val' if present; otherwise fall back to 'cov' then 'covariate'
  if ("val" %in% names(df)) {
    df$Variable <- df$val
  } else if ("cov" %in% names(df)) {
    df$Variable <- df$cov
  } else if (!("Variable" %in% names(df)) & ("covariate" %in% names(df))) {
    df$Variable <- df$covariate
  }
  # Capitalize labels for readability
  if ("Variable" %in% names(df)) df$Variable <- capitalize_words(df$Variable)
  if ("Category" %in% names(df)) df$Category <- capitalize_words(df$Category)
  return(df)
}

load_selected_cohorts <- function(selected_cohorts) {
  dfs <- list()
  for (cohort in selected_cohorts) {
    file_path <- file.path(dir, paste0(cohort, "_ref_ewas.csv"))
    if (file.exists(file_path)) {
      message("Loading: ", file_path)
      df <- read_and_standardize(file_path)
      df$dataset <- standardize_dataset_label(cohort)
      dfs[[length(dfs)+1]] <- df
    } else {
      message("No file found for cohort: ", cohort, " (looked for ", file_path, ")")
    }
  }
  if (length(dfs) == 0) {
    return(data.frame())
  }
  dplyr::bind_rows(dfs)
}
args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 0) {
  cohorts <- tolower(args)
  message("Cohorts specified via CLI: ", paste(cohorts, collapse = ", "))
}

combined_df <- load_selected_cohorts(cohorts)
if (nrow(combined_df) > 0) {
  pdf_out_e_comb <- file.path(save_dir, paste0("combined_effectsizes_", paste(toupper(cohorts), collapse="_"), ".pdf"))
  plot_effect_size(combined_df, pdf_path=pdf_out_e_comb, width=12, height=8)
  message("Saved combined effect-size plot: ", pdf_out_e_comb)

  pdf_out_p_comb <- file.path(save_dir, paste0("combined_pvalues_", paste(toupper(cohorts), collapse="_"), ".pdf"))
  plot_p_values(combined_df, pdf_path=pdf_out_p_comb, width=12, height=8)
  message("Saved combined p-values plot: ", pdf_out_p_comb)
} else {
  message("No data loaded. Exiting without plotting.")
}