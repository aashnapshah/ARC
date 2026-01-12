library("dplyr")
library("boot")
library('tidyverse')
library('reticulate')
library('ggplot2')
library('sysfonts')
library('showtext')

font_add_google("Open Sans", "Open Sans")
font_add_google("Work Sans", "Work Sans")

font = 'Work Sans'
font_size = 18
showtext_auto()

setwd("~/Desktop/ssh_mount/ARC")

# combine datasets
col_desc_dict <- read.csv("data/predictor_category_map.csv")
col_cat_dict <- read.csv("data/column_description_dictionary.csv")
custom_colors <- c('Anthropometrics' = '#3357FF', 'Race' = '#FF5733', 
                   'Sociodemographics' = '#A6D854', 'Exposures' = '#C94CC2', 'Other' = '#118AB2')

combine_cohorts <- function(ukb_data, nhanes_data) {
  df = rbind(ukb_data, nhanes_data)
  # Merge df with col_cat_dict based on Model
  df <- merge(df, col_desc_dict, by.x = "Model", by.y = "Predictor.Name", all.x = TRUE) %>%
    filter(Category %in% c('Anthropometrics', 'Exposures', 'Sociodemographics') | Model == 'RIDRETH') %>%
    mutate(Category = ifelse(is.na(Category), "Race" , Category)) %>% 
    filter(startsWith(Covariate, Model)) 
  
  #If covariate starts with ridreth_, extract what is after else, keep the same
  df$Covariate <- gsub("RIDRETH_", "", df$Covariate) 
  df$Column.Name <- ifelse(is.na(df$Column.Name), df$Covariate, df$Column.Name)
  df$dataset <- factor(df$dataset, levels = c('UKB', 'NHANES'))
  return(df)
}

plot_p_values <- function(df) {
  alpha <- 0.05
  num_tests <- nrow(df)  # Assuming each row represents a test
  bonferroni_threshold <- alpha / num_tests
  
  # Convert p-value to numeric (assuming it's in a column named 'Pr(>|t|)')
  df$P.value <- as.numeric(df[['Pr(>|t|)']])
  df$log_P_value <- -log10(df$P.value)
  
  df_sorted <- df[order(df$dataset, df$Category, -df$log_P_value), ]
  df_sorted$Category <- factor(df_sorted$Category, levels = c('Race', 'Anthropometrics', 
                                                              'Sociodemographics', 'Exposures'))
  df_sorted <- df_sorted[order(df_sorted$Category), ]
  df_sorted$Column.Name <- factor(df_sorted$Column.Name, levels = unique(df_sorted$Column.Name))
  df_sorted$log_P_value[is.infinite(df_sorted$log_P_value)] <- 400
  
  custom_shapes <- c('+' = 4, '-' = 16)
  
  # Create the plot
  p <- ggplot(df_sorted, aes(x = Column.Name, y = log_P_value, color = Category)) +
    geom_point(size = 2, alpha = 0.8) +
    # make scales free
    facet_grid(dataset ~ ., scales = 'free') +
    geom_hline(yintercept = -log10(bonferroni_threshold), linetype = "dashed", color = "red") +
    labs(x = "", y = "-log10(P-value)", color = "Category", fill = "Category") +
    scale_color_manual(values = custom_colors) +
    scale_fill_manual(values = custom_colors) +
    scale_shape_manual(values = custom_shapes) +
    coord_cartesian(ylim = c(0, 400)) +
    theme_minimal(base_family = font, base_size = font_size+2) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = font_size, color = 'black'),
          axis.text.y = element_text(hjust = 1, size = font_size, color = 'black'),
          legend.position = "bottom",
          legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
          legend.margin=margin(0,0,0,0), #rtical spacing between legend and plot panels
          legend.title = element_blank(), #element_text(size = font_size, family = font, color = 'black'),  # Centered title
          legend.text = element_text(size = font_size, family = font, color = 'black'),
          panel.border = element_rect(color = "gray", fill = NA, size = 1.0),
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.grid.major = element_line(color = "gray", size = 0.05),
          panel.spacing = unit(1, "lines"),  # Increase space between panels
          strip.text = element_text(size = font_size))
  
  return(p)
}

plot_effect_size <- function(df) {
  # Assuming you have already loaded the necessary libraries and have the required data in the df data frame
  alpha <- 0.05
  num_comparisons <- nrow(df)  # Assuming each row represents a comparison
  bonferroni_threshold <- alpha / num_comparisons
  df$P.value <- as.numeric(df[['Pr(>|t|)']])
  
  # Filter rows passing Bonferroni test
  filtered_df <- df[df$P.value < bonferroni_threshold, ]
  filtered_df$Category <- factor(filtered_df$Category, levels = c('Race', 'Anthropometrics', 'Sociodemographics',
                                                                  'Exposures'))
  filtered_df <- filtered_df[order(filtered_df$Category, filtered_df$Estimate), ]
  filtered_df$Column.Name <- factor(filtered_df$Column.Name, levels = unique(filtered_df$Column.Name))
  # Flip the plot and add standard error bars
  p <- ggplot(filtered_df, aes(x = Estimate, y = Column.Name, color = dataset)) +
    geom_errorbar(aes(xmin = Estimate - `Std. Error`, xmax = Estimate + `Std. Error`), width = 0.5, 
                  alpha = 0.3, position = position_dodge(width = 0.0)) + 
    geom_point(position = position_dodge(width = 0.0), size = 1.5, shape = 16, fill = "white") +  # Use filled circle with white color
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.7, alpha = 0.7) +  # Add line at y = 0 within the specified ylim
    facet_wrap(~ Category, scales = "free", ncol = 2) +  
    labs(x = "FVC Effect Size (mL)", y = "Variable", color = "Dataset") +
    scale_color_manual(values = c("#C94CC2", "#118AB2")) +  
    theme_minimal(base_family = font, base_size = font_size + 2) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = font_size + 3, color = 'black'),
          axis.text.y = element_text(hjust = 1, size = font_size + 3, color = 'black'),
          legend.position = "top",
          legend.title = element_blank(), #element_text(size = font_size, family = font, color = 'black'),  # Centered title
          legend.text = element_text(size = font_size+10, family = font, color = 'black'),
          legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
          legend.margin=margin(0,0,0,0), # the margin around each legendst vertical spacing between legend and plot panels
          panel.border = element_rect(color = "gray", fill = NA, size = 1.0),
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.grid.major = element_line(color = "gray", size = 0.05),
          panel.spacing = unit(1, "lines"),  # Increase space between panels
          strip.text = element_text(size = font_size+7))
  return(p)
}

