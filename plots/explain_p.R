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
font_size = 20
showtext_auto()

setwd("/n/groups/patel/aashna/PFT-ML-2024/")

# combine datasets
col_desc_dict <- read.csv("data/predictor_category_map.csv")
col_cat_dict <- read.csv("data/column_description_dictionary.csv")

# Define a custom color palette for race
custom_colors <- c('Asian' = '#FFD166', 'Black' = '#FF6F61', 'White' = '#A6D854', 'Hispanic' = '#C94CC2', 'Other' = '#118AB2')
custom_shapes <- c('+' = 4, '-' = 16)

model_labels <- c(
  'Baseline' = 'Age, Sex',
  '+ BMXHT' = '+ Height',
  '+ BMXWAIST' = '+ Waist Circumference',
  '+ BMXSIT' = '+ Sitting Height',
  '+ BMXHIP' = '+ Hip Circumference',
  '+ BMXARML' = '+ Arm Length',
  '+ DMD_IncPovRatio' = '+ Income to Poverty Ratio',
  '+ BMXARMC' = '+ Arm Circumference',
  '+ BMXBMI' = '+ Body Mass Index',
  '+ DMD_English_Lang' = '+ English Language',
  '+ DMD_Born_US' = '+ Native Born',
  '+ BMXLEG' = '+ Leg Length',
  '+ DMD_Health_Ins' = '+ Health Insurance',
  '+ DMD_Veteran' = '+ Military Veteran',
  '+ BMXWT' = '+ Body Weight',
  '+ DMD_PM2.5_10' = '+ PM2.5 Exposure',
  '+ DMD_Born_UK' = '+ Native Born',
  '+ DMD_Military_HC' = '+ Military Healthcare',
  '+ DMD_Married' = '+ Marital Status',
  '+ DMD_Home_Smoke' = '+ Home Smoking Exposure',
  '+ DMD_Finish_HS' = '+ High School Completion',
  '+ DMD_HH_Size' = '+ Household Size',
  '+ DMD_Work_Smoke' = '+ Work Smoking Exposure',
  '+ DMD_INC' = '+ Income',
  '+ OCQ_Work_Smoke' = '+ Work Smoking Exposure'
)

combine_cohorts <- function(df_ukb, df_nhanes, num=10){
  
  # add a new column specificying the dataset (ukb or nhanes) and then combine them
  df_ukb$dataset = 'UKB'
  df_nhanes$dataset = 'NHANES'
  # combine datasets
  df = rbind(df_ukb, df_nhanes)
  df$Covariate <- gsub("_spline", "", df$Covariate)
  df$Model <- gsub("_spline", "", df$Model)
  
  ridreth_df <- df %>%
    filter(startsWith(Covariate, 'RIDRETH')) %>%
    select(Covariate, Model, Estimate, `Std. Error`, dataset, MSE)
  
  # Remove "RIDRETH_" from Covariate names
  ridreth_df$Covariate <- gsub('RIDRETH_', '', ridreth_df$Covariate)
  ridreth_df$Model_ <- model_labels[ridreth_df$Model]
  ridreth_df$index <- paste(ridreth_df$Model_, ridreth_df$dataset, sep = " - ")
  ridreth_df$index <- factor(ridreth_df$index, levels = unique(ridreth_df$index))
  
  #group by covariate and get fraction explained after addition of new estimate then ungroup
  
  ridreth_df <- ridreth_df %>%
    group_by(Covariate, dataset) %>%
    mutate(Fraction = (Estimate[2] - Estimate) / Estimate[2]) %>%
    ungroup()
  
  # group by dataset and only keep the first 12 uniwue models
  ridreth_df$dataset <- factor(ridreth_df$dataset, levels = c('UKB', 'NHANES'))
  ridreth_df <- ridreth_df %>%
    group_by(Covariate, dataset) %>%
    distinct(Model_, .keep_all = TRUE) %>%
    slice_head(n=num) %>%
    ungroup()
  

  return(ridreth_df)
}


# Define a custom color palette for race
custom_colors <- c('Asian' = '#FFD166', 'Black' = '#FF6F61', 'White' = '#A6D854', 'Hispanic' = '#C94CC2', 'Other' = '#118AB2')
custom_shapes <- c('+' = 4, '-' = 16)

plot_race_adj <- function(df, pft) {
  # Create the plot
  p <- ggplot(df, aes(x = index, y = Estimate, ymin = Estimate - `Std. Error`, ymax = Estimate + `Std. Error`, 
                              color = Covariate, group = Covariate)) +
    geom_point(size = 2, alpha = 1) +
    geom_errorbar(width = 0.5, alpha = 0.3) +
    geom_line() +
    facet_wrap(dataset ~ ., nrow = 2, scales = "free_x", strip.position = "right") +
    scale_shape_manual(values = custom_shapes) +
    scale_color_manual(values = custom_colors) +
    labs(x = "Covariate", y = paste("Race-Adjusted Difference Explained in",pft,"(L)"), color = "Race") +
    theme_minimal(base_family = font, base_size = font_size+2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.95, size = font_size, color = 'black', 
                                     margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(hjust = 1, size = font_size, color = 'black'),
          legend.position = "top",
          legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
          legend.margin=margin(0,0,0,0), #rtical spacing between legend and plot panels
          legend.title = element_blank(), #element_text(size = font_size, family = font, color = 'black'),  # Centered title
          legend.text = element_text(size = font_size, family = font, color = 'black'),
          panel.border = element_rect(color = "gray", fill = NA, size = 1.0),
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.grid.major = element_line(color = "gray", size = 0.05),
          panel.spacing = unit(1, "lines"),  # Increase space between panels
          strip.text = element_text(size = font_size))
  
  
  p <- p + scale_x_discrete(labels = function(x) gsub(" - NHANES| - UKB", "", x))
  return(p)
} 


plot_race_adj_frac <- function(df, pft) {
  # Create the plot
  p <- ggplot(df, aes(x = index, y = Fraction*100,
                      color = Covariate, group = Covariate)) +
    geom_point(size = 2, alpha = 1) +
    geom_line() +
    facet_wrap(dataset ~ ., nrow = 2, scales = "free_x", strip.position = "right") +
    scale_shape_manual(values = custom_shapes) +
    scale_color_manual(values = custom_colors) +
    labs(x = "Covariate", y = paste("Race-Adjusted Difference Explained in", pft, "(%)"), color = "Race") +
    theme_minimal(base_family = font, base_size = font_size+2) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 0.95, size = font_size, color = 'black',
                        margin = margin(t = 0, r = 0, b = 0, l = 0)),
          axis.text.y = element_text(hjust = 1, size = font_size, color = 'black'),
          legend.position = "top",
          legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
          legend.margin=margin(0,0,0,0), #rtical spacing between legend and plot panels
          legend.title = element_blank(), #element_text(size = font_size, family = font, color = 'black'),  # Centered title
          legend.text = element_text(size = font_size, family = font, color = 'black'),
          panel.border = element_rect(color = "gray", fill = NA, size = 1.0),
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.grid.major = element_line(color = "gray", size = 0.05),
          panel.spacing = unit(1, "lines"),  # Increase space between panels
          strip.text = element_text(size = font_size))
  
  
  p <- p + scale_x_discrete(labels = function(x) gsub(" - NHANES| - UKB", "", x))
  return(p)
} 

