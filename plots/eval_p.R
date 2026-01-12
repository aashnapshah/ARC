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

process_df <- function(df) {
  colnames(df) <- colnames(df) <- gsub(" ", ".", colnames(df))
  
  # Extract the lower limit CI value from the string like "(0.499, 0.507)" and get 0.499
  #df$ci_lower <- as.numeric(gsub("\\((.*), .*\\)", "\\1", df$MAE.CI))
  # Extract the upper limit CI value from the string like "(0.499, 0.507)" and get 0.507
  #df$ci_upper <- as.numeric(gsub("\\(.*, (.*)\\)", "\\1", df$MAE.CI))
  
  df <- df %>% 
    filter(Race != "All")  %>%
    filter(Model.Description != "b") %>%
    filter(Model.Description != 'bhwsfu') 
  
  pattern <- "gli|bh.*_"
  
  race_order <- c("Overall", "White", "Black", "Asian", "Hispanic", "Other")
  dataset_order <- c("UKB Test", "NHANES")
  model_order <- c("GLI-Global", "Height", "+ Sitting\nHeight", "+ Waist\nCircumference", "+ Native\nBorn", "+ Smoke\nExposure")
  pft_order <- c("FEV1", "FVC", "FEV1/FVC")
  df <- df %>%
    mutate(model = ifelse(Model.Description == "gli", "GLI-Global", Model.Description))  %>%
    mutate(model = ifelse(model == "bh", "Height", model))  %>%
    mutate(model = ifelse(model == "bhs", "+ Sitting\nHeight", model))  %>%
    mutate(model = ifelse(model == "bhws", "+ Waist\nCircumference", model))  %>%
    mutate(model = ifelse(model == "bhwsu", "+ Native\nBorn", model))  %>%
    #mutate(model = ifelse(model == "bhwsfu", "+ Education", model))  %>%
    mutate(model = ifelse(model == "bhwsfus", "+ Smoke\nExposure", model))
  
  df <- df %>% 
    mutate(race_included = ifelse(Race_Included, "+", "-")) %>%
    mutate(model = factor(model, levels = model_order)) %>%
    mutate(race = factor(Race, levels = race_order)) %>%
    mutate(Dataset = factor(Dataset, levels = dataset_order)) %>%
    #drop null
    drop_na()
  return(df)
}

# Define a custom color palette for race
custom_colors <- c('Asian' = '#FFD166', 'Black' = '#FF6F61', 'White' = '#A6D854', 'Hispanic' = '#C94CC2', 'Other' = '#118AB2')
custom_shapes <- c('+' = 4, '-' = 16)

plot_maes <- function(df) {
  # Create the plot with improvements

  
  p <- ggplot(df, aes(x = race_included, y = Z, color = race, shape = race_included)) +
  geom_line(aes(group = interaction(race, Dataset)), size = 0.5, alpha = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 1, alpha = 0.8) + # Add lines
  geom_errorbar(aes(ymin = MAE_CI_Lower, ymax = MAE_CI_Upper), width = 0.5, alpha = 0.3, position = position_dodge(width = 0.5)) +
  facet_grid(Dataset ~ model, scales = "free_y") +
  scale_shape_manual(values = custom_shapes) +
  scale_color_manual(values = custom_colors) +
  labs(x = "Race Corrected", y = "Mean Absolute Error (L)", color = "Race", shape = 'Race Corrected') +
    theme_minimal(base_family = font, base_size = font_size+2) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.90, size = font_size+10, color = 'black',
                                    margin = margin(0,0,0,0)),
          axis.text.y = element_text(hjust = 1, size = font_size, color = 'black'),
          legend.position = "top",
          legend.box.spacing = unit(0, "pt"),# The spacing between the plotting area and the legend box (unit)
          legend.margin=margin(0,0,0,0), #rtical spacing between legend and plot panels
          legend.title = element_blank(), #element_text(size = font_size, family = font, color = 'black'),  # Centered title
          legend.text = element_text(size = font_size, family = font, color = 'black'),
          panel.border = element_rect(color = "gray", fill = NA, size = 1.0),
          panel.grid.minor = element_blank(),  # Remove minor grid lines
          panel.grid.major = element_line(color = "gray", size = 0.05),
          panel.spacing = unit(0.2, "lines"),  # Increase space between panels
          strip.text = element_text(size = font_size, lineheight = 0.4)) + #, margin = margin(t = 0.5, r = 0.5, b = 0.5, l = 0.5, unit = "pt"))) +
    guides(shape = FALSE)

  return(p)
}



fvc_metrics_data <- read_csv('/n/groups/patel/aashna/PFT-ML-2024/projects/reference_pfts/metrics/combined_mae_results_07.31.24.csv', show_col_types = FALSE) %>% filter(Target == 'FVC') #ukb_train/metrics_mae_gli_FVC_02052024.csv', show_col_types = FALSE)
df_fvc <- process_df(fvc_metrics_data)

fev1_metrics_data <- read_csv('/n/groups/patel/aashna/PFT-ML-2024/projects/reference_pfts/metrics/combined_mae_results_07.31.24.csv', show_col_types = FALSE) %>% filter(Target == 'FEV1')
df_fev1 <- process_df(fev1_metrics_data)

p_fvc <- plot_maes(df_fvc)
p_fev1 <- plot_maes(df_fev1)

print(p_fvc)
print(p_fev1)

base_dir <- '/n/groups/patel/aashna/PFT-ML-2024/'

ggsave(sprintf('%sprojects/reference_pfts/figures/mae_fvc_%s.png', base_dir, Sys.Date()), p_fvc, width = 6, height = 5, dpi = 300, units = 'in')
ggsave(sprintf('%sprojects/reference_pfts/figures/mae_fev1_%s.png', base_dir, Sys.Date()), p_fev1, width = 6, height = 5, dpi = 300, units = 'in')
# 
# for (dataset in c('UKB Test', 'NHANES')) {
#   for (race_ in c(TRUE, FALSE)) {
#     df_ <- df_fev1 %>% filter(Race_Included==race_) %>% filter(Dataset == dataset) 
#     df <- df_ %>%
#       mutate(across(where(is.numeric), round, 3)) %>%
#       mutate(Z_CI_Combined = paste0(Z, " (", MAE_CI_Lower, MAE_CI_Upper, ")")) %>%
#       select(model, race, Z_CI_Combined) %>%
#       pivot_wider(names_from = race, values_from = Z_CI_Combined) %>%
#       arrange(model)  
#     print(dataset)
#     print(race)
# 
#     base_dir <- '/n/groups/patel/aashna/PFT-ML-2024/'
#     file_name <- sprintf('%sprojects/reference_pfts/tables/zscore_fev1_%s_%s_%s.csv', base_dir, dataset, race_, Sys.Date())
#     write_csv(df, file_name)
# 
#   }
# }
# 
# 
# 
# 
