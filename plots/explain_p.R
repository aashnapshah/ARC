library("dplyr")
library("boot")
library('tidyverse')
library('reticulate')
library('ggplot2')
library('grid')
library('sysfonts')
library('showtext')

library(reticulate)
source_python(normalizePath("../data/dicts/col_dict.py"))
titles_dict <- py_to_r(titles_dict)


font_add_google("Work Sans", "Work Sans")
font = 'Work Sans'
font_size = 20
showtext_auto()

custom_colors <- c('Asian' = '#FFD166', 'Black' = '#FF6F61', 'White' = '#A6D854', 'Hispanic' = '#C94CC2', 'Other' = '#118AB2')
custom_shapes <- c('+' = 4, '-' = 16)

# Function to combine UKB and NHANES cohorts for race-specific plots
combine_cohorts <- function(df1, df2, num=10, ch1='NHANES', ch2='UKB'){
  df1$dataset = ch1
  df2$dataset = ch2
  df = rbind(df1, df2)

  race_df <- df %>%
    filter(startsWith(index, 'race_')) %>%
    filter(!grepl('_adj', index)) %>%
    select(group, cov, index, coef, std_err, dataset, target)  %>%
    filter(!(startsWith(index, 'race_3') & dataset == 'UKB')) %>%
    filter(!(startsWith(index, 'race_2') & dataset == 'NHANES-III'))

  race_df <- race_df %>%
    arrange(index, dataset) %>%
    group_by(index, target, dataset)

  race_df <- race_df %>%
    mutate(fraction = 100*(coef[1] - coef) / coef[1]) %>% # baseline in coef[1]
    ungroup() 

  race_df$dataset <- factor(race_df$dataset, levels = c(ch1, ch2))

  race_df <- race_df %>%
    group_by(index, dataset) %>%
    distinct(cov, .keep_all = TRUE) %>%
    slice_head(n=num) %>%
    ungroup()

  race_df <- race_df %>%
    mutate(
      cov_ds = paste0(cov, " - ", dataset),
      index_ds = paste0(index, " - ", dataset)
    ) %>%
    group_by(dataset) %>%
    mutate(
      cov_ds = factor(cov_ds, levels = unique(cov_ds)),
      index_ds = factor(index_ds, levels = unique(index_ds))
    ) %>%
    ungroup()
  
  print(race_df)
  return(race_df)
}

plot_race_adj <- function(df, pft, fraction=FALSE) {
  ch_levels <- unique(as.character(df$dataset))
  df$cov_ds <- titles_dict[df$cov]
  ch1 <- ch_levels[1]
  ch2 <- ch_levels[2]
  if (fraction == TRUE) {
    y <- "fraction"
    y_label <- paste0("Race-Adjusted Difference Explained in ", pft, " (%)")
    p <- ggplot(df, aes(x = cov_ds, y = fraction,  color = index, group = index))
  } else {
    y <- "coef"
    y_label <- paste0("Race-Adjusted Difference Explained in ", pft, " (L)")
    p <- ggplot(df, aes(x = cov_ds, y = coef, color = index, group = index)) + 
        geom_errorbar(aes(ymin = coef - std_err, ymax = coef + std_err), width = 0.5, alpha = 0.3, linewidth = 0.5) 
  } 
  p <- p +
    geom_point(size = 2, alpha = 1) +
    geom_line(linewidth = 0.6) + 
    facet_wrap(dataset ~ ., nrow = 2, scales = "free_x", strip.position = "right") +
    scale_shape_manual(values = custom_shapes) +
    #scale_color_manual(values = custom_colors) +
    labs(x = "Feature", y = y_label, color = "Race") +
    theme_minimal(base_family = font, base_size = font_size+2) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0.95, size = font_size, color = 'black',
                                margin = margin(t = 0, r = 0, b = 0, l = 0)),
      axis.text.y = element_text(hjust = 1, size = font_size, color = 'black'),
      legend.position = "top",
      legend.box.spacing = unit(0, "pt"),
      legend.margin = margin(0,0,0,0),
      legend.title = element_blank(),
      legend.text = element_text(size = font_size, family = font, color = 'black'),
      panel.border = element_rect(color = "gray", fill = NA, linewidth = 1.0), # updated
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray", linewidth = 0.05),        # updated
      panel.spacing = unit(1, "lines"),
      strip.text = element_text(size = font_size)
    )

  p <- p + scale_x_discrete(labels = function(x) gsub(paste0(" - ", ch1, "| - ", ch2), "", x))
  return(p)
}

# Helper: save a ggplot to PDF at 300 dpi, ensuring directory exists
save_plot_pdf <- function(plot, out_path, width = 12, height = 8, dpi = 300) {
  out_dir <- dirname(out_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  ggplot2::ggsave(
    filename = out_path,
    plot = plot,
    width = width,
    height = height,
    units = 'in',
    dpi = dpi,
    device = grDevices::cairo_pdf
  )
  print(paste0('Saved plot to ', out_path))
}

data_dir = '../results/explain/tables/'
plot_dir = '../results/explain/figures/'

df_ukb = read.csv(paste0(data_dir, 'ukb_cont.csv'))
df_nh3 = read.csv(paste0(data_dir, 'nh3_cont.csv'))
df_nh4 = read.csv(paste0(data_dir, 'nh4_cont.csv'))
df_nh = read.csv(paste0(data_dir, 'nh_cont.csv'))

df = combine_cohorts(df_ukb, df_nh, num = 10, ch1='UKB', ch2='NHANES')
df2 = combine_cohorts(df_nh3, df_nh4, num = 10, ch1='NHANES-III', ch2='NHANES-IV')

targets = c('fev1', 'fvc')
for (i in seq_along(list(df, df2))) {

  df <- list(df, df2)[[i]]
  name <- c('main', 'supplemental')[i]
  for (target in targets) {
    subset_df <- df %>% filter(target == target)
    p_l <- plot_race_adj(subset_df, target, fraction=FALSE)
    p_pct <- plot_race_adj(subset_df, target, fraction=TRUE)
    save_plot_pdf(p_l,   file.path(plot_dir, paste0('explain_mL_', target, '_', name, '.pdf')), width = 12, height = 8, dpi = 300)
    save_plot_pdf(p_pct, file.path(plot_dir, paste0('explain_pct_', target, '_', name, '.pdf')), width = 12, height = 8, dpi = 300)
  } 
}
