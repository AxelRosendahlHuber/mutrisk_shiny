# Variables needed for the mutrisk package
library(data.table)
library(scales)

# default number of cells - 70kg male individual
tissue_ncells = data.frame(tissue = c("blood", "colon", "lung", "skin"),
                           ncells = c(1e5, 6.60e7, 4.33e+9, 2.53e+9))

# default number of cells - low: female, high: male, mid: middle
tissue_ncells_ci = data.frame(tissue = c("blood", "colon", "lung", "skin"),
                              high_estimate = c(1.3e6, 6.60e7, 4.33e+9, 2.53e+9),
                              mid_estimate = c(1e5, NA, NA, NA),
                              low_estimate = c(2.5e4, 6.42e+7, 3.87e+9, 1.76e+9))
# for all values for which we do no have the 'mean
tissue_ncells_ci$mid_estimate[2:4] = (tissue_ncells_ci$high_estimate[2:4] + tissue_ncells_ci$low_estimate[2:4]) /2
# Take the 'exteme of the values to demonstrate that most estimates still hold
tissue_ncells_ci_wide = tissue_ncells_ci
tissue_ncells_ci_wide$high_estimate[2:4] = tissue_ncells_ci$high_estimate[2:4] * 5
tissue_ncells_ci_wide$low_estimate[2:4] = tissue_ncells_ci$low_estimate[2:4] / 5

# Default colors for the different tissues
blood_colors = c(normal ="#ff725c", chemotherapy = "lightgreen")
lung_colors = c(`non-smoker` = "#4269d0", `ex-smoker` = "#7c86a1", smoker = "#161459")
colon_colors = c(normal = "#3ca951", IBD = "#6cc5b0", POLE = "#145220", POLD1 = "#222e24")
skin_colors = c("#efb118" , "#ff8ab7", "#9c6b4e")
liver_colors = c("#800000", "#B87333", "#5C4033")

tissue_colors = list(blood = blood_colors,
                     lung = lung_colors,
                     colon = colon_colors)

color_df = lapply(tissue_colors, \(x)
                  data.table::data.table(category = names(x), color = x)) |>
  data.table::rbindlist(idcol = "tissue") |>
  dplyr::mutate(tissue_category = paste0(tissue, "_", category))
tissue_category_colors = setNames(color_df$color, color_df$tissue_category)

