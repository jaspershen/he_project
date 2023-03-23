no_function()

rm(list = ls())

library(tidymass)
library(tidyverse)

masstools::setwd_project()

load("data_analysis/data_preparation/metabolomics/metabolomics_data")

dir.create("data_analysis/metabolomics/marker", recursive = TRUE)

setwd("data_analysis/metabolomics/marker/")

metabolomics_data@variable_info$Tags

####A(sgLuc, control) VS B(sgPSTK-1, sample)
unique(metabolomics_data@sample_info$group)

control_sample_id <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "sgLuc") %>%
  pull(sample_id)

case_sample_id <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "sgPSTK-1") %>%
  pull(sample_id)

temp_data1 <-
  metabolomics_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  log(2) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    method = "t.test",
    p_adjust_methods = "fdr",
    return_mass_dataset = TRUE
  )

plot <-
  volcano_plot(
    object = temp_data1,
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_column_name = "p_value",
    labs_y = "-log(p-value, 10)",
    add_text = TRUE,
    text_from = "Compound.name",
    point_size_scale = "p_value"
  ) +
  scale_size_continuous(range = c(0.1, 3))

ggsave(plot,
       file = "volcano plot sgLuc vs sgPSTK-1.pdf",
       width = 8,
       height = 7)

marker <-
  temp_data1 %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value < 0.05) %>%
  extract_variable_info()

write.csv(marker, file = "marker sgLuc vs sgPSTK-1.csv", row.names = FALSE)





####A(sgLuc, control) VS B(sgPSTK-1, sample)
unique(metabolomics_data@sample_info$group)

control_sample_id <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "sgLuc") %>%
  pull(sample_id)

case_sample_id <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "sgPSTK-2") %>%
  pull(sample_id)

temp_data1 <-
  metabolomics_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  log(2) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    method = "t.test",
    p_adjust_methods = "fdr",
    return_mass_dataset = TRUE
  )

plot <-
  volcano_plot(
    object = temp_data1,
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_column_name = "p_value",
    labs_y = "-log(p-value, 10)",
    add_text = TRUE,
    text_from = "Compound.name",
    point_size_scale = "p_value"
  ) +
  scale_size_continuous(range = c(0.1, 3))
plot
ggsave(plot,
       file = "volcano plot sgLuc vs sgPSTK-2.pdf",
       width = 8,
       height = 7)

marker <-
  temp_data1 %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value < 0.05) %>%
  extract_variable_info()

write.csv(marker, file = "marker sgLuc vs sgPSTK-2.csv", row.names = FALSE)









####A(sgLuc, control) VS B(sgPSTK-1, sample)
unique(metabolomics_data@sample_info$group)

control_sample_id <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "sgLuc1") %>%
  pull(sample_id)

case_sample_id <-
  metabolomics_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(group == "sgSARS1") %>%
  pull(sample_id)

temp_data1 <-
  metabolomics_data %>%
  mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  log(2) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    method = "t.test",
    p_adjust_methods = "fdr",
    return_mass_dataset = TRUE
  )

plot <-
  volcano_plot(
    object = temp_data1,
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_column_name = "p_value",
    labs_y = "-log(p-value, 10)",
    add_text = TRUE,
    text_from = "Compound.name",
    point_size_scale = "p_value"
  ) +
  scale_size_continuous(range = c(0.1, 3))
plot
ggsave(plot,
       file = "volcano plot sgLuc1 vs sgSARS1.pdf",
       width = 8,
       height = 7)

marker <-
  temp_data1 %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value < 0.05) %>%
  extract_variable_info()

write.csv(marker, file = "marker sgLuc1 vs sgSARS1.csv", row.names = FALSE)
