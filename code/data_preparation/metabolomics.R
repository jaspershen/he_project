no_function()

rm(list = ls())

library(tidymass)
library(tidyverse)

masstools::setwd_project()

sample_info <-
  readxl::read_xlsx("raw_data/metabolomics/Metabolites Lingli He 02-23-2023.xlsx")

data <-
  readxl::read_xlsx("raw_data/metabolomics/SUB12576_MetabolomicsData-xiaotao.xlsx")

colnames(data)

data <-
  data[, -c(15:34, 65:67)]


variable_info <-
  data[, c(1:14)]

expression_data <-
  data[, -c(1:14)]

colnames(variable_info)

mz_rt <-
  variable_info[, 1] %>%
  separate(col = variable_id,
           sep = "@",
           into = c("mz", "rt")) %>%
  dplyr::mutate(mz = as.numeric(mz),
                rt = as.numeric(rt) * 60)

variable_info <-
  cbind(variable_info,
        mz_rt) %>%
  dplyr::rename(Compound.name = Name) %>%
  dplyr::mutate(Formula = stringr::str_replace_all(Formula, " ", ""))

expression_data <-
  as.data.frame(expression_data)

rownames(expression_data) <- variable_info$variable_id

colnames(expression_data)

plot(
  expression_data$`Area: SUB12576p3_SPL01.raw (F2)`,
  expression_data$`Norm. Area: SUB12576p3_SPL01.raw (F2)`
)

expression_data <-
  expression_data[, c(16:30)]

colnames(expression_data) <-
  colnames(expression_data) %>%
  stringr::str_replace("Norm. Area: ", "") %>%
  stringr::str_replace("\\.raw", "") %>%
  stringr::str_replace("SUB12576p3_SPL", "") %>%
  stringr::str_replace(" \\(F[0-9]{1,2}\\)", "") %>%
  stringr::str_replace("^0", "") %>%
  paste0("sample_", .)

sample_info <-
  sample_info %>%
  dplyr::mutate(sample_id = paste0("sample_", number)) %>%
  dplyr::select(sample_id, dplyr::everything()) %>%
  dplyr::rename(group = `Gene edited`) %>%
  dplyr::mutate(class = "Subject")

colnames(expression_data) == sample_info$sample_id

metabolomics_data <-
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

dir.create("data_analysis/data_preparation/metabolomics", recursive = TRUE)

metabolomics_data <-
  metabolomics_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(!is.na(Compound.name)) %>%
  dplyr::filter(MS2 != "No MS2") %>%
  dplyr::filter(Tags != "Weak/poor match")

###add HMDB or KEGG ID
name <-
  metabolomics_data@variable_info$Compound.name

name <-
  name %>%
  stringr::str_replace("β-", "") %>%
  stringr::str_replace("Υ-", "") %>%
  stringr::str_replace("Similar to ", "") %>%
  stringr::str_replace(" \\(GABA\\)", "") %>%
  stringr::str_replace(" \\(PFBS\\)", "") %>%
  stringr::str_replace(" \\(UMP\\)", "") %>%
  stringr::str_replace(" \\(UDP\\)", "") %>%
  stringr::str_replace(" \\(NAD\\+\\)", "") %>%
  stringr::str_replace(" \\(GTP\\)", "") %>%
  stringr::str_replace(" \\(IMP\\)", "") %>%
  stringr::str_replace(" \\(GMP\\)", "") %>%
  stringr::str_replace(" \\(GDP\\)", "") %>%
  stringr::str_replace(" \\(FAD\\)", "") %>%
  stringr::str_replace(" \\(ATP\\)", "") %>%
  stringr::str_replace(" \\(dGMP\\)", "") %>%
  stringr::str_replace("FISH candidate: ", "") %>%
  stringr::str_replace("similar to ", "") %>%
  stringr::str_replace("Compound possibly containing ", "") %>%
  stringr::str_replace(" or similar", "") %>%
  stringr::str_replace(" possibly", "")

KEGG.ID <-
  HMDB.ID <-
  vector(mode = "list", length = length(name))

for (i in 1:length(HMDB.ID)) {
  cat(i, " ")
  HMDB.ID[[i]] <-
    masstools::trans_ID(query = name[i],
                        from = "Chemical Name",
                        to = "Human Metabolome Database")
  
  KEGG.ID[[i]] <-
    masstools::trans_ID(query = name[i],
                        from = "Chemical Name",
                        to = "KEGG")
}

HMDB.ID %>%
  lapply(nrow) %>%
  unlist()

KEGG.ID %>%
  lapply(nrow) %>%
  unlist()

HMDB.ID <-
  HMDB.ID %>%
  do.call(rbind, .) %>%
  as.data.frame()

KEGG.ID <-
  KEGG.ID %>%
  do.call(rbind, .) %>%
  as.data.frame()

id_df <-
  data.frame(HMDB.ID = HMDB.ID$`Human Metabolome Database`,
        KEGG.ID = KEGG.ID$KEGG)

idx <-
  apply(id_df, 1, function(x) {
    sum(is.na(x))
  }) %>%
  `==`(2) %>%
  which()

name[idx]

id_df

idx1 <- which(!is.na(id_df$HMDB.ID) &
        is.na(id_df$KEGG.ID))

idx2 <- which(is.na(id_df$HMDB.ID) &
        !is.na(id_df$KEGG.ID))


metabolomics_data@variable_info <-
  cbind(metabolomics_data@variable_info, id_df)

save(metabolomics_data, file = "data_analysis/data_preparation/metabolomics/metabolomics_data")
