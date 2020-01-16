## Format GTEx data and exctract muscle RNA
library(vroom)
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(tidyverse))

## Get gene counts ------------------------------------------------------------
gene_counts <- vroom(file = "./data/gtex/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_reads.gct", skip = 2, delim = "\t")
row_data <- gene_counts[, c("Name", "Description")]
gene_counts <- data.matrix(gene_counts[, -(1:2)])

## Get sample attributes ------------------------------------------------------
col_data <- vroom(file = "./data/gtex/GTEx_v7_Annotations_SampleAttributesDS.txt", delim = "\t")
matchvec <- match(colnames(gene_counts), col_data$SAMPID)
col_data <- col_data[matchvec, ]
stopifnot(all(col_data$SAMPID == colnames(gene_counts)))

## Get subject data -----------------------------------------------------------
subj_data <- vroom(file = "./data/gtex/GTEx_v7_Annotations_SubjectPhenotypesDS.txt", delim = "\t")
col_data %>%
  mutate(SUBJID = str_extract(SAMPID, "^[^-]+-[^-]+")) %>%
  left_join(subj_data, by = "SUBJID") ->
  col_data

## Validate primary key -------------------------------------------------------
col_data %>%
  count(SMTSD, SUBJID) %>%
  filter(n > 1) %>%
  nrow() -> keynum
stopifnot(keynum == 0)

## Set as SummarizedExperiment object -----------------------------------------
se <- SummarizedExperiment(assays = gene_counts, rowData = row_data, colData = col_data)

## Extract various tissues ----------------------------------------------------
tissue_vec <- sort(unique(colData(se)$SMTSD))
tissue_vec %>%
  str_replace_all("-", "") %>%
  str_replace_all("\\(", "") %>%
  str_replace_all("\\)", "") %>%
  str_squish() %>%
  str_to_lower() %>%
  str_replace_all("\\s+", "_") %>%
  str_c("./output/tissue_data/", ., ".RDS") ->
  name_vec

for (index in seq_along(tissue_vec)) {
  sesmall <- se[, colData(se)$SMTSD == tissue_vec[index]]
  saveRDS(object = sesmall, file = name_vec[index])
}
