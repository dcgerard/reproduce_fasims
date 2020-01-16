######################
## Extract pbmc data
## We use same filtering criteria as here: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
######################

library(dplyr)
library(Seurat)

## Extract data ----
pbmc.data <- Read10X(data.dir = "./data/pbmc/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

## Filter based on given criteria ----
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cdat <- as.matrix(GetAssay(pbmc)@counts)

## Find variable genes ----
pbmc_scale <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_scale <- NormalizeData(pbmc_scale)
pbmc_scale <- FindVariableFeatures(pbmc_scale, selection.method = "vst", nfeatures = 4000)
var_features <- VariableFeatures(pbmc_scale)

## Keep only variable genes ----
cdat <- cdat[var_features, ]

## Save data ----
saveRDS(object = cdat, "./output/sc/pbmc_cleaned.RDS")
