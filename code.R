library(dplyr)
library(Seurat)
library(SeuratWrappers)
library(SingleR)


## Read data, define sample by 10x barcode index
data <- Read10X("filtered_feature_bc_matrix")
seurat_obj <- CreateSeuratObject(counts = data)
seurat_obj$batch <- as.character(sapply(strsplit(rownames(seurat_obj@meta.data), split="-"), "[[", 2))

## Quality control, median was set by nFeature median
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
nFeature_Q3 <- as.numeric(quantile(seurat_obj$nFeature_RNA, 0.75))
nFeature_Q1 <- as.numeric(quantile(seurat_obj$nFeature_RNA, 0.25))
IQR <- nFeature_Q3 - nFeature_Q1
upthred <- round(nFeature_Q3 + 1.5 * IQR)
seurat_obj <- seurat_obj[, which(seurat_obj$nFeature_RNA < upthred)]
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 20)

Feature_median <- median(seurat_obj$nFeature_RNA)
if (Feature_median < 2000){
	Feature_median <- 2000
	}else{
		Feature_median <- Feature_median
	}

## Normalized
seurat_obj_mnn <- NormalizeData(seurat_obj) %>%
    FindVariableFeatures(nfeatures = round(Feature_median)) 
vargene <- seurat_obj_mnn[["RNA"]]@var.features
## Batch effect correction
seurat_obj_mnn <- RunFastMNN(object.list = SplitObject(seurat_obj_mnn, split.by = "batch")) 
## Dimensional reduction and clustering
seurat_obj_mnn <- ScaleData(seurat_obj_mnn) %>%
    RunPCA(features=vargene) %>%
    RunTSNE(reduction="mnn", dims = 1:20, perplexity = 30, seed.use = 100) %>%
    RunUMAP(reduction="mnn", dims = 1:20, n.components = 2, seed.use = 100) %>%
    FindNeighbors(dims = 1:20) %>%
    FindClusters(resolution = 0.4)

## SingleR	
singler <- CreateSinglerObject(seurat_obj_mnn[['RNA']]@counts,
	species = "Human", project.name="PBMC", 
	clusters = seurat_obj_mnn@active.ident)
    
