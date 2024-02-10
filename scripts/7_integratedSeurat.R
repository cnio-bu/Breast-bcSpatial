# conda activate pub-beyondcell
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(limma) # limma_3.56.2
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "results/all/"

# --- Data ---
# Merged sensitivity BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Paths to all Seurat objects
seuratobjs.path <-
  list.files("data/visium", pattern = ".rds", full.names = TRUE)
patients <- str_remove(basename(seuratobjs.path), pattern = "_.*$")
seuratobjs.path <- setNames(seuratobjs.path, patients)

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Create a list of Seurat objects
seurats <- lapply(patients, FUN = function(p) {
  seuratobj <- readRDS(seuratobjs.path[p])
  return(seuratobj)
})

# Select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(seurats, nfeatures = 3000)

# Identify anchors
seurats <- PrepSCTIntegration(seurats, anchor.features = features)
anchors <- FindIntegrationAnchors(seurats, normalization.method = "SCT",
                                  anchor.features = features)
saveRDS(anchors, file = "results/all/integration_anchors.rds")

# Create integrated Seurat object
seurat.integrated <- IntegrateData(anchorset = anchors,
                                   normalization.method = "SCT")
saveRDS(seurat.integrated, file = "results/all/integrated_seurat.rds")

# Add metadata
new.metadata <- bcmerged@meta.data %>%
  select(tcs, Tumour, major.tcs, main.cell.type)
seurat.integrated <- AddMetaData(seurat.integrated, new.metadata)

# UMAP reduction
seurat.integrated <- RunPCA(seurat.integrated, npcs = 50)
ElbowPlot(seurat.integrated, ndims = 50, reduction = "pca")
seurat.integrated <- RunUMAP(seurat.integrated, reduction = "pca", dims = 1:30,
                             n.components = 2)

# Plot metadata
DimPlot(seurat.integrated, group.by = "patient") |
  DimPlot(seurat.integrated, group.by = "subtype")
DimPlot(seurat.integrated, group.by = "Tumour", split.by = "major.tcs")
DimPlot(seurat.integrated, group.by = "main.cell.type")

# Log normalise counts
DefaultAssay(seurat.integrated) <- "Spatial"
seurat.integrated <- NormalizeData(seurat.integrated)
seurat.integrated <-
  FindVariableFeatures(seurat.integrated, selection.method = "vst",
                       nfeatures = 2000)
seurat.integrated <- ScaleData(seurat.integrated)

# Differential expression between major TCs
Idents(seurat.integrated) <- "major.tcs"
major.tcs <- unique(seurat.integrated@meta.data$major.tcs)

dea <- lapply(major.tcs, FUN = function(x) {
  markers <- FindMarkers(seurat.integrated, ident.1 = x, min.pct = 0,
                         logfc.threshold = 0, test.use = "wilcox")
  markers <- markers %>%
    rownames_to_column("gene") %>%
    mutate(major.tcs = x)
  return(markers)
}) %>%
  bind_rows()

write.table(dea, file = paste0(out.dir, "DEA_majorTCs.tsv"), sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Save
saveRDS(seurat.integrated, file = "results/all/integrated_seurat.rds")
