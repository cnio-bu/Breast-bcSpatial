# conda activate pub-semla
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(semla) # semla_1.1.6
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "results/tumour/ROI/"

# --- Data ---
# Seurat object
seuratobj <- readRDS("data/visium/V19L29_spatial_SCTnormalised_clones.rds")

# Tumour sensitivity BCS with TCs
bctumour <- readRDS("results/tumour/sensitivity/tcs_tumour.rds")

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Dimensionality reduction: run PCA
seuratobj <- RunPCA(seuratobj, assay = "SCT", npcs = 50,
                    features = VariableFeatures(seuratobj))

# Elbow plot
elbow <- ElbowPlot(seuratobj, ndims = 50, reduction = "pca")
elbow

# Clustering
seuratobj <- FindNeighbors(seuratobj, reduction = "pca", dims = 1:30)
seuratobj <- FindClusters(seuratobj, resolution = 0.25, verbose = FALSE)

# Run UMAP
seuratobj <- RunUMAP(seuratobj, reduction = "pca", dims = 1:30,
                     n.components = 2)

# Plot clusters
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(SCT_snn_res.0.25 = factor(as.numeric(as.character(SCT_snn_res.0.25)) + 1))
Idents(seuratobj) <- "SCT_snn_res.0.25"
DimPlot(seuratobj, group.by = c("ident", "Phase"))
DimPlot(seuratobj, group.by = c("ident", "slide"))
SpatialDimPlot(seuratobj, group.by = "ident")

# Subcluster cluster 7
seuratobj <- FindSubCluster(seuratobj, "7", resolution = 0.075,
                            graph.name = "SCT_snn")
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(sub.cluster = as.character(sub.cluster),
         suffix = str_extract(sub.cluster, pattern = "_[0-9]{1}$"),
         suffix = as.numeric(str_remove(suffix, pattern = "_")) + 1,
         sub.cluster = if_else(is.na(suffix), sub.cluster,
                               paste(SCT_snn_res.0.25, suffix, sep = "_")))
SpatialDimPlot(seuratobj, group.by = "sub.cluster")

# Define ROIs from expression clusters
ROI <- seuratobj@meta.data %>%
  mutate(ROI = case_when(sub.cluster == "7_2" ~ "ROI1",
                         sub.cluster == "4" ~ "ROI2",
                         sub.cluster == "8" ~ "ROI3",
                         sub.cluster == "6" ~ "ROI4",
                         sub.cluster == "7_1" ~ "ROI5",
                         sub.cluster == "3" ~ "ROI6",
                         sub.cluster == "9" ~ "ROI7",
                         sub.cluster == "5" ~ "ROI8")) %>%
  select(ROI)
seuratobj <- AddMetaData(seuratobj, ROI)
SpatialDimPlot(seuratobj, group.by = "ROI")

# Disconnect regions
semlaobj <- UpdateSeuratForSemla(seuratobj, image_type = "tissue_lowres",
                                 verbose = FALSE)
semlaobj <- DisconnectRegions(semlaobj, column_name = "ROI",
                              selected_groups = paste0("ROI", 1:8))

# Add metadata to Seurat object
new.metadata <- semlaobj@meta.data %>%
  select(ends_with("_split"))
seuratobj <- AddMetaData(seuratobj, new.metadata)
for (i in 1:8) {
  SpatialDimPlot(seuratobj, group.by = paste0("ROI", i, "_split")) |>
    print()
}

# Keep main regions
ROI <- seuratobj@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, ends_with("_split"))

for (i in 1:8) {
  main.regions <- c("S1_region1", "S2_region1")
  ROI <- ROI %>%
    rename(ROI_split := !!paste0("ROI", i, "_split")) %>%
    mutate(ROI_split = if_else(ROI_split %in% main.regions, paste0("ROI", i),
                               NA_character_)) %>%
    rename(!!paste0("ROI", i, "_split") := ROI_split)
}

ROI <- ROI %>%
  pivot_longer(ends_with("_split"), names_to = "split", values_to = "ROI") %>%
  select(-split) %>%
  na.omit() %>%
  column_to_rownames("spots")

# Update metadata
seuratobj <- AddMetaData(seuratobj, ROI)
seuratobj@meta.data <- seuratobj@meta.data %>%
  select(-ends_with("split"))
SpatialDimPlot(seuratobj, group.by = "ROI")

# Keep tumour TCs and subclone annotations that overlap with ROI annotations and
# ROIs with at least 100 spots from the same subclone
tcs.subclones <- bctumour@meta.data %>%
  rownames_to_column("spots") %>%
  filter(patient == "V19L29") %>%
  select(spots, tumour.tcs, subclone)

metadata <- ROI %>%
  rownames_to_column("spots") %>%
  left_join(tcs.subclones, by = "spots") %>%
  mutate(ROI.tcs = if_else(is.na(ROI), NA_character_, tumour.tcs),
         ROI.subclone =
           case_when(subclone == "diploid" ~ NA_character_,
                     is.na(ROI) | is.na(subclone) ~ NA_character_,
                     TRUE ~ subclone)) %>%
  select(spots, ROI, ROI.tcs, ROI.subclone) %>%
  column_to_rownames("spots")

final.ROIs <- metadata %>%
  count(ROI, ROI.subclone) %>%
  na.omit() %>%
  filter(n >= 100) %>%
  pull(ROI) %>%
  unique() %>%
  sort()

metadata <- metadata %>%
  filter(ROI %in% final.ROIs) %>%
  mutate(ROI = factor(paste0("ROI", match(ROI, final.ROIs))))

# Update metadata
seuratobj <- AddMetaData(seuratobj, metadata)
SpatialDimPlot(seuratobj, group.by = "ROI")

# Edge and center of ROIs
semlaobj <- UpdateSeuratForSemla(seuratobj, image_type = "tissue_lowres",
                                 verbose = FALSE)

roi.lvls <- sort(na.omit(unique(semlaobj@meta.data$ROI)))
for (r in roi.lvls) {
  semlaobj <- RegionNeighbors(semlaobj, column_name = "ROI", column_labels = r,
                              mode = "inner")
}

inner.edge <- semlaobj@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, ROI, starts_with("inner_border")) %>%
  pivot_longer(starts_with("inner_border"), names_to = "names",
               values_to = "edge") %>%
  select(-ROI, -names) %>%
  na.omit() %>%
  column_to_rownames("spots")

seuratobj <- AddMetaData(seuratobj, inner.edge)
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(edge = case_when(!is.na(edge) ~ "edge",
                          is.na(edge) & !is.na(ROI) ~ "inner",
                          TRUE ~ NA_character_))

# Subset ROI spots
rois <- levels(seuratobj@meta.data$ROI)
seurat.ROI <- subset(seuratobj, ROI %in% rois)

# Log normalise counts
DefaultAssay(seurat.ROI) <- "Spatial"
seurat.ROI <- NormalizeData(seurat.ROI)
seurat.ROI <- FindVariableFeatures(seurat.ROI, selection.method = "vst",
                                   nfeatures = 2000)
seurat.ROI <- ScaleData(seurat.ROI)

# Differential expression between ROIs
Idents(seurat.ROI) <- "ROI"

dea <- lapply(rois, FUN = function(x) {
  markers <- FindMarkers(seurat.ROI, ident.1 = x, min.pct = 0,
                         logfc.threshold = 0, test.use = "wilcox")
  markers <- markers %>%
    rownames_to_column("gene") %>%
    mutate(ROI = x)
  return(markers)
}) %>%
  bind_rows()

write.table(dea, file = paste0(out.dir, "DEA_ROIs.tsv"), sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Save
saveRDS(seuratobj,
        file = paste0(out.dir, "V19L29_spatial_SCTnormalised_clones_ROI.rds"))
