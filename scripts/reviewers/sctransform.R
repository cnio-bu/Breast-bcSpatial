# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/reviewers/"

# --- Functions --
# Log normalises Seurat object
logNorm <- function(seuratobj) {
  seuratobj <- NormalizeData(seuratobj, assay = "Spatial")
  seuratobj <- ScaleData(seuratobj, features = rownames(seuratobj), assay = "Spatial")
  return(seuratobj)
}

# Unifies pathologist annotations from different studies
unifyMetadata <- function(seuratobj) {
  lvls <- c("Cancer", "Stroma", "Immune cells", "Fibrous tissue",
            "Adipose tissue", "Normal tissue", "Necrosis")
  invasive.others <-
    c("Invasive cancer + lymphocytes",
      "Invasive cancer + adipose tissue + lymphocytes",
      "Invasive cancer + stroma + lymphocytes", "Invasive cancer + stroma")
  normal.tissue <- c("Normal duct", "Normal glands + lymphocytes")
  seuratobj@meta.data <- seuratobj@meta.data %>%
    mutate(Pathologist =
             case_when(is.na(Pathologist) ~ "",
                       Pathologist == "Invasive_carcinoma" ~ "Cancer",
                       Pathologist == "Cancer trapped in lymphocyte aggregation" ~
                         "Cancer",
                       Pathologist %in% invasive.others ~
                         "Cancer",
                       Pathologist %in% c("Lymphocytes", "TLS") ~
                         "Immune cells",
                       Pathologist == "Fat" ~ "Adipose tissue",
                       Pathologist %in% normal.tissue ~ "Normal tissue",
                       TRUE ~ Pathologist),
           Pathologist = str_replace(Pathologist, pattern = "_",
                                     replacement = " "),
           Pathologist = factor(Pathologist, levels = lvls))
  return(seuratobj)
}

# Computes the correlation between normalised and raw counts
countCorrelation <- function(seuratobj) {
  seuratobj <-
    GroupCorrelation(seuratobj, group.assay = "Spatial", assay = "Spatial",
                     slot = "data", do.plot = FALSE)
  seuratobj <-
    GroupCorrelation(seuratobj, group.assay = "Spatial", assay = "SCT",
                     slot = "scale.data", do.plot = FALSE)
  return(seuratobj)
}

# Plots the correlation between normalised and raw counts,
# grouped by expression group
plotCorrelation <- function(seuratobj, title = FALSE) {
  if (title) {
    title.log <- "Log normalisation"
    title.sct <- "SCTransform normalisation"
  } else title.log <- title.sct <- NULL
  cor.log <- seuratobj@assays$Spatial@meta.features %>%
    pull(nCount_Spatial_cor) %>%
    na.omit()
  cor.sct <- seuratobj@assays$SCT@meta.features %>%
    pull(nCount_Spatial_cor) %>%
    na.omit()
  limits <- pretty(c(cor.log, cor.sct))
  limits <- limits[c(1, length(limits))]
  colour.genes <- scales::brewer_pal(palette = 'YlOrRd')(n = 7)
  box.log <- 
    GroupCorrelationPlot(seuratobj, assay = "Spatial", cor = "nCount_Spatial_cor") +
    scale_fill_manual(labels = 1:6, values = rev(colour.genes)) +
    ggtitle(title.log) +
    ylim(limits) +
    labs(y = "Correlation with raw counts", fill = "Gene bin") +
    theme(text = element_text(size = 18), legend.position = "right",
          plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
  box.sct <- 
    GroupCorrelationPlot(seuratobj, assay = "SCT", cor = "nCount_Spatial_cor") +
    scale_fill_manual(labels = 1:6, values = rev(colour.genes)) +
    ggtitle(title.sct) +
    ylim(limits) +
    labs(y = NULL, fill = "Gene bin") +
    theme(text = element_text(size = 18), legend.position = "right",
          plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
  return((box.log | box.sct) + plot_layout(guides = "collect"))
}

# --- Data ---
# Seurat objects
seuratobj.lum <- readRDS("data/visium/CID4290_spatial_SCTnormalised_clones.rds")
seuratobj.tnbc1 <- readRDS("data/visium/CID44971_spatial_SCTnormalised_clones.rds")
seuratobj.tnbc2 <- readRDS("data/visium/1160920F_spatial_SCTnormalised_clones.rds")

# --- Code ---
# Create outdirs
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Aesthetics
colour.pathologist <-
  c("Cancer" = "#F0027F", "Stroma" = "#7FC97F", "Immune cells" = "#386CB0", 
    "Fibrous tissue" = "#FDC086", "Adipose tissue" = "#FFFF99",
    "Normal tissue" = "#BF5B17", "Necrosis" = "black")

# Log-normalisation
seuratobj.lum <- logNorm(seuratobj.lum)
seuratobj.tnbc1 <- logNorm(seuratobj.tnbc1)
seuratobj.tnbc2 <- logNorm(seuratobj.tnbc2)

# Unify pathologist annotations
seuratobj.lum <- unifyMetadata(seuratobj.lum)
seuratobj.tnbc1 <- unifyMetadata(seuratobj.tnbc1)
seuratobj.tnbc2 <- unifyMetadata(seuratobj.tnbc2)

# Keep non-NA annotations
seuratobj.lum.subset <- subset(seuratobj.lum, Pathologist != "")
seuratobj.tnbc1.subset <- subset(seuratobj.tnbc1, Pathologist != "")
seuratobj.tnbc2.subset <- subset(seuratobj.tnbc2, Pathologist != "")

# Pathologist annotations
patho.lum <- SpatialDimPlot(seuratobj.lum.subset, group.by = "Pathologist") +
  scale_fill_manual(values = scales::hue_pal()(7), drop = FALSE)
patho.tnbc1 <- SpatialDimPlot(seuratobj.tnbc1.subset, group.by = "Pathologist") +
  scale_fill_manual(values = scales::hue_pal()(7), drop = FALSE)
patho.tnbc2 <- SpatialDimPlot(seuratobj.tnbc2.subset, group.by = "Pathologist") +
  scale_fill_manual(values = scales::hue_pal()(7), drop = FALSE)

pathologist <- (patho.lum /patho.tnbc1 / patho.tnbc2) + plot_layout(guides = "collect")
pathologist

# Spatial projection of raw counts
spatialcounts.lum <- 
  SpatialFeaturePlot(seuratobj.lum, features = "nCount_Spatial", 
                     image.alpha = 0.7, pt.size = 1.7) +
  ggtitle("Raw counts") +
  theme(text = element_text(size = 18), legend.position = "bottom", legend.title = element_blank(),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
spatialcounts.tnbc1 <- 
  SpatialFeaturePlot(seuratobj.tnbc1, features = "nCount_Spatial", 
                     image.alpha = 0.7, pt.size = 1.7) +
  theme(text = element_text(size = 18), legend.position = "bottom", legend.title = element_blank())
spatialcounts.tnbc2 <- 
  SpatialFeaturePlot(seuratobj.tnbc2, features = "nCount_Spatial", 
                     image.alpha = 0.7, pt.size = 2.5) +
  theme(text = element_text(size = 18), legend.position = "bottom", legend.title = element_blank())

spatialcounts <- spatialcounts.lum / spatialcounts.tnbc1 /spatialcounts.tnbc2
spatialcounts

# Violin plots of raw counts
vln.lum <- 
  VlnPlot(seuratobj.lum.subset, features = "nCount_Spatial",
          group.by = "Pathologist", same.y.lims = TRUE) +
  scale_fill_manual(values = colour.pathologist, drop = FALSE) +
  ggtitle("Raw counts") +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 18), legend.position = "none",
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
vln.tnbc1 <- 
  VlnPlot(seuratobj.tnbc1.subset, features = "nCount_Spatial",
          group.by = "Pathologist", same.y.lims = TRUE) +
  scale_fill_manual(values = colour.pathologist, drop = FALSE) +
  ggtitle(NULL) +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 18), legend.position = "none")
vln.tnbc2 <- 
  VlnPlot(seuratobj.tnbc2.subset, features = "nCount_Spatial",
          group.by = "Pathologist", same.y.lims = TRUE) +
  scale_fill_manual(values = colour.pathologist, drop = FALSE) +
  ggtitle(NULL) +
  labs(x = NULL, y = NULL) +
  theme(text = element_text(size = 18), legend.position = "none")

vln <- vln.lum / vln.tnbc1 / vln.tnbc2
vln

# Compute correlation between raw and normalised counts
seuratobj.lum <- countCorrelation(seuratobj.lum)
seuratobj.tnbc1 <- countCorrelation(seuratobj.tnbc1)
seuratobj.tnbc2 <- countCorrelation(seuratobj.tnbc2)

# Correlation plots
box.lum <- plotCorrelation(seuratobj.lum, title = TRUE) & theme(legend.position = "none")
box.tnbc1 <- plotCorrelation(seuratobj.tnbc1, title = FALSE) & theme(legend.position = "none")
box.tnbc2 <- plotCorrelation(seuratobj.tnbc2, title = FALSE) & theme(legend.position = "bottom")

box <- box.lum / box.tnbc1 / box.tnbc2
box

# Final plot
fig.sctransform <- (spatialcounts | vln | box) +
  plot_layout(widths = c(1.5, 1, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))
fig.sctransform

ggsave(plot = fig.sctransform, width = 21, height = 16,
       filename = paste0(out.dir, "SCTransform.pdf"))
