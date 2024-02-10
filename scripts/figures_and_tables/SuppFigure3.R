# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/"

# --- Functions ---
# Add metadata to Seurat object
custom.AddMetaData <- function(seuratobj, bcobj, colnames) {
  metadata <- bcobj@meta.data[, colnames, drop = FALSE] %>%
      rownames_to_column("spots")
  seuratobj@meta.data <- seuratobj@meta.data %>%
    rownames_to_column("spots") %>%
    left_join(metadata, by = "spots") %>%
    column_to_rownames("spots")
  return(seuratobj)
}

# Return Seurat's SpatialFeaturePlot default colourscale
colourscale <- function(n = 100, rev = TRUE) {
  colours <- RColorBrewer::brewer.pal(n = 11, name = "Spectral")
  if (rev) {
    colours <- colours |>
      rev()
  }
  colours <- colours |>
    grDevices::colorRampPalette()
  return(colours(n = n))
}

# --- Data ---
# Merged BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "SuppFigure3"), recursive = TRUE, 
           showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Modify metadata
bcmerged@meta.data <- bcmerged@meta.data %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs = factor(major.tcs, levels = c("Tumour-rich", "Interphase",
                                                  "TME-rich"))) %>%
  column_to_rownames("spots")

# TME spots
spots.tme <- bcmerged@meta.data %>%
  filter(Tumour != "Tumour") %>%
  rownames()

# ER+ example
patient <- "CID4535"
subtype <- "Luminal"
pt.size <- 1.85
bcobj.sensitivity <-
  readRDS(paste0("data/beyondcell/sensitivity/", patient,
                 "_spatial_bcs_sensitivity.rds"))
bcobj.tme <- 
  readRDS(paste0("results/all/tme/", patient, "_spatial_tme_status.rds")) |>
  bcSubset(cells = spots.tme) |>
  bcRecompute(slot = "normalized")
seuratobj <-
  readRDS(paste0("data/visium/", patient, "_spatial_SCTnormalised_clones.rds"))
seuratobj <- custom.AddMetaData(seuratobj, bcmerged, colnames = "major.tcs")

er.compartments <-
  SpatialDimPlot(seuratobj, group.by = "major.tcs", pt.size.factor = pt.size,
                 image.alpha = 0.7) +
  scale_fill_manual(values = colour.tctypes) +
  ggtitle(paste0(patient, " (", subtype, ")")) +
  theme(text = element_text(size = 18), legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold"))
er.compartments

er.biomarker <-
  SpatialFeaturePlot(seuratobj, features = "ESR1", pt.size.factor = pt.size,
                     image.alpha = 0.7) +
  ggtitle("ESR1") +
  theme(text = element_text(size = 18), legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
er.biomarker

er.tumour.response <-
  bcSignatures(bcobj.sensitivity,
               signatures = list(values = "tamoxifen_CTRP_K93754473"),
               spatial = TRUE, pt.size = pt.size, image.alpha = 0.7) +
  ggtitle("Tamoxifen (CTRP) sensitivity", subtitle = "Hormonal therapy") +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
er.tumour.response

lims <- range(bcobj.tme@meta.data$TME.score, bounds = FALSE)
er.tme <- 
  bcClusters(bcobj.tme, idents = "TME.score", spatial = TRUE,
             factor.col = FALSE, pt.size = pt.size) +
  scale_fill_gradientn(colours = c("#7BC57D", "white", "#EB5528"), 
                       values = scales::rescale(c(lims[1], 0, lims[2])),
                       guide = "colorbar", limits = lims) +
  ggtitle("Antitumor TME") +
  labs(fill = "TME status") +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold"))
er.tme

er <- er.compartments / er.biomarker / er.tumour.response / er.tme

# HER2+ example
patient <- "738811QB"
subtype <- "HER2+"
pt.size <- 1.85
bcobj.sensitivity <-
  readRDS(paste0("data/beyondcell/sensitivity/", patient,
                 "_spatial_bcs_sensitivity.rds"))
bcobj.tme <- 
  readRDS(paste0("results/all/tme/", patient, "_spatial_tme_status.rds")) |>
  bcSubset(cells = spots.tme) |>
  bcRecompute(slot = "normalized")
seuratobj <-
  readRDS(paste0("data/visium/", patient, "_spatial_SCTnormalised_clones.rds"))
seuratobj <- custom.AddMetaData(seuratobj, bcmerged, colnames = "major.tcs")

her2.compartments <-
  SpatialDimPlot(seuratobj, group.by = "major.tcs", pt.size.factor = pt.size,
                 image.alpha = 0.7) +
  scale_fill_manual(values = colour.tctypes) +
  ggtitle(paste0(patient, " (", subtype, ")")) +
  labs(fill = "Major TC") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(text = element_text(size = 18), legend.key = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5, face = "bold"))
her2.compartments

her2.biomarker <-
  SpatialFeaturePlot(seuratobj, features = "ERBB2", pt.size.factor = pt.size,
                     image.alpha = 0.7) +
  ggtitle("ERBB2") +
  theme(text = element_text(size = 18), legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))
her2.biomarker

her2.tumour.response <-
  bcSignatures(bcobj.sensitivity,
               signatures = list(values = "afatinib_PRISM_K66175015"),
               spatial = TRUE, pt.size = pt.size, image.alpha = 0.7) +
  ggtitle("Afatinib (PRISM) sensitivity", subtitle = "HER2 inhibitor") +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
her2.tumour.response

lims <- range(bcobj.tme@meta.data$TME.score, bounds = FALSE)
her2.tme <-
  bcClusters(bcobj.tme, idents = "TME.score", spatial = TRUE,
             factor.col = FALSE, pt.size = pt.size) +
  scale_fill_gradientn(colours = c("#7BC57D", "white", "#EB5528"), 
                       values = scales::rescale(c(lims[1], 0, lims[2])),
                       guide = "colorbar", limits = lims) +
  ggtitle("Dual TME") +
  labs(fill = "TME status") +
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, face = "bold"))
her2.tme

her2 <- her2.compartments / her2.biomarker / her2.tumour.response / her2.tme

# Save figure
figSupp3 <- (er | her2) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = figSupp3, width = 21, height = 29.7,
       filename = paste0(out.dir, "SuppFigure3/SuppFigure3.pdf"))
