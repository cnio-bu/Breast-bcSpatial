# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/reviewers/"

# --- Data ---
# Seurat object
seuratobj.her2 <-
  readRDS("results/tumour/ROI/V19L29_spatial_SCTnormalised_clones_ROI.rds")

# --- Code ---
# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Compute trastuzumab response
trastuzumab <- GenerateGenesets("data/signatures/immunotherapy.gmt")
bcobj.her2.trastuzumab <- bcScore(seuratobj.her2, trastuzumab, expr.thres = 0.1)

response.trastuzumab <-
  bcSignatures(bcobj.her2.trastuzumab, signatures = list(values = "LIU_TRASTUZUMAB_RESPONSE"),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) |>
  lapply(FUN = function(x) {
    x + ggtitle(NULL)
  }) |>
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Trastuzumab sensitivity",
                  subtitle = "HER2 inhibitor, Monoclonal antibody") &
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5),
        legend.position = "bottom")

response.trast.lapatinib <-
  bcSignatures(bcobj.her2.trastuzumab, signatures = list(values = "DIAZ-GIL_TRASTUZUMAB:LAPATINIB_RESPONSE"),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) |>
  lapply(FUN = function(x) {
    x + ggtitle(NULL)
  }) |>
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Trastuzumab + Lapatinib sensitivity",
                  subtitle = "HER2 inhibitor, Monoclonal antibody and Tyrosine kinase inhibitor") &
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5),
        legend.position = "bottom")
response.trast.lapatinib

fig.response.her2 <- response.trastuzumab / response.trast.lapatinib +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))
fig.response.her2

ggsave(plot = fig.response.her2, width = 21, height = 16,
       filename = paste0(out.dir, "HER2_response.pdf"))
