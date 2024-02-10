# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(eulerr) # eulerr_7.0.0
library(RColorBrewer) # RColorBrewer_1.1-3
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
set.seed(1)
out.dir <- "figures/"

# --- Functions --
# Read all metadata files in a folder of Seurat objects
read.all.metadata <- function(pathdir) {
  seurat.files <- list.files(pathdir, pattern = ".rds", full.names = TRUE)
  lapply(seurat.files, FUN = function(x) {
    seuratobj <- readRDS(x)
    seuratobj@meta.data %>%
      rownames_to_column("spots")
  }) %>%
    bind_rows()
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
# All metadata
metadata <- read.all.metadata("data/visium/")

# Example Seurat object
patient <- "CID44971"
subtype <- "(TNBC)"
seuratobj <-
  readRDS(paste0("data/visium/", patient, "_spatial_SCTnormalised_clones.rds"))

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "Figure1"), recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Workflow
fig1A <- fig("data/images/workflow.pdf")
ggsave(plot = fig1A, filename = paste0(out.dir, "Figure1/workflow.png"))

# Create new metadata column indicating whether the spot is in a tumour region
# according to the pathologist
metadata <- metadata %>%
  mutate(tumour.patho = Pathologist == "DCIS" |
           str_detect(tolower(Classification), pattern = "cancer") |
           str_detect(tolower(Classification), pattern = "carcinoma"))

# Overlap between deconvolution, ESTIMATE, pathologist and SCEVAN tumour 
# annotations
annotations <- metadata %>%
  select(spots, patient, Cancer.Epithelial, ESTIMATEScaled, tumour.patho,
         subclone) %>%
  mutate(spots.patient = paste(spots, patient, sep = "_"),
         tumour.estimate = ESTIMATEScaled <= 0.4,
         tumour.deconvolution = Cancer.Epithelial >= 0.6,
         tumour.scevan = !is.na(subclone))

overlap <- list(Pathologist = annotations %>%
                  filter(tumour.patho) %>%
                  pull(spots.patient),
                ESTIMATE = annotations %>%
                  filter(tumour.estimate) %>%
                  pull(spots.patient),
                Deconvolution = annotations %>%
                  filter(tumour.deconvolution) %>%
                  pull(spots.patient),
                SCEVAN = annotations %>%
                  filter(tumour.scevan) %>%
                  pull(spots.patient))

# Venn diagram
fit <- euler(overlap)
fig1B <- plot(fit, labels =list(fontsize = 18), quantities = list(fontsize = 15),
              fills = c("#CEBCDB", "darkgoldenrod1", "#A3BEE3", "#AEF067"))
fig1B
ggsave(plot = fig1B, filename = paste0(out.dir, "Figure1/venn.png"))

# Edit patient metadata
seuratobj@meta.data <- seuratobj@meta.data %>%
  mutate(Pathologist = str_replace_all(Pathologist, pattern = "\\+ ",
                                       replacement = "\\+\n"),
         Pathologist = if_else(startsWith(Pathologist, prefix = "Normal"),
                               "Normal + stroma +\nlymphocytes", Pathologist),
         subclone = factor(str_to_title(subclone)),
         Tumour = if_else(Tumour == "Non-tumour", "TME", Tumour),
         Tumour = factor(Tumour, levels = c("Tumour", "TME")))

# Histopathological annotations
fig1C <- SpatialDimPlot(seuratobj, group.by = "Pathologist",
                        image.alpha = 0.7, pt.size = 1.7) +
  scale_fill_manual(values = colour.histo) +
  ggtitle(paste(patient, subtype)) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(text = element_text(size = 18), legend.key = element_rect(fill = NA),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
fig1C
ggsave(plot = fig1C, filename = paste0(out.dir, "Figure1/pathologist_anno.png"))

# Tumour purity
fig1D <- SpatialFeaturePlot(seuratobj, features = "ESTIMATEScaled",
                            image.alpha = 0.7, pt.size = 1.7) +
  scale_fill_gradientn(colours = colourscale(n = 100, rev = FALSE),
                       limits = c(0, 1)) +
  labs(fill = "ESTIMATE") +
  theme(text = element_text(size = 18), legend.position = "right")
fig1D
ggsave(plot = fig1D, filename = paste0(out.dir, "Figure1/tumour_purity.png"))

# Deconvolution results
fig1E <- SpatialFeaturePlot(seuratobj, features = "Cancer.Epithelial",
                            image.alpha = 0.7, pt.size = 1.7) +
  scale_fill_gradientn(colours = colourscale(n = 100, rev = TRUE),
                       limits = c(0, 1)) +
  labs(fill = "Cancer Epithelial") +
  theme(text = element_text(size = 18), legend.position = "right")
fig1E
ggsave(plot = fig1E,
       filename = paste0(out.dir, "Figure1/cancer_deconvolution.png"))

# SCEVAN results
fig1F <- SpatialDimPlot(seuratobj, group.by = "subclone", image.alpha = 0.7, 
                        pt.size = 1.7) +
  scale_fill_manual(values = colour.subclones, na.value = NA) +
  labs(fill = "Subclone") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(text = element_text(size = 18), legend.key = element_rect(fill = NA))
fig1F
ggsave(plot = fig1F, filename = paste0(out.dir, "Figure1/subclones.png"))

# Tumour and TME annotations
fig1G <- SpatialDimPlot(seuratobj, group.by = "Tumour", image.alpha = 0.7, 
                        pt.size = 1.7) +
  scale_fill_manual(values = colour.annotations) +
  labs(fill = "Annotation") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(text = element_text(size = 18), legend.key = element_rect(fill = NA))
fig1G
ggsave(plot = fig1G, filename = paste0(out.dir, "Figure1/tumour_tme.png"))

# Save figure
bottom <- (ggplotify::as.ggplot(fig1B) | fig1C) /
  (fig1D | fig1E) /
  (fig1F | fig1G)

fig1 <- fig1A / bottom +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = fig1, width = 21, height = 28,
       filename = paste0(out.dir, "Figure1/Figure1.pdf"))
