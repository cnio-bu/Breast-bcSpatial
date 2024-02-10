# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(SPOTlight) # SPOTlight_1.4.1
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/"

# --- Functions ---
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

# --- Data ---
# All metadata
metadata <- read.all.metadata("data/visium/")

# Merged single-cell Seurat object
sc.raw <- readRDS("data/single-cell/single_cell.rds")

# Merged and filtered single-cell Seurat object
sc.postqc <- readRDS("data/single-cell/single_cell_filtered.rds")

# Paths to all Seurat objects
seuratobjs.path <-
  list.files("data/visium", pattern = ".rds", full.names = TRUE)
patients <- str_remove(basename(seuratobjs.path), pattern = "_.*$")
seuratobjs.path <- setNames(seuratobjs.path, patients)

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "SuppFigure1"), recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Number of tumour and TME spots per patient
figSupp1A <- metadata %>%
  mutate(patient = factor(patient, levels = patient.order),
         Tumour = if_else(Tumour == "Non-tumour", "TME", Tumour),
         Tumour = factor(Tumour, levels = c("Tumour", "TME"))) %>%
  ggplot(aes(x = patient, fill = Tumour)) +
  geom_bar(position = position_dodge()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.annotations) +
  labs(x = "Patient ID", y = "Number of spots", fill = "Annotation") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
figSupp1A
ggsave(plot = figSupp1A,
       filename = paste0(out.dir, "SuppFigure1/spot_number.png"))

# Assign each spot to the cell type with the maximum proportion
spot.type <- metadata %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "cell.type",
               values_to = "cell.prop") %>%
  filter(cell.type != "Cancer.Epithelial") %>%
  group_by(spots) %>%
  filter(cell.prop == max(cell.prop)) %>%
  ungroup() %>%
  mutate(main.cell.type = case_when(Tumour == "Tumour" ~ "Cancer.Epithelial",
                                    TRUE ~ cell.type),
         main.cell.type = str_replace(main.cell.type, pattern = "\\.",
                                      replacement = " "))

# Proportion of main cell types per patient
figSupp1B <- spot.type %>%
  mutate(patient = factor(patient, levels = patient.order)) %>%
  count(patient, main.cell.type) %>%
  ggplot(aes(x = patient, y = n, fill = main.cell.type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.celltypes) +
  labs(x = "Patient ID", y = "Proportion", fill = "Cell type") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), axis.line = element_blank())
figSupp1B
ggsave(plot = figSupp1B,
       filename = paste0(out.dir, "SuppFigure1/st_celltypes.png"))

# Number of patients with single-cell data, per subtype
patients.sc <- sc.raw@meta.data %>%
  rename(patient = orig.ident) %>%
  mutate(subtype = if_else(subtype == "ER+", "Luminal", subtype)) %>%
  select(patient, subtype) %>%
  unique()
table(patients.sc$subtype)

# Number of filtered single-cells types per breast cancer subtype, patient and
# cell type
n.cells <- sc.postqc@meta.data %>%
  rename(patient = orig.ident) %>%
  add_count(patient, celltype_major) %>%
  select(patient, subtype, celltype_major, n) %>%
  unique() %>%
  group_by(patient) %>%
  mutate(total = sum(n)) %>%
  ungroup() %>%
  mutate(celltype_major = str_replace(celltype_major, pattern = "-",
                                      replacement = " "),
         subtype = if_else(subtype == "ER+", "Luminal", subtype))

n.cells <- n.cells %>%
  mutate(celltype_major = factor(celltype_major),
         subtype = factor(subtype, levels = c("Luminal", "TNBC", "HER2+")))

figSupp1C <- n.cells %>%
  arrange(subtype, total) %>%
  mutate(patient = factor(patient, levels = unique(patient))) %>%
  ggplot(aes(x = patient, y = n, fill = celltype_major)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.celltypes) +
  labs(x = NULL, y = "Number of cells", fill = "Cell type") +
  facet_grid(cols = vars(subtype), scales = "free_x", space = "free_x") +
  pub_theme_facet +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
figSupp1C
ggsave(plot = figSupp1C,
       filename = paste0(out.dir, "SuppFigure1/sc_celltypes.png"))

# SpatialDimplots with deconvolution metrics, for each patient
spot.type.metadata <- spot.type %>%
  column_to_rownames("spots") %>%
  select(main.cell.type)

deconv.plots <- lapply(seuratobjs.path[patient.order], FUN = function(x) {
  seuratobj <- readRDS(x)
  patient <- unique(seuratobj@meta.data$patient)
  slides <- names(seuratobj@images)
  pt.size <- 2.5
  if (patient == "V19L29") pt.size <- 1.5
  else if (patient %in% c("CID44971", "CID4290")) pt.size <- 1.7
  else if (patient %in% c("738811QB", "CID4535", "CID4465")) pt.size <- 1.85
  lapply(slides, FUN = function(s) {
    if (length(slides) == 1) title <- patient
    else title <- paste0(patient, ", ", s)
    seuratobj <- AddMetaData(seuratobj, spot.type.metadata)
    spatial.plot <-
      SpatialDimPlot(seuratobj, group.by = "main.cell.type", images = s,
                     image.alpha = 0, pt.size.factor = pt.size) +
      scale_fill_manual(values = colour.celltypes,
                        labels = levels(n.cells$celltype_major),
                        guide = guide_legend(reverse = TRUE)) +
      ggtitle(title) +
      theme(text = element_text(size = 18), legend.position = "none",
            plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
    return(spatial.plot)
  })
}) |>
  unlist(recursive = FALSE)

figSupp1D <- wrap_plots(deconv.plots, nrow = 2, widths = 1) +
  plot_layout(guides = "collect") & theme(legend.position = "none")
figSupp1D
ggsave(plot = figSupp1D, width = 21, height = 15,
       filename = paste0(out.dir, "SuppFigure1/deconvolution_plots.png"))

# Save figure
figSupp1 <- ((figSupp1A | figSupp1B) / figSupp1C / figSupp1D) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.25, 2)) &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = figSupp1, width = 21, height = 29.7,
       filename = paste0(out.dir, "SuppFigure1/SuppFigure1.pdf"))
