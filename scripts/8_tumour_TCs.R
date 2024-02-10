# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
library(ComplexHeatmap) # ComplexHeatmap_2.16.0
library(circlize) # circlize_0.4.15
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "results/tumour/"

# --- Functions ---
# Computes correlations between mean BCS per signature, TC and patient
bcCorr <- function(bc, tcs = "tumour.tcs", n.spots = 50, method = "pearson") {
  # Metadata
  metadata <- bc@meta.data %>%
    rownames_to_column("spots") %>%
    select(spots, patient, all_of(tcs), subtype) %>%
    rename(group := !!tcs) %>%
    mutate(id = paste(patient, group, sep = "_")) %>%
    add_count(patient, group) %>%
    filter(n >= n.spots) %>%
    select(id, spots, patient, group, subtype)
  # Long normalised BCS for correlation
  normalised <- data.frame(t(bc@normalized))
  normalised.long <- normalised %>%
    rownames_to_column("spots") %>%
    mutate(spots = str_replace(spots, pattern = "\\.", replacement = "-")) %>%
    pivot_longer(colnames(normalised), names_to = "signature",
                 values_to = "enrichment") %>%
    mutate(enrichment = if_else(is.na(enrichment), 0, enrichment)) %>%
    right_join(metadata, by = c("spots"))
  # Edit metadata
  metadata <- metadata %>%
    select(-spots) %>%
    unique() %>%
    as.data.frame() %>%
    rename(!!tcs := group)
  rownames(metadata) <- metadata$id
  # Mean BCS per patient, major TC and signature
  mean.bcs <- normalised.long %>%
    group_by(patient, group, signature) %>%
    summarise(mean.bcs = mean(enrichment)) %>%
    mutate(id = paste(patient, group, sep = "_")) %>%
    ungroup()
  # Create mean BCS matrix
  matrix.bcs <- mean.bcs %>%
    select(-patient, -group) %>%
    pivot_wider(names_from = signature, values_from = mean.bcs) %>%
    column_to_rownames("id")
  matrix.bcs[is.na(matrix.bcs)] <- 0
  # Create correlation matrix
  matrix.corr <- cor(t(matrix.bcs), method = method)
  # Output
  return(list(corr = matrix.corr, metadata = metadata))
}

# --- Data ---
# Merged BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "sensitivity/"), recursive = TRUE,
           showWarnings = FALSE)
dir.create("figures/SuppFigure3/", recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Format metadata
bcmerged@meta.data <- bcmerged@meta.data %>%
  mutate(subtype = if_else(subtype == "ER", "Luminal", subtype),
         subtype = factor(subtype, levels = c("Luminal", "TNBC", "HER2+")))

# Subset tumour spots from merged BCS
tumour.rich <- bcmerged@meta.data %>%
  filter(Tumour == "Tumour") %>%
  rownames()

bctumour <- bcSubset(bcmerged, cells = tumour.rich)

# Restore original Beyondcell object and regress the patient effect
bctumour <- bcRecompute(bctumour, slot = "data")
bctumour@regression <- list(order = c("", ""), vars = NULL)
bctumour <- bcRegressOut(bctumour, vars.to.regress = "patient")

# Recluster
bctumour <- bcUMAP(bctumour, pc = 10, k.neighbors = 20, npcs = 50, res = 0.15)
bctumour@meta.data <- bctumour@meta.data %>%
  mutate(bc_clusters_res.0.15 = as.numeric(as.character(bc_clusters_res.0.15)),
         initial.tumour.tcs = factor(paste0("TC", bc_clusters_res.0.15 + 1)))

initial.tcs <- bcClusters(bctumour, idents = "initial.tumour.tcs") +
  scale_color_manual(values = colour.initial.tcs) +
  theme_classic()
initial.tcs

umap.initial.tumour <- initial.tcs |
  (bcClusters(bctumour, idents = "patient") +
     scale_colour_manual(values = scales::hue_pal()(length(patient.order)),
                         breaks = patient.order) +
     theme_classic())
umap.initial.tumour

ggsave(plot = umap.initial.tumour,  width = 18,
       filename = paste0("figures/UMAP/UMAPinitial_tumour.png"))

# Pearson correlation between tumour TCs
# (for samples with > 50 spots per TCs and patient)
corr.results <- bcCorr(bctumour, tcs = "initial.tumour.tcs", n.spots = 50,
                       method = "pearson")
meta.heat <- corr.results$metadata
corr <- corr.results$corr

# Heatmap
col.fun <- colorRamp2(c(-1, 0, 1), c("#9832C2", "white", "orange"))
col.anno <- HeatmapAnnotation("Initial tumour TC" =
                                meta.heat[colnames(corr), "initial.tumour.tcs"],
                              "Subtype" = meta.heat[colnames(corr), "subtype"],
                              col = list("Initial tumour TC" = colour.initial.tcs,
                                         "Subtype" = colour.subtype))

heat <- Heatmap(corr, name = "Pearson\nCorrelation", col = col.fun,
                show_column_dend = TRUE, show_row_dend = TRUE,
                top_annotation = col.anno,
                column_split = 4, row_split = 4,
                column_title = NULL, row_title = NULL,
                column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
                clustering_method_columns = "ward.D2",
                clustering_method_rows = "ward.D2",
                column_labels = str_remove(colnames(corr), pattern = "_.*$"),
                row_labels = str_remove(rownames(corr), pattern = "_.*$"))

heat

# Relabel tumour TCs according to heatmap results
patients.subcluster2 <- c("1168993F", "1160920F", "1142243F", "CID4290")

bctumour@meta.data <- bctumour@meta.data %>%
  mutate(tumour.tcs = case_when(
    initial.tumour.tcs %in% c("TC4", "TC5") ~ "TC3",
    initial.tumour.tcs == "TC1" & patient %in% patients.subcluster2 ~ "TC1.2",
    initial.tumour.tcs == "TC1" ~ "TC1.1",
    TRUE ~ initial.tumour.tcs),
    tumour.tcs = factor(tumour.tcs))

# Replot UMAP
final.tcs <- bcClusters(bctumour, idents = "tumour.tcs") +
  scale_color_manual(values = colour.tcs) +
  theme_classic()
final.tcs

umap.final.tumour <- initial.tcs | final.tcs |
  (bcClusters(bctumour, idents = "patient") +
     scale_colour_manual(values = scales::hue_pal()(length(patient.order)),
                         breaks = patient.order) +
     theme_classic())
umap.final.tumour

umap.final.patients <- bcClusters(bctumour, idents = "patient",
                                  split.by = "patient") +
  facet_grid(~factor(patient, levels = patient.order)) +
  scale_colour_manual(values = scales::hue_pal()(length(patient.order)),
                      breaks = patient.order) +
  theme_classic()
umap.final.patients

ggsave(plot = umap.final.tumour,  width = 28,
       filename = paste0("figures/UMAP/UMAPfinal_tumour.png"))
ggsave(plot = umap.final.patients,  width = 28,
       filename = paste0("figures/UMAP/UMAPfinal_tumour_patients.png"))

# Redo Pearson correlation
corr.results <- bcCorr(bctumour, tcs = "tumour.tcs", n.spots = 50,
                       method = "pearson")
meta.heat <- bctumour@meta.data %>%
  select(patient, initial.tumour.tcs, tumour.tcs) %>%
  right_join(meta.heat, by = c("patient", "initial.tumour.tcs")) %>%
  select(-id, -subtype) %>%
  unique() %>%
  right_join(corr.results$metadata, by = c("patient", "tumour.tcs")) %>%
  column_to_rownames("id")
corr <- corr.results$corr

col.anno <- HeatmapAnnotation("Initial tumour TC" =
                                meta.heat[colnames(corr), "initial.tumour.tcs"],
                              "Final tumour TC" =
                                meta.heat[colnames(corr), "tumour.tcs"],
                              "Subtype" = meta.heat[colnames(corr), "subtype"],
                              col = list("Initial tumour TC" = colour.initial.tcs,
                                         "Final tumour TC" = colour.tcs,
                                         "Subtype" = colour.subtype))

heat <- Heatmap(corr, name = "Pearson\nCorrelation", col = col.fun,
                show_column_dend = TRUE, show_row_dend = TRUE,
                top_annotation = col.anno,
                column_split = 4, row_split = 4,
                column_title = NULL, row_title = NULL,
                column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
                clustering_method_columns = "ward.D2",
                clustering_method_rows = "ward.D2",
                column_labels = str_remove(colnames(corr), pattern = "_.*$"),
                row_labels = str_remove(rownames(corr), pattern = "_.*$"))

pdf(paste0("figures/Figure4/correlation_drugs.pdf"), width = 20, height = 15)
heat
dev.off()

# Save
saveRDS(bctumour, file = paste0(out.dir, "sensitivity/tcs_tumour.rds"))
