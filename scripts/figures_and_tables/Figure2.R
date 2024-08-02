# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (1823-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(ggpattern) # ggpattern_1.0.1
library(ComplexHeatmap) # ComplexHeatmap_2.16.0
library(circlize) # circlize_0.4.15
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
set.seed(1)
out.dir <- "figures/"

# --- Functions ---
# Compute correlations between mean BCS per signature, TC and patient
bcCorr <- function(bc, tcs = "major.tcs", n.spots = 50, method = "pearson") {
  # Metadata
  metadata <- bc@meta.data %>%
    rownames_to_column("spots") %>%
    select(spots, patient, all_of(tcs), subtype) %>%
    rename(group := !!tcs) %>%
    filter(group != "Mixed") %>%
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

# Neighbourhood enrichment
neighbours <- read.table("results/all/neighbourhood_enrichment.tsv",
                         sep = "\t", header = TRUE)

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "Figure2"), recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Format metadata
bcmerged@meta.data <- bcmerged@meta.data %>%
  mutate(patient = factor(patient, levels = patient.order),
         Tumour = if_else(Tumour == "Non-tumour", "TME", Tumour),
         Tumour = factor(Tumour, levels = c("Tumour", "TME")),
         subtype = if_else(subtype == "ER", "Luminal", subtype),
         subtype = factor(subtype, levels = c("Luminal", "TNBC", "HER2+")),
         major.tcs = factor(major.tcs, levels = c("Tumour-rich", "Mixed",
                                                  "TME-rich")))

# Pearson correlation between Tumour and TME-rich TCs
# (for samples with > 50 spots per TCs and patient)
corr.results <- bcCorr(bcmerged, tcs = "major.tcs", n.spots = 50,
                       method = "pearson")
meta.heat <- corr.results$metadata %>%
  mutate(major.tcs = factor(major.tcs, levels = c("Tumour-rich", "TME-rich")))
corr <- corr.results$corr

# UMAP coloured by TC, Tumour and TME labels and major TC
fig2A <- (bcClusters(bcmerged, idents = "tcs") + theme_classic() +
            theme(text = element_text(size = 18),
                  axis.text = element_text(size = 16))) |
  (bcClusters(bcmerged, idents = "Tumour") + theme_classic() +
     theme(text = element_text(size = 18),
           axis.text = element_text(size = 16))) |
  (bcClusters(bcmerged, idents = "major.tcs") + theme_classic() +
     scale_color_manual(values = colour.tctypes) +
     theme(text = element_text(size = 18),
           axis.text = element_text(size = 16)))
fig2A
ggsave(plot = fig2A, filename = paste0(out.dir, "Figure2/umap.png"),
       width = 28)

# Proportion of Tumour and TME spots per TC
metadata <- bcmerged@meta.data
fig2B <- metadata %>%
  add_count(tcs, Tumour) %>%
  select(tcs, Tumour, major.tcs, n) %>%
  unique() %>%
  ggplot(aes(x = tcs, y = n, fill = major.tcs, )) +
  geom_bar_pattern(position = "fill", stat = "identity",
                   mapping = aes(pattern = Tumour),
                   pattern_density = 0.15, pattern_spacing = 0.01,
                   pattern_fill = "black", pattern_color = NA)  +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.tctypes) +
  scale_pattern_manual(values = c("none", "crosshatch")) +
  labs(x = "Therapeutic clusters", y = "Proportion", fill = "Major TC",
       alpha = "Annotation") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), axis.line = element_blank())
#fig2B <- metadata %>%
#  add_count(tcs, Tumour) %>%
#  select(tcs, Tumour, major.tcs, n) %>%
#  unique() %>%
#  ggplot(aes(x = tcs, y = n, fill = major.tcs, alpha = Tumour)) +
#  geom_bar(position = "fill", stat = "identity") +
#  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
#  scale_alpha_discrete(range = c(1, 0.6)) +
#  scale_fill_manual(values = colour.tctypes) +
#  labs(x = "Therapeutic clusters", y = "Proportion", fill = "Major TC",
#       alpha = "Annotation") +
#  theme_classic() +
#  theme(text = element_text(size = 18),
#        axis.text = element_text(size = 16),
#        axis.text.x = element_text(angle = 45, hjust = 1),
#        axis.ticks.x = element_blank(), axis.line = element_blank())
fig2B
ggsave(plot = fig2B, filename = paste0(out.dir, "Figure2/proportion.png"))

# Distribution of major TCs accross patients
fig2C <- metadata %>%
  count(patient, major.tcs) %>%
  ggplot(aes(x = patient, y = n, fill = major.tcs)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.tctypes) +
  labs(x = "Patient ID", y = "Proportion", fill = "Major TC") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), axis.line = element_blank())
fig2C
ggsave(plot = fig2C, filename = paste0(out.dir, "Figure2/patient.png"))

# Create major TC proportion table
proportion <- metadata %>%
  add_count(patient) %>%
  add_count(patient, major.tcs, name = "n.major.tcs") %>%
  mutate(prop = round(n.major.tcs/n, digits = 2)) %>%
  select(patient, subtype, major.tcs, prop) %>%
  unique() %>%
  pivot_wider(id_cols = c(patient, subtype), names_from = major.tcs,
              values_from = prop)
proportion <- proportion[match(patient.order, proportion$patient), ]
proportion <- proportion[c("patient", "subtype", "Tumour-rich", "Mixed",
                           "TME-rich")]
write.table(proportion, file = "tables/SuppTable4.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Neighbourhood enrichment
major.lvls <- c("Tumour-rich", "Mixed", "TME-rich")
mean.zscore <- neighbours %>%
  filter(edges >= 50) %>%
  mutate(label_1 = str_remove(label_1, pattern = "Label_"),
         label_2 = str_remove(label_2, pattern = "Label_")) %>%
  add_count(label_1, label_2) %>%
  group_by(label_1, label_2) %>%
  mutate(mean.zscore = mean(z_score),
         label_1 = factor(label_1, levels = major.lvls),
         label_2 = factor(label_2, levels = major.lvls)) %>%
  ungroup()

fig2D <- mean.zscore %>%
  ggplot(aes(x = label_1, y = label_2, fill = mean.zscore)) +
  geom_tile(color = "white", lwd = 1.5, linetype = 1) +
  scale_y_discrete(limits = rev) +
  coord_fixed() +
  scale_fill_gradient2(low = "#0474BA", mid = "grey90", high = "#F17720") +
  labs(x = NULL, y = NULL, fill = "Mean\nzscore") +
  theme_neighbours +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16))
fig2D
ggsave(plot = fig2D, filename = paste0(out.dir, "Figure2/neighbourhood.png"),
       width = 10)

# Three examples
metadata.filtered <- metadata %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs = factor(major.tcs, levels = c("Tumour-rich", "Interphase",
                                                  "TME-rich"))) %>%
  select(spots, tcs, major.tcs) %>%
  unique()

patient <- "1168993F"
subtype <- "HER2+"
bcobj <- readRDS(paste0("data/beyondcell/sensitivity/", patient,
                        "_spatial_bcs_sensitivity.rds"))
bcobj@meta.data <- bcobj@meta.data %>%
  rownames_to_column("spots") %>%
  left_join(metadata.filtered, by = "spots") %>%
  column_to_rownames("spots")
fig.her2 <- bcClusters(bcobj, idents = "major.tcs", spatial = TRUE,
                       pt.size = 2.5, image.alpha = 0.7) +
  scale_fill_manual(values = colour.tctypes) +
  ggtitle(paste0(patient, " (", subtype, ")")) +
  theme(text = element_text(size = 18), legend.position = "none",
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
fig.her2

patient <- "1160920F"
subtype <- "TNBC"
bcobj <- readRDS(paste0("data/beyondcell/sensitivity/", patient,
                        "_spatial_bcs_sensitivity.rds"))
bcobj@meta.data <- bcobj@meta.data %>%
  rownames_to_column("spots") %>%
  left_join(metadata.filtered, by = "spots") %>%
  column_to_rownames("spots")
fig.tnbc <- bcClusters(bcobj, idents = "major.tcs", spatial = TRUE,
                       pt.size = 2.5, image.alpha = 0.7) +
  scale_fill_manual(values = colour.tctypes) +
  ggtitle(paste0(patient, " (", subtype, ")")) +
  theme(text = element_text(size = 18), legend.position = "none",
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
fig.tnbc

patient <- "CID4535"
subtype <- "Luminal"
bcobj <- readRDS(paste0("data/beyondcell/sensitivity/", patient,
                        "_spatial_bcs_sensitivity.rds"))
bcobj@meta.data <- bcobj@meta.data %>%
  rownames_to_column("spots") %>%
  left_join(metadata.filtered, by = "spots") %>%
  column_to_rownames("spots")
fig.er <- bcClusters(bcobj, idents = "major.tcs", spatial = TRUE,
                     pt.size = 1.85, image.alpha = 0.7) +
  scale_fill_manual(values = colour.tctypes) +
  ggtitle(paste0(patient, " (", subtype, ")")) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
fig.er

fig2E <- (fig.er | fig.tnbc | fig.her2) + plot_layout(guides = "collect")
fig2E
ggsave(plot = fig2E, filename = paste0(out.dir, "Figure2/slice_examples.png"))

# Correlation heatmap
col.fun <- colorRamp2(c(-1, 0, 1), c("#9832C2", "white", "orange"))
col.anno <- HeatmapAnnotation("Major TC" = meta.heat[colnames(corr), "major.tcs"],
                              "Subtype" = meta.heat[colnames(corr), "subtype"],
                              col = list("Major TC" = colour.tctypes,
                                         "Subtype" = colour.subtype))
heat <- Heatmap(corr, name = "Pearson\nCorrelation", col = col.fun,
                show_column_dend = TRUE, show_row_dend = TRUE,
                top_annotation = col.anno,
                column_split = 2, row_split = 2,
                column_title = NULL, row_title = NULL,
                column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
                clustering_method_columns = "ward.D2",
                clustering_method_rows = "ward.D2",
                column_labels = str_remove(colnames(corr), pattern = "_.*$"),
                row_labels = str_remove(rownames(corr), pattern = "_.*$"))

pdf(paste0(out.dir, "Figure2/correlation_drugs.pdf"), width = 18, height = 15)
heat
dev.off()

fig2F <- fig(paste0(out.dir, "Figure2/correlation_drugs.pdf"))

# Save figure
fig2 <- fig2A / (fig2B | fig2C | fig2D) / fig2E / fig2F +
  plot_layout(heights = c(1.2, 1.2, 1.5, 2)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = fig2, width = 21, height = 29.7,
       filename = paste0(out.dir, "Figure2/Figure2.pdf"))
