# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(ggsankey) # ggsankey_0.0.99999
library(rstatix) # rstatix_0.7.2
library(RColorBrewer) # RColorBrewer_1.1-3
library(viridis) # viridis_0.6.4
library(ggpubr) # ggpubr_0.6.0
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/"

# --- Functions ---
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
# Seurat object with ROI annotation
seuratobj <-
  readRDS("results/tumour/ROI/V19L29_spatial_SCTnormalised_clones_ROI.rds")

# Beyondcell object
bcobj <-
  readRDS(paste0("data/beyondcell/sensitivity/V19L29_spatial_bcs_sensitivity.rds"))

# CNA annotation
cna <- read.table("results/tumour/ROI/CNA_ROIs.tsv", header = TRUE, sep = "\t")

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "Figure5"), recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")
colour.all <- c(setNames(colour.subclones, paste0("SC", names(colour.subclones))),
                colour.tcs, colour.rois)

# ROIs
fig5A <- SpatialDimPlot(seuratobj, group.by = "ROI", pt.size = 1.5,
                        image.alpha = 0.7, combine = FALSE) |>
  lapply(FUN = function(x) {
    title <- str_remove(x$labels$title, pattern = "slide")
    x + scale_fill_manual(values = colour.rois, na.value = NA) +
      ggtitle(paste0("V19L29, slide", title)) +
      guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3))) +
      theme(text = element_text(size = 18), 
            legend.key = element_rect(fill = NA),
            plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
  }) |>
  wrap_plots(nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig5A
ggsave(plot = fig5A, filename = paste0(out.dir, "/Figure5/ROI.png"))

# Subclones
fig5B <- SpatialDimPlot(seuratobj, group.by = "ROI.subclone", pt.size = 1.5,
                        image.alpha = 0.7, combine = FALSE) |>
  lapply(FUN = function(x) {
    x + scale_fill_manual(values = colour.subclones, na.value = NA) +
      ggtitle(NULL) +
      labs(fill = "Subclone") +
      guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3))) +
      theme(text = element_text(size = 18),
            legend.key = element_rect(fill = NA))
  }) |>
  wrap_plots(nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig5B
ggsave(plot = fig5B, filename = paste0(out.dir, "/Figure5/subclones.png"))

# TCs
fig5C <- SpatialDimPlot(seuratobj, group.by = "ROI.tcs", pt.size = 1.5,
                        image.alpha = 0.7, combine = FALSE) |>
  lapply(FUN = function(x) {
    x + scale_fill_manual(values = colour.tcs, na.value = NA) +
      ggtitle(NULL) +
      labs(fill = "Tumour TC") +
      guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3))) +
      theme(text = element_text(size = 18),
            legend.key = element_rect(fill = NA))
  }) |>
  wrap_plots(nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
fig5C
ggsave(plot = fig5C, filename = paste0(out.dir, "/Figure5/tumourTC.png"))

# Sankey diagram
fig5D <- seuratobj@meta.data %>%
  select(ROI, ROI.subclone, ROI.tcs) %>%
  na.omit() %>%
  mutate(ROI.subclone = paste0("SC", ROI.subclone)) %>%
  make_long(ROI, ROI.subclone, ROI.tcs) %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node,
             fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE) +
  geom_sankey_label(size = 3, color = "black") +
  scale_fill_manual(values = colour.all) +
  theme_void() +
  theme(text = element_text(size = 18),
        legend.position = "none")
fig5D
ggsave(plot = fig5D, filename = paste0(out.dir, "/Figure5/sankey.png"))

# Edge and center ROI proportions (Fisher's exact test)
metadata <- seuratobj@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, ROI, edge, ROI.subclone, ROI.tcs)

stat.list <- metadata %>%
  na.omit() %>%
  add_count(ROI.subclone, ROI) %>%
  filter(n >= 100) %>%
  group_by(ROI.subclone, ROI, edge, ROI.tcs) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  spread(ROI.tcs, n) %>%
  mutate(TC1.1 = if_else(is.na(TC1.1), 0, TC1.1),
         TC2 = if_else(is.na(TC2), 0, TC2),
         TC3 = if_else(is.na(TC3), 0, TC3)) %>%
  group_split(ROI.subclone, ROI)

stat.test <- lapply(stat.list, FUN = function(x) {
  contingency <- x %>%
    select(-ROI.subclone, -ROI) %>%
    column_to_rownames("edge")
  x %>%
    mutate(pval = fisher_test(contingency, simulate.p.value = TRUE)$p) %>%
    select(ROI.subclone, ROI, pval) %>%
    unique()
}) |>
  bind_rows() %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  mutate(ROI.subclone = str_remove(ROI.subclone, pattern = "subclone"),
         ROI.ROI.subclone = paste(ROI,"\nSubclone", ROI.subclone),
         p.adj.signif = paste("p.adj =", pval.adj.signif))

fig5E <- metadata %>%
  mutate(ROI.ROI.subclone =  paste(ROI,"\nSubclone", ROI.subclone),
         edge = str_to_title(edge)) %>%
  filter(ROI.ROI.subclone %in% stat.test$ROI.ROI.subclone) %>%
  ggplot(aes(x = edge, fill = ROI.tcs)) +
  geom_bar(position = "fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) +
  scale_fill_manual(values = colour.tcs) +
  labs(x = NULL, y = "Proportion", fill = "Tumour TC") +
  facet_grid(~ROI.ROI.subclone, scales = "free_x") +
  pub_theme_facet +
  geom_text(data = stat.test, size = 3.5, inherit.aes = FALSE,
            aes(x = 1.5, y = 1.05, label = p.adj.signif, vjust = 0.5)) +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, vjust = 0.5),
        legend.position = "bottom")
fig5E
ggsave(plot = fig5E, paste0(out.dir, "/Figure5/fisher.png"))

# CNA
erbb2 <- cna %>%
  filter(gene_name == "ERBB2") %>%
  mutate(subclone = as.character(subclone),
         CN = as.character(CN)) %>%
  rename(ERBB2.CN = CN) %>%
  select(subclone, ERBB2.CN)

erbb2.metadata <- seuratobj@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, subclone) %>%
  left_join(erbb2, by = "subclone") %>%
  select(-subclone) %>%
  column_to_rownames("spots")

seuratobj <- AddMetaData(seuratobj, erbb2.metadata)

fig5F <-
  SpatialDimPlot(seuratobj, group.by = "ERBB2.CN", pt.size.factor = 1.5,
                 image.alpha = 0.7, combine = FALSE) |>
  lapply(FUN = function(x) {
    slide <- str_remove(x$labels$title, pattern = ": ERBB2")
    x + scale_fill_manual(values = colour.cna, na.value = NA) +
      ggtitle(paste("V19L29,", slide)) +
      labs(fill = "ERBB2 CNA") +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      theme(text = element_text(size = 18),
            legend.key = element_rect(fill = NA))
  }) |>
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect")
fig5F
ggsave(plot = fig5F, filename = paste0(out.dir, "/Figure5/ERBB2_CNA.png"))

# Biomarker expression
fig5G <-
  bcSignatures(bcobj, genes = list(values = "ERBB2", share.limits = TRUE),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) |>
  lapply(FUN = function(x) {
    x +  scale_fill_gradientn(colours = colourscale(n = 100, rev = TRUE),
                              limits = c(0, 3.2)) +
      ggtitle(NULL) +
      labs(fill = "ERBB2\nexpression") +
      theme(text = element_text(size = 18))
  }) |>
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect")
fig5G
ggsave(plot = fig5G, filename = paste0(out.dir, "/Figure5/ERBB2_biomarker.png"))

# Response to HER2 inhibitor and alterantive treatment
fig5H <-
  bcSignatures(bcobj,
               signatures = list(values = "BRD-K46386702-001-02-1_PRISM_K46386702"),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) |>
  lapply(FUN = function(x) {
    x + ggtitle(NULL)
  }) |>
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Varlitinib (PRISM) sensitivity",
                  subtitle = "HER2 inhibitor") &
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5))
fig5H
ggsave(plot = fig5H, 
       filename = paste0(out.dir, "/Figure5/varlitinib_sensitivity.png"))

fig5I <-
  bcSignatures(bcobj,
               signatures = list(values = "doxorubicin_PRISM_K92093830"),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) |>
  lapply(FUN = function(x) {
    x + ggtitle(NULL)
  }) |>
  wrap_plots(ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(title = "Doxorubicin (PRISM) sensitivity",
                  subtitle = "DNA related agent") &
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5))
fig5I
ggsave(plot = fig5I, 
       filename = paste0(out.dir, "/Figure5/doxorubicin_sensitivity.png"))

# Save figure
left <- (fig5A / fig5B / fig5C)
right <- (fig5D / fig5E) +
  plot_layout(heights = c(2, 1))
bottom <- (wrap_elements(fig5F) | wrap_plots(fig5G)) / 
               (wrap_elements(fig5H) | wrap_elements(fig5I)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

fig5 <- (((left | right) + plot_layout(widths = c(1.75, 2.25))) / bottom) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = fig5, width = 21, height = 29.7,
       filename = paste0(out.dir, "Figure5/Figure5.pdf"))
