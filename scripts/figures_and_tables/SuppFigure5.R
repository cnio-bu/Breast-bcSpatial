# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(fgsea) # fgsea_1.26.0
library(ggseabubble) # ggseabubble_0.1.0
library(beyondcell) # beyondcell_2.2.0
library(RColorBrewer) # RColorBrewer_1.1-3
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
set.seed(1)
out.dir <- "figures/"

# --- Data ---
# Seurat object with ROI annotation
seuratobj <-
  readRDS("results/tumour/ROI/V19L29_spatial_SCTnormalised_clones_ROI.rds")

# Merged BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Differential gene expression result
dea <- read.table("results/tumour/ROI/DEA_ROIs.tsv", header = TRUE, sep = "\t")

# Signatures
reactome <- readGMT("data/signatures/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmarks <- readGMT("data/signatures/h.all.v2023.1.Hs.symbols.gmt")
interesting <- readGMT("data/signatures/functional_signatures.gmt")
functional.sigs <- c(reactome, hallmarks, interesting)

# Functional signature references
functional.refs <- read.table("tables/SuppTable6.tsv", header = TRUE,
                              sep = "\t")

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "SuppFigure5"), recursive = TRUE,
           showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# H&E image
figSupp5A <- SpatialDimPlot(seuratobj, group.by = "ROI", pt.size = 1.5,
                            image.alpha = 0.7, combine = FALSE) |>
  lapply(FUN = function(x) {
    title <- str_remove(x$labels$title, pattern = "slide")
    x + scale_fill_manual(values = colour.annotations, na.value = NA) +
      ggtitle(paste0("V19L29, slide", title)) +
      theme(text = element_text(size = 18), legend.position = "none",
            plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))
  }) |>
  wrap_plots(nrow = 1)
figSupp5A
ggsave(plot = figSupp5A,
       filename = paste0(out.dir, "SuppFigure5/HE_image.png"))

# Expression clusters
figSupp5B <- SpatialDimPlot(seuratobj, group.by = "sub.cluster", pt.size = 1.5,
                            image.alpha = 0.7, combine = FALSE) |>
  lapply(FUN = function(x) {
    x + scale_fill_manual(values = colour.expression, na.value = NA) +
      ggtitle(NULL) +
      labs(fill = "Expression\nclusters") +
      guides(fill = guide_legend(override.aes = list(size = 3))) +
      theme(text = element_text(size = 18),
            legend.key = element_rect(fill = NA))
  }) |>
  wrap_plots(nrow = 1) +
  plot_layout(guides = "collect")
figSupp5B
ggsave(plot = figSupp5B,
       filename = paste0(out.dir, "SuppFigure5/expression_clusters.png"))

# Compute fgsea
comparisons <- unique(dea$ROI)
gsea <- lapply(comparisons, FUN = function(x) {
  markers <- dea %>%
    filter(ROI == x) %>%
    arrange(desc(avg_log2FC))
  rnk <- setNames(markers$avg_log2FC, markers$gene)
  res <- fgsea(pathways = functional.sigs, stats = rnk, minSize = 15,
               maxSize = 500, nPermSimple = 10000) %>%
    mutate(COMPARISON = x)
})

# Format table
gsea.table <- gsea %>%
  bind_rows() %>%
  rename(sigID = pathway) %>%
  left_join(functional.refs, by = "sigID") %>%
  select(-sigID, -Reference) %>%
  rename(NAME = Signature,
         FDR.q.val = padj) %>%
  unnest(leadingEdge) %>%
  group_by(NAME, COMPARISON) %>%
  mutate(leadingEdge = paste(leadingEdge, collapse = ",")) %>%
  ungroup() %>%
  unique() %>%
  filter(!is.na(COMPARISON))

supp.table <- gsea.table %>%
  rename(Signature = NAME,
         Comparison = COMPARISON) %>%
  mutate(Comparison = paste(Comparison, "vs rest")) %>%
  select(Signature, Comparison, pval:leadingEdge)

write.table(supp.table, file = "tables/SuppTable12.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Get the top 15 signatures per ROI
chosen.fun.sigs <- gsea.table %>%
  filter(FDR.q.val < 0.25) %>%
  slice_max(order_by = NES, n = 15, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()

gsea.table <- gsea.table %>%
  filter(NAME %in% chosen.fun.sigs) %>%
  mutate(COMPARISON = factor(COMPARISON))

# Plot bubbleheatmap
bubbleheat <-
  ggbubbleHeatmap(gsea.table, cluster.cols = FALSE, FDR.threshold = 0.25,
                  bubble.size.range = c(4, 10))
bubbleheat <- bubbleheat + theme(text = element_text(size = 18))

pdf(paste0(out.dir, "SuppFigure5/bubble_functional.pdf"), width = 18.78,
    height = 9.47)
bubbleheat
dev.off()

figSupp5C <- fig(paste0(out.dir, "SuppFigure5/bubble_functional.pdf"))

# Save figure
figSupp5 <- figSupp5A / figSupp5B / figSupp5C +
  plot_layout(heights = c(2, 2, 5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = figSupp5, width = 21, height = 29.7,
       filename = paste0(out.dir, "SuppFigure5/SuppFigure5.pdf"))
