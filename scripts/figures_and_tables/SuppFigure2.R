# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (1823-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
library(fgsea) # fgsea_1.26.0
library(ggseabubble) # ggseabubble_0.1.0
library(ComplexHeatmap) # ComplexHeatmap_2.16.0
library(circlize) # circlize_0.4.15
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
set.seed(1)
out.dir <- "figures/"

# --- Data ---
# Merged BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Differential gene expression result
dea <- read.table("results/all/DEA_majorTCs.tsv", header = TRUE, sep = "\t")

# Functional signatures
reactome <- readGMT("data/signatures/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmarks <- readGMT("data/signatures/h.all.v2023.1.Hs.symbols.gmt")
interesting <- readGMT("data/signatures/functional_signatures.gmt")
functional.sigs <- c(reactome, hallmarks, interesting)

# Drug ranking
drugrank <- read.table("tables/SuppTable8.tsv", header = TRUE, sep = "\t")

# Functional signature references
functional.refs <- read.table("tables/SuppTable6.tsv", header = TRUE,
                              sep = "\t")

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "SuppFigure2"), recursive = TRUE,
           showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Compute fgsea
comparisons <- unique(dea$major.tcs)
gsea <- lapply(comparisons, FUN = function(x) {
  markers <- dea %>%
    filter(major.tcs == x) %>%
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
  mutate(COMPARISON = if_else(COMPARISON == "Mixed", "Interphase",
                              COMPARISON)) %>%
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

write.table(supp.table, file = "tables/SuppTable7.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Get the top 15 signatures per major TC
chosen.fun.sigs <- gsea.table %>%
  filter(FDR.q.val < 0.25) %>%
  slice_max(order_by = NES, n = 15, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()

gsea.plot <- gsea.table %>%
  filter(NAME %in% chosen.fun.sigs) %>%
  mutate(COMPARISON = factor(COMPARISON,
                             levels = c("Tumour-rich", "Interphase",
                                        "TME-rich")))

# Plot bubbleheatmap
bubbleheat <-
  ggbubbleHeatmap(gsea.plot, cluster.cols = FALSE, FDR.threshold = 0.25,
                  bubble.size.range = c(4, 10))

bubbleheat <- bubbleheat + theme(text = element_text(size = 18))

pdf(paste0(out.dir, "SuppFigure2/functional_bubbleheatmap.pdf"), width = 19.78,
    height = 8.47)
bubbleheat
dev.off()

figSupp2A <- fig(paste0(out.dir, "SuppFigure2/functional_bubbleheatmap.pdf"))

# Format figure MoAs
moas <- drugrank %>%
  mutate(drug.name = toupper(drug_name)) %>%
  select(drugID, drug.name, MoA) %>%
  unique() %>%
  as.data.frame()

rownames(moas) <- moas$drugID

# Chosen drugs
chosen.drugs <- drugrank %>%
  pull(drugID) %>%
  unique()

# Get metadata
metadata <- bcmerged@meta.data %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs =
           factor(major.tcs,
                  levels = c("Tumour-rich", "Interphase", "TME-rich")),
         subtype = if_else(subtype == "ER", "Luminal", subtype),
         subtype = factor(subtype, levels = c("Luminal", "TNBC", "HER2+"))) %>%
  arrange(as.numeric(major.tcs), main.cell.type, subtype)

# Heatmap of different response patterns
scores <- bcmerged@normalized
col.anno <- HeatmapAnnotation("Major TC" = metadata$major.tcs,
                              "Cell type" = metadata$main.cell.type,
                              "Subtype" = metadata$subtype,
                              col = list("Major TC" = colour.tctypes,
                                         "Cell type" = colour.celltypes,
                                         "Subtype" = colour.subtype))
row.anno <- rowAnnotation("MoA" = moas[chosen.drugs, "MoA"],
                          col = list("MoA" = colour.moas))

heat <- Heatmap(scores[chosen.drugs, metadata$spots, drop = FALSE],
                name = "Normalised\nBCS", cluster_columns = FALSE,
                row_dend_reorder = TRUE,
                top_annotation = col.anno, right_annotation = row.anno,
                column_split = metadata$major.tcs, row_split = 2,
                column_title = NULL, row_title = NULL,
                column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
                clustering_method_rows = "ward.D2",
                show_column_names = FALSE,
                row_labels = moas[chosen.drugs, "drug.name"])

pdf(paste0(out.dir, "SuppFigure2/sensitivity_heatmap.pdf"), width = 20,
    height = 15)
heat
dev.off()

figSupp2B <- fig(paste0(out.dir, "SuppFigure2/sensitivity_heatmap.pdf"))

# Mean BCS per group
drugrank <- drugrank %>%
  select(drugID, major_TC) %>%
  rename(signature = drugID, response = major_TC)

scores[chosen.drugs, ] %>%
  as.data.frame() %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = colnames(scores), names_to = "spots",
               values_to = "enrichment") %>%
  left_join(metadata[, c("spots", "major.tcs")], by = "spots") %>%
  left_join(drugrank, by = "signature") %>%
  group_by(major.tcs, response) %>%
  summarise(mean.enrichment = mean(enrichment),
            sd.enrichment = sd(enrichment)) %>%
  arrange(response, major.tcs)

# Save figure
figSupp2 <- (figSupp2A / figSupp2B) +
  plot_layout(widths = c(4, 5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = figSupp2, width = 21, height = 29.7,
       filename = paste0(out.dir, "SuppFigure2/SuppFigure2.pdf"))
