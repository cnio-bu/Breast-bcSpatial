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
# Tumour sensitivity BCS with TCs
bctumour <- readRDS("results/tumour/sensitivity/tcs_tumour.rds")

# Differential gene expression result
dea <- read.table("results/tumour/DEA_tumourTCs.tsv", header = TRUE, sep = "\t")

# Functional signatures
reactome <- readGMT("data/signatures/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmarks <- readGMT("data/signatures/h.all.v2023.1.Hs.symbols.gmt")
interesting <- readGMT("data/signatures/functional_signatures.gmt")
functional.sigs <- c(reactome, hallmarks, interesting)

# Drug ranking
drugrank <- read.table("tables/SuppTable12.tsv", header = TRUE, sep = "\t")

# Functional signature references
functional.refs <- read.table("tables/SuppTable7.tsv", header = TRUE,
                              sep = "\t")

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "Figure4"), recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/chosenIDs", recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# UMAP with initial TCs
fig4A <- bcClusters(bctumour, idents = "initial.tumour.tcs") +
  scale_color_manual(values = colour.initial.tcs) +
  theme_classic() +
  theme(text = element_text(size = 18))
fig4A
ggsave(plot = fig4A, width = 10,
       filename = paste0(out.dir, "Figure4/UMAPinitial.png"))

# Correlation heatmap
fig4B <- fig(paste0(out.dir, "Figure4/correlation_drugs.pdf"))

# UMAP with final TCs
fig4C <- bcClusters(bctumour, idents = "tumour.tcs") +
  scale_color_manual(values = colour.tcs) +
  theme_classic() +
  theme(text = element_text(size = 18))
fig4C
ggsave(plot = fig4C, width = 10,
       filename = paste0(out.dir, "Figure4/UMAPfinal.png"))

# Compute fgsea
comparisons <- unique(dea$tumour.tcs)
gsea <- lapply(comparisons, FUN = function(x) {
  markers <- dea %>%
    filter(tumour.tcs == x) %>%
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
  mutate(Comparison = paste(Comparison, "vs rest"),
         pval = formatC(pval, format = "e", digits = 2),
         FDR.q.val = formatC(FDR.q.val, format = "e", digits = 2),
         log2err = round(log2err, digits = 2),
         ES = round(ES, digits = 2),
         NES = round(NES, digits = 2)) %>%
  select(Signature, Comparison, pval:leadingEdge)

write.table(supp.table, file = "tables/SuppTable11.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Get the top 15 signatures per major TC
chosen.fun.sigs <- gsea.table %>%
  filter(FDR.q.val < 0.25) %>%
  slice_max(order_by = NES, n = 15, by = COMPARISON) %>%
  pull(NAME) %>%
  unique()

# Get the figure signatures
figure.sigs <-
  c("HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION_UP", "HALLMARK_COMPLEMENT",
    "REACTOME_INTERFERON_GAMMA_SIGNALING_UP", "REACTOME_PD_1_SIGNALING_UP",
    "HALLMARK_DNA_REPAIR", "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING_UP",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "REACTOME_CELLULAR_RESPONSE_TO_STARVATION_UP",
    "REACTOME_REGULATION_OF_EXPRESSION_OF_SLITS_AND_ROBOS_UP",
    "REACTOME_NONSENSE_MEDIATED_DECAY_NMD_UP",
    "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION_UP",
    "REACTOME_COLLAGEN_FORMATION_UP", "sig_invasion_CancerSEA_UP",
    "SCHUETZ_BREAST_CANCER_DUCTAL_VS_INVASIVE_DOWN",
    "REACTOME_FCGR_ACTIVATION_UP", "REACTOME_CD22_MEDIATED_BCR_REGULATION_UP",
    "REACTOME_FCERI_MEDIATED_MAPK_ACTIVATION_UP",
    "REACTOME_ROLE_OF_PHOSPHOLIPIDS_IN_PHAGOCYTOSIS_UP",
    "REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_LEADING_TO_GENERATION_OF_SECOND_MESSENGERS_UP")

chosen.fun.sigs <- functional.refs %>%
  filter(sigID %in% figure.sigs) %>%
  pull(Signature)
write.table(chosen.fun.sigs, file = "tables/chosenIDs/pathways_Figure4D.txt",
            sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

gsea.plot <- gsea.table %>%
  filter(NAME %in% chosen.fun.sigs) %>%
  mutate(COMPARISON = factor(COMPARISON))

# Plot bubbleheatmap
bubbleheat <-
  ggbubbleHeatmap(gsea.plot, cluster.cols = FALSE, FDR.threshold = 0.25,
                  bubble.size.range = c(4, 10))

bubbleheat <- bubbleheat + theme(text = element_text(size = 18))

pdf(paste0(out.dir, "Figure4/functional_bubbleheatmap.pdf"), width = 19.78,
    height = 8.47)
bubbleheat
dev.off()

fig4D <- fig(paste0(out.dir, "Figure4/functional_bubbleheatmap.pdf"))

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
write.table(chosen.drugs, file = "tables/chosenIDs/drugs_Figure4E.txt",
            sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

# Get metadata
metadata <- bctumour@meta.data %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs =
           factor(major.tcs,
                  levels = c("Tumour-rich", "Interphase", "TME-rich")),
         subtype = if_else(subtype == "ER", "Luminal", subtype),
         subtype = factor(subtype, levels = c("Luminal", "TNBC", "HER2+")),
         mol.subtype = if_else(mol.subtype == "Her2", "HER2+", mol.subtype),
         mol.subtype =
           factor(mol.subtype, levels = c("LumA", "LumB", "HER2+", "Basal"))) %>%
  arrange(as.numeric(major.tcs), main.cell.type, subtype)

# Heatmap of different response patterns
scores <- bctumour@normalized
col.anno <- HeatmapAnnotation("Tumour TC" = metadata$tumour.tcs,
                              "Major TC" = metadata$major.tcs,
                              "Intrinsic subtype" = metadata$mol.subtype,
                              col = list("Tumour TC" = colour.tcs,
                                         "Major TC" = colour.tctypes,
                                         "Intrinsic subtype" = colour.molsubtypes))
row.anno <- rowAnnotation("MoA" = moas[chosen.drugs, "MoA"],
                          col = list("MoA" = colour.moas))

heat <- Heatmap(scores[chosen.drugs, metadata$spots, drop = FALSE],
                name = "Normalised\nBCS", cluster_columns = FALSE,
                row_dend_reorder = TRUE,
                top_annotation = col.anno, right_annotation = row.anno,
                column_split = metadata$tumour.tcs, row_split = 5,
                column_title = NULL, row_title = NULL,
                column_gap = unit(1, "mm"), row_gap = unit(1, "mm"),
                clustering_method_rows = "ward.D2",
                show_column_names = FALSE,
                row_labels = moas[chosen.drugs, "drug.name"])

pdf(paste0(out.dir, "Figure4/sensitivity_heatmap.pdf"), width = 20,
    height = 15)
heat
dev.off()

fig4E <- fig(paste0(out.dir, "Figure4/sensitivity_heatmap.pdf"))

# Mean BCS per group
drugrank <- drugrank %>%
  filter(response == "TOP-Differential-HighSensitivity") %>%
  mutate(tumour_TC = str_remove(tumour_TC, pattern = "\\..*$")) %>%
  select(drugID, tumour_TC) %>%
  rename(signature = drugID, response = tumour_TC)

scores[chosen.drugs, ] %>%
  as.data.frame() %>%
  rownames_to_column("signature") %>%
  pivot_longer(cols = colnames(scores), names_to = "spots",
               values_to = "enrichment") %>%
  left_join(metadata[, c("spots", "tumour.tcs")], by = "spots") %>%
  inner_join(drugrank, by = "signature") %>%
  group_by(tumour.tcs, response) %>%
  summarise(mean.enrichment = median(enrichment),
            sd.enrichment = sd(enrichment)) %>%
  arrange(response, tumour.tcs)

# Save figure
top <- (fig4A | fig4B | fig4C) +
  plot_layout(widths = c(1, 2, 1))
fig4 <- (top / fig4D / fig4E) +
  plot_layout(heights = c(2, 2.5, 4)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = fig4, width = 21, height = 29.7,
       filename = paste0(out.dir, "Figure4/Figure4.pdf"))
