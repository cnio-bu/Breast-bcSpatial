# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(corrplot) # corrplot_0.92
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
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
# Merged sensitivity BCS with TCs
bcsensitivity <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Merged functional BCS
bcfunctional <- readRDS("results/all/functional/functional_bcs_regressed.rds")

# Functional signatures
interesting <- readGMT("data/signatures/functional_signatures.gmt")

# Drug ranking
drugrank <- read.table("tables/SuppTable9.tsv", header = TRUE, sep = "\t")

# Radial distances to the tumour core
radial.dist <- read.table("results/all/radial_distance.tsv", header = TRUE,
                          sep = "\t")

# Sensitivity enrichment correlation with distances
sensitivity.corr <-
  read.table("results/all/sensitivity/sensitivity_correlation.tsv",
             header = TRUE, sep = "\t")

# Functional enrichment correlation with distances
functional.corr <-
  read.table("results/all/functional/functional_correlation.tsv",
             header = TRUE, sep = "\t")

# Functional signature references
functional.refs <- read.table("tables/SuppTable7.tsv", header = TRUE,
                              sep = "\t")

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "Figure3"), recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)
dir.create("tables/chosenIDs", recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Get metadata
metadata <- bcsensitivity@meta.data %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs =
           factor(major.tcs,
                  levels = c("Tumour-rich", "Interphase", "TME-rich")),
         patient = factor(patient, levels = patient.order))

# Cell type proportions in each compartment, per patient
fig3A <- metadata %>%
  ggplot(aes(x = major.tcs, fill = main.cell.type)) +
  geom_bar(position = "fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  facet_grid(~patient) +
  scale_fill_manual(values = colour.celltypes) +
  labs(y = "Proportion", fill = "Cell type") +
  pub_theme_facet +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 18))
fig3A
ggsave(plot = fig3A, height = 5,
       filename = paste0(out.dir, "Figure3/cell_proportions.png"))

# Example
patient <- "CID44971"
subtype <- "TNBC"
bcobj.sensitivity <-
  readRDS(paste0("data/beyondcell/sensitivity/", patient,
                 "_spatial_bcs_sensitivity.rds"))
bcobj.functional <-
  readRDS(paste0("data/beyondcell/functional/", patient,
                 "_spatial_bcs_functional.rds"))
seuratobj <-
  readRDS(paste0("data/visium/", patient, "_spatial_SCTnormalised_clones.rds"))

spots.patient <- colnames(bcobj.sensitivity@normalized)
dist <- radial.dist %>%
  filter(spots %in% spots.patient) %>%
  select(spots, r_dist_tumour_squared) %>%
  column_to_rownames("spots")
bcobj.sensitivity <- bcAddMetadata(bcobj.sensitivity, dist)

bcobj.sensitivity@meta.data <- bcobj.sensitivity@meta.data %>%
  mutate(r_dist_tumour_squared = if_else(region == "A", NA_real_, r_dist_tumour_squared))

gs.biomarkers <- GenerateGenesets("data/signatures/biomarkers.gmt")
bcobj.biomarkers <- bcScore(seuratobj, gs.biomarkers)
bcobj.biomarkers@normalized[is.na(bcobj.biomarkers@normalized)] <- 0
bcobj.biomarkers <- bcRecompute(bcobj.biomarkers, slot = "normalized")
biomarkers <-
  data.frame(t(bcobj.biomarkers@scaled["biomarkers", , drop = FALSE]))
seuratobj <- AddMetaData(seuratobj, biomarkers)

cid44971.biomarker <-
  SpatialFeaturePlot(seuratobj, features = "biomarkers", pt.size.factor = 1.7,
                     image.alpha = 0.7) +
  ggtitle("ESR1 + PGR + ERBB2") +
  theme(text = element_text(size = 18), legend.position = "right",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"))

cid44971.doxorubicin <-
  bcSignatures(bcobj.sensitivity,
               signatures = list(values = "doxorubicin_PRISM_K92093830"),
               pt.size = 1.7, spatial = TRUE, image.alpha = 0.7) +
  ggtitle("Allitinib (PRISM) sensitivity", subtitle = "HER2 inhibitor") +
  theme(text = element_text(size = 18), legend.position = "bottom",
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5))

cid44971.cellcycle <-
  bcSignatures(bcobj.functional,
               signatures = list(values = "sig_cell_cycle_CancerSEA"),
               pt.size = 1.7, spatial = TRUE, image.alpha = 0.7) +
  ggtitle("EMT (Groger et al. 2012) enrichment") +
  theme(text = element_text(size = 18), legend.position = "bottom",
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))

cid44971.dist <-
  bcClusters(bcobj.sensitivity, idents = "r_dist_tumour_squared",
             factor.col = FALSE, pt.size = 1.7, spatial = TRUE,
             image.alpha = 0.7) +
  scale_fill_gradientn(colours = colourscale(n = 100, rev = TRUE),
                       values = scales::rescale(c(min(dist$r_dist_tumour), 0,
                                                  max(dist$r_dist_tumour))),
                       na.value = NA) +
  labs(fill = NULL) +
  ggtitle("Radial distance") +
  theme(text = element_text(size = 18), legend.position = "bottom",
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"))

fig3B <- cid44971.biomarker / cid44971.doxorubicin / cid44971.cellcycle /
  cid44971.dist
fig3B
ggsave(plot = fig3B, width = 8, height = 16,
       filename = paste0(out.dir, "Figure3/example.png"))

# Distance correlation with biological functions
interesting.sigs <-
  str_remove(names(interesting), pattern = "(_UP$)|(_DOWN$)") %>%
  unique()

functional.refs <- functional.refs %>%
  filter(!sigID %in% paste0("sig_EMT_GROGER_2012_", c("UP", "DOWN"))) %>%
  filter(!Signature %in%
    paste(c("DUCTAL", "INVASIVE"), "BREAST CANCER (Schuetz et al. 2006)")) %>%
  mutate(sigID = str_remove(sigID, pattern = "(_UP$)|(_DOWN$)"))

functional.corr <- functional.corr %>%
  filter(signature %in% interesting.sigs) %>%
  rename(sigID = signature) %>%
  left_join(functional.refs, by = "sigID") %>%
  select(-sigID) %>%
  rename(signature = Signature)

functional.supp.table <- functional.corr %>%
  rename(Signature = signature, Patient = patient, Slide = slide,
         Region = region, Pearson.correlation = cor, Statistic = statistic,
         pval = p, FDRpval = p.adj) %>%
  mutate(Pearson.correlation = round(Pearson.correlation, digits = 2),
         Statistic = round(Statistic, digits = 2),
         pval = formatC(pval, format = "e", digits = 2),
         FDRpval = formatC(FDRpval, format = "e", digits = 2),
         conf.low = round(conf.low, digits = 2),
         conf.high = round(conf.high, digits = 2),
         CI95 = paste(conf.low, conf.high, sep = ",")) %>%
  select(Signature, Patient, Slide, Region, Pearson.correlation, Statistic,
         pval, FDRpval, CI95)

write.table(functional.supp.table, file = "tables/SuppTable5.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Create correlation matrix
functional.corr.wide <- functional.corr %>%
  mutate(patient = factor(patient, levels = patient.order)) %>%
  arrange(as.numeric(patient)) %>%
  mutate(ID = paste(patient, slide, region, sep = "_")) %>%
  filter(!ID %in% c("CID44971_slide1_A", "CID4535_slide_A")) %>%
  pivot_wider(id_cols = "ID", names_from = "signature",
              values_from = "cor") %>%
  column_to_rownames("ID")

functional.pval.wide <- functional.corr %>%
  mutate(patient = factor(patient, levels = patient.order)) %>%
  arrange(as.numeric(patient)) %>%
  mutate(ID = paste(patient, slide, region, sep = "_")) %>%
  filter(!ID %in% c("CID44971_slide1_A", "CID4535_slide_A")) %>%
  pivot_wider(id_cols = "ID", names_from = "signature",
              values_from = "p.adj") %>%
  column_to_rownames("ID")

## Compute correlation mean per patient
functional.corr.mean <- colMeans(functional.corr.wide)
write.table(names(functional.corr.mean), file = "tables/chosenIDs/pathways_Figure3C.txt",
            sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Patient subtype
subtypes <- bcfunctional@meta.data %>%
  rownames_to_column("spots") %>%
  select(patient, subtype) %>%
  mutate(subtype = if_else(subtype == "ER", "Luminal", subtype)) %>%
  select(patient, subtype) %>%
  unique() %>%
  arrange(patient, subtype) %>%
  column_to_rownames("patient")

functional.corr.subtype <-
  subtypes[str_remove(rownames(functional.corr.wide), pattern = "_.*"), ]

## Enrichment correlation with distance
functional.colour.labels <- colour.subtype[functional.corr.subtype]

pdf(width = 10, height = 7,
    file = paste0(out.dir, "Figure3/functional_correlation.pdf"))
corrplot(t(as.matrix(functional.corr.wide[, order(functional.corr.mean)])),
         p.mat = t(as.matrix(functional.pval.wide[, order(functional.corr.mean)])),
         sig.level = 0.05, insig = "label_sig", pch.cex = 1,
         col = colorRampPalette(c("#FF5D5DFF", "white", "#45B7C4"))(200),
         tl.col = functional.colour.labels[functional.corr.subtype],
         tl.cex = 0.6, tl.srt = 45)
dev.off()

fig3C <- fig(paste0(out.dir, "Figure3/functional_correlation.pdf"))

# Distance correlation with drug sensitivities
chosen.moas <- c("Hormonal therapy", "HER2 inhibitor",
                 "WNT signaling inhibitor", "VEGFR inhibitor",
                 "TGFBR inhibitor", "NFkB signaling inhibitor",
                 "JAK-STAT signaling inhibitor",
                 "Anti-inflammatory / Immunosuppressant", "Statin",
                 "BRAF inhibitor", "EGFR inhibitor")
other.drugs <- c("Olaparib_GDSC_1017", "piperlongumine:MST-312_CTRP_M37545453",
                 "MST-312_CTRP_K19894101", "BI-2536_CTRP_K64890080",
                 "BI-2536_GDSC_1086", "KW-2449_CTRP_K21718444",
                 "PI-103_CTRP_K67868012", "PF-4708671_GDSC_1129",
                 "Trametinib_GDSC_1372", "Ulixertinib_GDSC_1908",
                 "Selumetinib_GDSC_1736", "trifluridine_PRISM_K03243820",
                 "Mitoxantrone_GDSC_1810", "Oxaliplatin_GDSC_1806",
                 "Teniposide_GDSC_180")
chosen.drugs <- drugrank %>%
  filter(MoA %in% chosen.moas | drugID %in% other.drugs) %>%
  pull(drugID) %>%
  unique()

sensitivity.corr <- sensitivity.corr %>%
  filter(signature %in% chosen.drugs)

## Format figure MoAs
moas <- drugrank %>%
  filter(drugID %in% chosen.drugs) %>%
  mutate(drug.name = toupper(drug_name)) %>%
  select(drugID, drug.name, MoA) %>%
  unique() %>%
  as.data.frame()

rownames(moas) <- moas$drugID

sensitivity.supp.table <- sensitivity.corr %>%
  rename(drugID = signature, Patient = patient, Slide = slide,
         Region = region, Pearson.correlation = cor, Statistic = statistic,
         pval = p, FDRpval = p.adj) %>%
  mutate(Pearson.correlation = round(Pearson.correlation, digits = 2),
         Statistic = round(Statistic, digits = 2),
         pval = formatC(pval, format = "e", digits = 2),
         FDRpval = formatC(FDRpval, format = "e", digits = 2),
         conf.low = round(conf.low, digits = 2),
         conf.high = round(conf.high, digits = 2),
         CI95 = paste(conf.low, conf.high, sep = ",")) %>%
  left_join(moas, by = "drugID") %>%
  select(drugID, drug.name, MoA, Patient, Slide, Region, Pearson.correlation,
         Statistic, pval, FDRpval, CI95)

write.table(sensitivity.supp.table, file = "tables/SuppTable6.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

## Create correlation matrix
sensitivity.corr.wide <- sensitivity.corr %>%
  mutate(patient = factor(patient, levels = patient.order)) %>%
  arrange(as.numeric(patient)) %>%
  mutate(ID = paste(patient, slide, region, sep = "_")) %>%
  filter(!ID %in% c("CID44971_slide1_A", "CID4535_slide_A")) %>%
  pivot_wider(id_cols = "ID", names_from = "signature",
              values_from = "cor") %>%
  column_to_rownames("ID")

sensitivity.pval.wide <- sensitivity.corr %>%
  mutate(patient = factor(patient, levels = patient.order)) %>%
  arrange(as.numeric(patient)) %>%
  mutate(ID = paste(patient, slide, region, sep = "_")) %>%
  filter(!ID %in% c("CID44971_slide1_A", "CID4535_slide_A")) %>%
  pivot_wider(id_cols = "ID", names_from = "signature",
              values_from = "p.adj") %>%
  column_to_rownames("ID")

## Compute correlation mean per patient
sensitivity.corr.mean <- colMeans(sensitivity.corr.wide)
write.table(names(sensitivity.corr.mean), file = "tables/chosenIDs/drugs_Figure3C.txt",
            sep = "\n", col.names = FALSE, row.names = FALSE, quote = FALSE)

## Patient subtype
sensitivity.corr.subtype <-
  subtypes[str_remove(rownames(sensitivity.corr.wide), pattern = "_.*"), ]

## Enrichment correlation with distance
sensitivity.colour.labels <-
  colour.moas[moas[colnames(sensitivity.corr.wide), "MoA"]]

colnames(sensitivity.corr.wide) <-
  moas[colnames(sensitivity.corr.wide), "drug.name"]
names(sensitivity.corr.mean) <-
  moas[names(sensitivity.corr.mean), "drug.name"]
colnames(sensitivity.pval.wide) <-
  moas[colnames(sensitivity.pval.wide), "drug.name"]

pdf(width = 10, height = 15,
    file = paste0(out.dir, "Figure3/sensitivity_correlation.pdf"))
corrplot(t(as.matrix(sensitivity.corr.wide[, order(sensitivity.corr.mean)])),
         p.mat = t(as.matrix(sensitivity.pval.wide[, order(sensitivity.corr.mean)])),
         sig.level = 0.05, insig = "label_sig", pch.cex = 1,
         col = colorRampPalette(c("#FF5D5DFF", "white", "#45B7C4"))(200),
         tl.col = sensitivity.colour.labels[order(sensitivity.corr.mean)],
         tl.cex = 0.6, tl.srt = 45)
dev.off()

fig3D <- fig(paste0(out.dir, "Figure3/sensitivity_correlation.pdf"))

# Save figure
bottom.right <- (fig3C / fig3D) +
  plot_layout(heights = c(3, 6))
bottom <- (fig3B | bottom.right) +
  plot_layout(widths = c(1, 2))

fig3 <- (fig3A / bottom) +
  plot_layout(heights = c(1, 7)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = fig3, width = 21, height = 29.7,
       filename = paste0(out.dir, "Figure3/Figure3.pdf"))
