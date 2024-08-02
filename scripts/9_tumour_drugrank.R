# conda activate pub-beyoncell
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
set.seed(1)
out.dir <- "results/tumour/"

# --- Data ---
# Merged tumour BCS with TCs
bctumour <- readRDS("results/tumour/sensitivity/tcs_tumour.rds")

# Paths to all Seurat objects
seuratobjs.path <-
  list.files("data/visium", pattern = ".rds", full.names = TRUE)
patients <- str_remove(basename(seuratobjs.path), pattern = "_.*$")
seuratobjs.path <- setNames(seuratobjs.path, patients)

# MoAs
moas <- read.table("data/signatures/MoAs_tumour.tsv", header = TRUE,
                   sep = "\t")

# Integrated Seurat object
seurat.integrated <- readRDS("results/all/integrated_seurat.rds")

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "sensitivity/"), recursive = TRUE,
           showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Molecular stratification of tumour spots
mol.strata <- GenerateGenesets("data/signatures/molecular_stratification.gmt")

strata.metadata <- lapply(seuratobjs.path, FUN = function(x) {
  seuratobj <- readRDS(x)
  bcobj <- bcScore(seuratobj, mol.strata, expr.thres = 0.1)
  bcobj@normalized[is.na(bcobj@normalized)] <- 0
  spots.tumour <- bcobj@meta.data %>%
    rownames_to_column("spots") %>%
    filter(Tumour == "Tumour") %>%
    pull(spots)
  mol.subtype <- bcobj@normalized[, spots.tumour] %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spots") %>%
    pivot_longer(cols = rownames(bcobj@normalized), names_to = "mol.subtype",
                 values_to = "enrichment") %>%
    mutate(mol.subtype = str_remove(mol.subtype, pattern = "_SC$"),
           mol.subtype = case_when(mol.subtype == "LUMA" ~ "LumA",
                                   mol.subtype == "LUMB" ~ "LumB",
                                   mol.subtype == "HER2E" ~ "Her2",
                                   TRUE ~ str_to_title(mol.subtype))) %>%
    group_by(spots) %>%
    filter(enrichment == max(enrichment)) %>%
    ungroup() %>%
    add_count(spots) %>%
    mutate(mol.subtype = if_else(n > 1, "Undetermined", mol.subtype)) %>%
    select(spots, mol.subtype) %>%
    unique()
  return(mol.subtype)
}) |>
  bind_rows() %>%
  column_to_rownames("spots")

bctumour <- bcAddMetadata(bctumour, strata.metadata)

# Retrieve the most specific drugs to target each tumour TC
bctumour <- bcRanks(bctumour, idents = "tumour.tcs",
                    resm.cutoff = c(0.01, 0.99), extended = FALSE)

diff.response <- c("TOP-Differential-HighSensitivity",
                   "TOP-Differential-LowSensitivity")

response <- bctumour@ranks$tumour.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("group"), names_to = "tumour_TC",
               values_to = "response", names_prefix = "group.") %>%
  filter(response %in% diff.response) %>%
  select(IDs, tumour_TC, response)

switch.point <- bctumour@ranks$tumour.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("switch.point"), names_to = "tumour_TC",
               values_to = "switch_point", names_prefix = "switch.point.") %>%
  select(IDs, tumour_TC, switch_point)

mean.bcs <- bctumour@ranks$tumour.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("mean"), names_to = "tumour_TC",
               values_to = "mean", names_prefix = "mean.") %>%
  select(IDs, tumour_TC, mean)

residuals.mean <- bctumour@ranks$tumour.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("residuals.mean"), names_to = "tumour_TC",
               values_to = "residuals_mean", names_prefix = "residuals.mean.") %>%
  select(IDs, tumour_TC, residuals_mean)

ranking <- response %>%
  left_join(switch.point, by = c("IDs", "tumour_TC")) %>%
  left_join(mean.bcs, by = c("IDs", "tumour_TC")) %>%
  left_join(residuals.mean, by = c("IDs", "tumour_TC")) %>%
  rename(drugID = IDs) %>%
  slice_max(order_by = residuals_mean, n = 15, by = tumour_TC)

# Add drug names and MoAs
ranking <- ranking %>%
  left_join(moas, by = "drugID", relationship = "many-to-many") %>%
  mutate(study = str_remove_all(str_extract(drugID, pattern= "_.*_"),
                                pattern = "_"),
        drug_name = paste0(drug_name, " (", study, ")")) %>%
  select(drugID, drug_name, MoA, tumour_TC, response, switch_point, mean,
         residuals_mean)

write.table(ranking, file = "tables/SuppTable12.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Save
saveRDS(bctumour, file = paste0(out.dir, "sensitivity/tcs_tumour.rds"))

# Add metadata and subset tumour spots from integrated Seurat object
new.metadata <- bctumour@meta.data %>%
  select(Tumour, tumour.tcs)
seurat.integrated <- AddMetaData(seurat.integrated, new.metadata)
seurat.tumour <- subset(seurat.integrated, Tumour == "Tumour")

# Log normalise counts
DefaultAssay(seurat.tumour) <- "Spatial"
seurat.tumour <- NormalizeData(seurat.tumour)
seurat.tumour <-
  FindVariableFeatures(seurat.tumour, selection.method = "vst",
                       nfeatures = 2000)
seurat.tumour <- ScaleData(seurat.tumour)

# Differential expression between tumour TCs (for GSEA pre-ranked)
Idents(seurat.tumour) <- "tumour.tcs"
tumour.tcs <- unique(seurat.tumour@meta.data$tumour.tcs)

dea <- lapply(tumour.tcs, FUN = function(x) {
  markers <- FindMarkers(seurat.tumour, ident.1 = x, min.pct = 0,
                         logfc.threshold = 0, test.use = "wilcox")
  markers <- markers %>%
    rownames_to_column("gene") %>%
    mutate(tumour.tcs = x)
  return(markers)
}) %>%
  bind_rows()

write.table(dea, file = paste0(out.dir, "DEA_tumourTCs.tsv"), sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
