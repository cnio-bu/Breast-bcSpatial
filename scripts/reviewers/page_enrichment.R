# conda activate pub-giotto
rm(list = ls()) # R version 4.3.3 (2024-02-29)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(Giotto) # Giotto_4.1.0
library(rstatix) # rstatix_0.7.2
library(patchwork) # patchwork_1.2.0
set.seed(1)
out.dir <- "figures/reviewers/"

# --- Data ---
# Merged sensitivity BCS with TCs
bcsensitivity <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Paths to all Seurat objects
seuratobjs.path <-
  list.files("data/visium", pattern = ".rds", full.names = TRUE)

# --- Code ---
# SSc breast
ssc.breast <- GenerateGenesets("data/signatures/drug_signatures_classic_nodup.gmt")
ssc.breast.list <- readGMT("data/signatures/drug_signatures_classic_nodup.gmt")

# Patient IDs
patient.id <- bcsensitivity@meta.data %>%
  select(patient, slide) %>%
  unique() %>%
  mutate(ID = paste(patient, slide, sep = "/")) %>%
  pull(ID)

# Compute PAGE and BCS enrichments for each patient and slide
BCS <- PAGE <- metadata <- vector(mode = "list", length = length(patient.id))

for (i in 1:length(patient.id)) {
  # Patient and slide
  p <- str_remove(patient.id[i], pattern = "/.*$")
  s <- str_remove(patient.id[i], pattern = "^.*/")
  # Read Seurat object
  seuratobj <- readRDS(seuratobjs.path[min(i, 9)])
  # Split Seurat object
  seuratobjs <- SplitObject(seuratobj, split.by = "slide")
  # Select desired slide
  seuratobj.splitted <- seuratobjs[[s]]
  seuratobj.splitted@images <- seuratobj.splitted@images[s]
  # Convert to Giotto object
  giottoobj <- seuratToGiottoV4(seuratobj.splitted, spatial_assay = "Spatial", dim_reduction = NULL)
  # Normalise
  giottoobj <- normalizeGiotto(gobject = giottoobj, verbose = TRUE)
  # PAGE enrichment
  splits <- seq(1, length(ssc.breast.list), by = 50)
  splits <- lapply(1:length(splits), FUN = function(i) {
    c(splits[i], min(na.omit(splits[i+1]-1), length(ssc.breast.list)))
  })
  PAGE[[i]] <- lapply(splits, FUN = function(x) {
    sig.matrix <- makeSignMatrixPAGE(sign_names = names(ssc.breast.list)[x[1]:x[2]],
                                     sign_list = ssc.breast.list[x[1]:x[2]])
    runPAGEEnrich(gobject = giottoobj, sign_matrix = sig.matrix, min_overlap_genes = 25,
                  output_enrichment = "zscore", return_gobject = FALSE)$DT
  }) |>
    bind_rows() %>%
    mutate(patient = p, slide = s)
  # Normalised matrix and BCS
  normalized.m <- giottoobj@expression$cell$rna$normalized@exprMat |>
    as.matrix()
  BCS[[i]] <- bcScore(normalized.m, ssc.breast)@normalized %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spots") %>%
    mutate(patient = p, slide = s)
  # Metadata
  metadata[[i]] <- giottoobj@cell_metadata$cell$rna@metaDT %>%
    select(cell_ID) %>%
    rename(spots = cell_ID) %>%
    mutate(patient = p, slide = s)
}

# Merge lists
BCS <- bind_rows(BCS)
PAGE <- bind_rows(PAGE)
metadata <- bind_rows(metadata)

# Operate PAGE enrichments
PAGE <- PAGE %>%
  mutate(sigID = str_remove(cell_type, pattern = "(_UP)|(_DOWN)"),
         mode = str_remove(str_extract(cell_type, pattern = "(_UP)|(_DOWN)"),
                           pattern = "_")) %>%
  rename(spots = cell_ID) %>%
  select(spots, patient, slide, sigID, mode, zscore) %>%
  pivot_wider(id_cols = c("spots", "patient", "slide", "sigID"),
              names_from = mode, values_from = zscore) %>%
  group_by(spots, patient, slide, sigID) %>%
  summarise(zscore = UP - DOWN)

# Pivot longer BCS
BCS <- BCS %>%
  pivot_longer(cols = rownames(bcsensitivity@normalized), names_to = "sigID",
               values_to = "BCS")

# All enrichments
enrich <- BCS %>%
  left_join(PAGE, by = c("spots", "patient", "slide", "sigID"))

# Spearman correlation between PAGE enrichment and BCS
sensitivity.corr <- enrich %>%
  group_by(sigID, patient) %>%
  cor_test(vars = zscore, vars2 = BCS, method = "spearman") %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup()

# Plot
fig.page.bcs.corr <- sensitivity.corr %>%
  filter(p.adj < 0.05) %>%
  ggplot(aes(x = patient, y = cor, fill = patient)) +
  geom_boxplot() +
  labs(x = "Patient", y = "Spearman Correlation") +
  theme_classic() +
  theme(text = element_text(size = 18), legend.position = "none",
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
fig.page.bcs.corr

ggsave(plot = fig.page.bcs.corr,
       filename = paste0(out.dir, "PAGE_BCS_correlation.pdf"))
