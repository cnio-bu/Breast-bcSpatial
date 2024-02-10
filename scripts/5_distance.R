# conda activate pub-semla
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(semla) # semla_1.1.6
library(rstatix)# rstatix_0.7.2
set.seed(1)
out.dir <- "results/all/"

# --- Functions ---
# Computes radial distances to the tumour core
radial.dist <- function(patient.id) {
  # Read Seurat object
  seuratobj <- readRDS(seuratobjs.path[patient.id])
  # Add metadata
  metadata <- bcsensitivity@meta.data %>%
    filter(patient == patient.id) %>%
    select(major.tcs) %>%
    mutate(major.tcs = if_else(major.tcs == "Tumour-rich", "tumour", major.tcs))
  seuratobj <- AddMetaData(seuratobj, metadata)
  # Semla object
  semlaobj <- UpdateSeuratForSemla(seuratobj, image_type = "tissue_lowres",
                                   verbose = FALSE)
  # Radial distances
  semlaobj <- RadialDistance(semlaobj, column_name = "major.tcs",
                             selected_groups = "tumour")
  # Transform distances to Seurat object
  distances <- semlaobj@meta.data %>%
    rownames_to_column("spots") %>%
    select(spots, subtype, slide, region, r_dist_tumour, major.tcs) %>%
    mutate(r_dist_tumour_squared = sign(r_dist_tumour) * sqrt(abs(r_dist_tumour)),
           major.tcs = if_else(major.tcs == "tumour", "Tumour-rich", major.tcs),
           patient = patient.id)
  return(distances)
}

# --- Data ---
# Merged sensitivity BCS with TCs
bcsensitivity <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Merged functional BCS
bcfunctional <- readRDS("results/all/functional/functional_bcs_regressed.rds")

# Paths to all Seurat objects
seuratobjs.path <-
  list.files("data/visium", pattern = ".rds", full.names = TRUE)
patients <- str_remove(basename(seuratobjs.path), pattern = "_.*$")
seuratobjs.path <- setNames(seuratobjs.path, patients)

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "sensitivity"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(paste0(out.dir, "functional"), recursive = TRUE,
           showWarnings = FALSE)

# Compute radial distances to the tumour core
radial.dist <- lapply(patients, FUN = function(p) {
  radial.dist(p)
}) |>
  bind_rows()

# Merge BCS, radial distances and metadata
sensitivity.normalised <- bcsensitivity@normalized %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("spots") %>%
  pivot_longer(cols = rownames(bcsensitivity@normalized),
               names_to = "signature", values_to = "enrichment") %>%
  left_join(radial.dist, by = "spots")

functional.normalised <- bcfunctional@normalized %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("spots") %>%
  pivot_longer(cols = rownames(bcfunctional@normalized),
               names_to = "signature", values_to = "enrichment") %>%
  left_join(radial.dist, by = "spots")

# Pearson correlation between distances and enrichment scores
sensitivity.corr <- sensitivity.normalised %>%
  group_by(patient, slide, region, signature) %>%
  cor_test(vars = enrichment, vars2 = r_dist_tumour, method = "pearson") %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup()

functional.corr <- functional.normalised %>%
  group_by(patient, slide, region, signature) %>%
  cor_test(vars = enrichment, vars2 = r_dist_tumour, method = "pearson") %>%
  adjust_pvalue(method = "fdr") %>%
  ungroup()

# Save
write.table(radial.dist, file = paste0(out.dir, "radial_distance.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(sensitivity.corr,
            file = paste0(out.dir, "sensitivity/sensitivity_correlation.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

write.table(functional.corr,
            file = paste0(out.dir, "functional/functional_correlation.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
