# conda activate pub-semla
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(semla) # semla_1.1.6
set.seed(1)
out.dir <- "results/all/"

# --- Functions ---
# Performs edge neighbourhood enrichment
edge.neighbourhood <- function(patient.id, column = "main.cell.type",
                               n.perm = 1000) {
  # Read Seurat object
  seuratobj <- readRDS(seuratobjs.path[patient.id])
  # Add metadata
  metadata <- bcmerged@meta.data %>%
    filter(patient == patient.id) %>%
    select(major.tcs)
  seuratobj <- AddMetaData(seuratobj, metadata)
  # Semla object
  semlaobj <- UpdateSeuratForSemla(seuratobj,
                                   image_type = "tissue_lowres",
                                   verbose = FALSE)
  # Region neighbours
  semlaobj <- RegionNeighbors(semlaobj, column_name = "major.tcs",
                              column_labels = "Mixed",
                              mode = "inner_outer")
  # Create ID column
  ids <- semlaobj@meta.data %>%
    rownames_to_column("spots") %>%
    mutate(id = case_when(is.na(nb_to_Mixed) ~ NA_character_,
                          nb_to_Mixed == "nb_to_Mixed" ~ major.tcs,
                          TRUE ~ nb_to_Mixed)) %>%
    select(spots, id) %>%
    column_to_rownames("spots")
  semlaobj <- AddMetaData(semlaobj, ids)
  # Enrichment test
  res <- RunNeighborhoodEnrichmentTest(semlaobj, column_name = "id",
                                       n_permutations = 1000) %>%
    mutate(patient = patient.id)
  return(res)
}

# --- Data ---
# Merged sensitivity BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Paths to all Seurat objects
seuratobjs.path <-
  list.files("data/visium", pattern = ".rds", full.names = TRUE)
patients <- str_remove(basename(seuratobjs.path), pattern = "_.*$")
seuratobjs.path <- setNames(seuratobjs.path, patients)

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Perform edge neighbourhood enrichment analysis for each patient and major TC
enrichment.compartments <- lapply(patients, FUN = function(p) {
  edge.neighbourhood(p, n.perm = 1000)
}) |>
  bind_rows()

# Save
write.table(enrichment.compartments,
            file = paste0(out.dir, "neighbourhood_enrichment.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
