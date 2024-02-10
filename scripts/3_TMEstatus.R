# conda activate pub-singscore
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(singscore) # singscore_1.22.0
set.seed(1)
out.dir <- "results/all/tme/"

# --- Data ---
# Patient IDs
bcobjs.path <- list.files("data/beyondcell/sensitivity", full.names = TRUE)
patients <- str_remove(basename(bcobjs.path), pattern = "_.*$")

# Bagaev et al.
protumour <- c("PROTUMOR_CYTOKINES", "TUMOR-ASSOCIATED_MACROPHAGES",
               "TH2_SIGNATURE", "NEUTROPHIL_SIGNATURE", "CHECKPOINT_MOLECULES",
               "TREG", "GRANULOCYTE_TRAFFIC", "TREG_AND_TH2_TRAFFIC",
               "MYELOID_CELLS_TRAFFIC", "IMMUNE_SUPPRESSION_BY_MYELOID_CELLS",
               "MACROPHAGE_AND_DC_TRAFFIC")
antitumour <- c("MHCI", "MHCII", "ANTITUMOR_CYTOKINES", "M1_SIGNATURE",
                "T_CELLS", "NK_CELLS", "CO-ACTIVATION_MOLECULES", "B_CELLS",
                "TH1_SIGNATURE", "EFFECTOR_CELL_TRAFFIC", "EFFECTOR_CELLS")

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# TME genesets
tme.subtypes <- GenerateGenesets("data/signatures/TME_subtypes.gmt")

# For each patient...
for (p in patients) {
  seuratobj <-
    readRDS(paste0("data/visium/", p, "_spatial_SCTnormalised_clones.rds"))
  # Compute scores
  bcobj.tme <- bcScore(seuratobj, tme.subtypes)
  # Replace NAs by 0s
  bcobj.tme@normalized[is.na(bcobj.tme@normalized)] <- 0
  bcobj.tme <- bcRecompute(bcobj.tme, slot = "normalized")
  # Apply singscore to compute a TME status score
  rankdata <- rankGenes(bcobj.tme@normalized)
  scoredf <- simpleScore(rankdata, upSet = protumour, downSet = antitumour)
  TME.score <- scoredf %>%
    select(TotalScore) %>%
    rename(TME.score = TotalScore)
  # Add metadata
  bcobj.tme <- bcAddMetadata(bcobj.tme, TME.score)
  # Save
  saveRDS(bcobj.tme, file = paste0(out.dir, p, "_spatial_tme_status.rds"))
}
