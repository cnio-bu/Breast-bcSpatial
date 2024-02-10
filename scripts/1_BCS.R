# conda activate pub-beyondcell
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
set.seed(1)
out.dir <- "results/all/"

# --- Functions ---
# Merges Beyondcell objects computed from different samples
bcJoin <- function(bc1, bc2) {
  genes <- unique(c(rownames(bc1@expr.matrix), rownames(bc2@expr.matrix)))
  spots <- unique(c(colnames(bc1@expr.matrix), colnames(bc2@expr.matrix)))
  sigs <- unique(c(rownames(bc1@normalized), rownames(bc2@normalized)))
  # Create new object
  bc <- bc1
  # Fill missing signatures with 0s
  in.bc1 <- sigs %in% rownames(bc1@data)
  if (!all(in.bc1)) {
    zeros.bc1 <- matrix(rep(0, times = ncol(bc1@data) * sum(!in.bc1)),
                        ncol = ncol(bc1@data),
                        dimnames = list(sigs[!in.bc1], colnames(bc1@data)))
    bc1@data <- rbind(bc1@data, zeros.bc1)
    bc1 <- bcRecompute(bc1, slot = "data")
  }
  in.bc2 <- sigs %in% rownames(bc2@data)
  if (!all(in.bc2)) {
    zeros.bc2 <- matrix(rep(0, times = ncol(bc2@data) * sum(!in.bc2)),
                        ncol = ncol(bc2@data),
                        dimnames = list(sigs[!in.bc2], colnames(bc2@data)))
    bc2@data <- rbind(bc2@data, zeros.bc2)
    bc2 <- bcRecompute(bc2, slot = "data")
  }
  # Update matrices
  bc@data <- cbind(bc1@data[sigs, , drop = FALSE],
                   bc2@data[sigs, , drop = FALSE])
  bc@expr.matrix <- matrix(nrow = length(genes), ncol = length(spots),
                           dimnames = list(genes, spots))
  bc@expr.matrix[rownames(bc1@expr.matrix), colnames(bc1@expr.matrix)] <- bc1@expr.matrix
  bc@expr.matrix[rownames(bc2@expr.matrix), colnames(bc2@expr.matrix)] <- bc2@expr.matrix
  bc@expr.matrix[is.na(bc@expr.matrix)]<- 0
  bc@background <- matrix(nrow = 0, ncol = 0)
  # Update metadata
  bc@meta.data <- rbind(bc1@meta.data, bc2@meta.data)
  # Update parameters
  bc@n.genes <- unique(bc1@n.genes, bc2@n.genes)
  bc@mode <- unique(c(bc1@mode, bc2@mode))
  bc@thres <- unique(bc1@thres, bc2@thres)
  # Update rest of slots
  bc@ranks <- list()
  bc@reductions <- list()
  bc@regression <- list(order = c("", ""), vars = NULL,
                        order.background = c("", ""))
  bc@SeuratInfo <- c(bc1@SeuratInfo, bc2@SeuratInfo)
  # Recompute
  bc <- bcRecompute(bc, slot = "data")
  return(bc)
}

# --- Data ---
# Patient IDs
bcobjs.path <- list.files("data/beyondcell/sensitivity", full.names = TRUE)
patients <- str_remove(basename(bcobjs.path), pattern = "_.*$")

# --- Code ---
# Create outdirs
dir.create(paste0(out.dir, "sensitivity"), recursive = TRUE,
           showWarnings = FALSE)
dir.create(paste0(out.dir, "functional"), recursive = TRUE,
           showWarnings = FALSE)

# Sensitivity and functional BCS
bcs.file.suffix <- list(sensitivity = "_spatial_bcs_sensitivity.rds",
                        functional = "_spatial_bcs_functional.rds")

# For each BCS...
for (bc in names(bcs.file.suffix)) {
  # Create a list of Beyondcell objects
  bcs <- lapply(patients, FUN = function(p) {
    bcobj <- readRDS(paste0("data/beyondcell/", bc, "/", p,
                            bcs.file.suffix[[bc]]))
    bcobj@meta.data <- bcobj@meta.data %>%
      select(patient, subtype, slide, region, x:z, nCount_SCT, nFeature_SCT,
             Phase, Pathologist, Tumour, ESTIMATEScaled, subclone,
             B.cells:T.cells)
    return(bcobj)
  })
  # Merge Beyondcell objects
  bcmerged <- reduce(bcs, bcJoin)
  # Type of sample
  bcmerged@meta.data <- bcmerged@meta.data %>%
    mutate(sample.type = if_else(patient == "738811QB", "FFPE", "Fresh Frozen"))
  # Replace NAs by 0s
  bcmerged@normalized[is.na(bcmerged@normalized)] <- 0
  bcmerged <- bcRecompute(bcmerged, slot = "normalized")
  # Regress patient effect
  bcregressed <- bcRegressOut(bcmerged, vars.to.regress = "patient")
  # Save
  if (bc == "sensitivity") {
    saveRDS(bcmerged, file = paste0(out.dir, bc,  "/", bc, "_bcs.rds"))
  }
  saveRDS(bcregressed,
          file = paste0(out.dir, bc, "/", bc, "_bcs_regressed.rds"))
}
