# conda activate pub-cellchat
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(CellChat) # CellChat_2.1.0
set.seed(1)
out.dir <- "results/tumour/CellChat/"

# --- Data ---
# Integrated Seurat object
seurat.integrated <- readRDS("results/all/integrated_seurat.rds")

# Tumour sensitivity BCS with TCs
bctumour <- readRDS("results/tumour/sensitivity/tcs_tumour.rds")

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Tumour TCs
tumour.tcs <- bctumour@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, tumour.tcs)

# Normalised counts
norm.counts <- GetAssayData(seurat.integrated, slot = "data", assay = "SCT")

# Metadata
metadata <- seurat.integrated@meta.data %>%
  rownames_to_column("spots") %>%
  left_join(tumour.tcs, by = "spots") %>%
  mutate(slices = factor(paste(patient, slide, sep = "_")),
         tumour.tcs = if_else(is.na(tumour.tcs), main.cell.type, tumour.tcs),
         tumour.tcs = factor(tumour.tcs)) %>%
  select(spots, slices, tumour.tcs) %>%
  rename(labels = tumour.tcs) %>%
  column_to_rownames("spots")

metadata <- metadata[colnames(norm.counts), ]

# Scale factors of spatial coordinates
# For 10X Visium, the conversion factor of converting spatial coordinates from
# Pixels to Micrometers can be computed as the ratio of the theoretical spot
# size (i.e. 65um) over the number of pixels that span the diameter of a
# theoretical spot size in the full-resolution image
# (i.e. 'spot_diameter_fullres' in pixels in the 'scalefactors_json.json'
# file)
# Of note, the 'spot_diameter_fullres' factor is different from the `spot`
# in Seurat object and thus users still need to get the value from the original
# json file.
spot.size <- 65 # the theoretical spot size (um) in 10X Visium
path <- list.files("data/visium/scalefactors", full.names = TRUE)
scalefactors <- lapply(path, FUN = function(x) {
  scale.factors <- jsonlite::fromJSON(txt = x)
  conversion.factor <- spot.size/scale.factors$spot_diameter_fullres
  scale.factors <- data.frame(ratio = conversion.factor, tol = spot.size/2)
  return(scale.factors)
}) |>
  bind_rows()
patients <- str_extract(basename(path), pattern = ".*_slide") %>%
  str_remove(pattern = "_slide")
slides <- str_extract(path, pattern = "_slide[0-9]_") %>%
  str_remove_all(pattern = "_")
rownames(scalefactors) <- paste(patients, slides, sep = "_")

# Coordinates
images <- names(seurat.integrated@images)
coordinates <- lapply(images, FUN = function(x) {
  GetTissueCoordinates(seurat.integrated, scale = NULL,
                       cols = c("imagerow", "imagecol"), image = x)
}) |>
  bind_rows()

# Create CellChat object
cellchat <- createCellChat(object = norm.counts, meta = metadata,
                           group.by = "labels", datatype = "spatial",
                           coordinates = coordinates,
                           scale.factors = scalefactors)

# Set ligand-receptor interaction database
CellChatDB <- CellChatDB.human

# Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <-
  subsetDB(CellChatDB, key = "annotation",
           search = c("Secreted Signaling", "ECM-Receptor",
                      "Cell-Cell Contact"))

# Set the used database in the object
cellchat@DB <- CellChatDB.use

# Pre-process the expression data
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = FALSE, interaction.range = 250,
                              contact.knn = TRUE, contact.knn.k = 6,
                              population.size = TRUE)

# Filter out the cell-cell communication if there are only few number of cells
# in certain groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Infer cell communication at the signalling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate the aggregated cell communication network
cellchat <- aggregateNet(cellchat)

# Network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat)

# Scan all possible signalling pathways
netAnalysis_signalingRole_heatmap(
  cellchat, signaling = cellchat@netP$pathways[1:64], pattern = "all",
  width = 20, height = 20)
netAnalysis_signalingRole_heatmap(
  cellchat, signaling = cellchat@netP$pathways[65:128], pattern = "all",
  width = 20, height = 20)

# Plot interesting signalling pathways
netVisual_heatmap(cellchat, signaling = "PD-L1", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "IGF", color.heatmap = "Reds")
netVisual_heatmap(cellchat, signaling = "DESMOSOME", color.heatmap = "Reds")

# Retrieve communication scores for interesting signalling pathways and cell
# types
tcs <- levels(tumour.tcs$tumour.tcs)
cell.types <- c("T cells", "CAFs", "B cells")
df.net <-
  subsetCommunication(cellchat, slot.name = "net",
                      sources.use = c(tcs, cell.types),
                      targets.use = c(tcs, cell.types),
                      signaling = c("PD-L1","IGF", "DESMOSOME"),
                      thresh = 0.05) %>%
  filter(source != target | (source %in% tcs & target %in% tcs))

# Scale the scores by interaction, between 0 and 1
df.net <- df.net %>%
  group_by(interaction_name_2) %>%
  mutate(prob = case_when(length(unique(prob)) > 1 ~
                            scales::rescale(prob, to = c(0, 1)),
                          TRUE ~ 1)) %>%
  ungroup()

# Save
write.table(df.net, file = paste0(out.dir, "signalling_scores_scaled.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
saveRDS(cellchat, file = paste0(out.dir, "spatial_cellchat.rds"))
