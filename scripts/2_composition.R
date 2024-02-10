# conda activate pub-beyondcell
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "results/all/sensitivity/"

# --- Functions ---
# Reads all deconvolution files in a folder
read.all.deconv <- function(pathdir) {
  deconv.files <- list.files(pathdir, pattern = "deconvolution.tsv",
                             full.names = TRUE)
  lapply(deconv.files, FUN = function(x) {
    patient <- str_remove(basename(x), pattern = "_.*$")
    read.table(x, header = TRUE, sep = "\t") %>%
      rownames_to_column("spots") %>%
      mutate(spots = paste(spots, patient, sep = "_"))
  }) %>%
    bind_rows()
}

# --- Data ---
# Merged sensitivity BCS without regression of the patient effect
bcnoregress <- readRDS("results/all/sensitivity/sensitivity_bcs.rds")

# Merged sensitivity BCS with the patient effect regressed out
bcmerged <- readRDS("results/all/sensitivity/sensitivity_bcs_regressed.rds")

# --- Code ---
# Create outdirs
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
dir.create("figures/UMAP/", recursive = TRUE, showWarnings = FALSE)
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# Patient order
patient.order <- c("CID4290", "CID4535", "1142243F", "1160920F", "CID4465",
                   "CID44971", "738811QB", "1168993F", "V19L29")

# Number of signatures and unique drugs
nrow(bcmerged@normalized)
drugs <- rownames(bcmerged@normalized) |>
  str_remove(pattern = "_.*$") |>
  toupper()
length(unique(drugs))

# Compute TCs at res = 0.5
bcnoregress <- bcUMAP(bcnoregress, pc = 20, k.neighbors = 20, npcs = 50,
                      res = 0.5)
bcnoregress@meta.data <- bcnoregress@meta.data %>%
  rename(tcs = bc_clusters_res.0.5) %>%
  mutate(tcs = as.numeric(as.character(tcs)),
         tcs = factor(paste0("TC", tcs + 1), levels = paste0("TC", 1:18)))

bcmerged <- bcUMAP(bcmerged, pc = 20, k.neighbors = 20, npcs = 50, res = 0.5)

# Define major TCs
metadata <- bcmerged@meta.data %>%
  rownames_to_column("spots") %>%
  rename(tcs = bc_clusters_res.0.5)

prop <- metadata %>%
  add_count(tcs) %>%
  add_count(tcs, Tumour, name = "n.tumour") %>%
  mutate(prop = round(n.tumour/n, digits = 2)) %>%
  select(tcs, Tumour, prop) %>%
  unique() %>%
  mutate(zero.prop = Tumour == "Non-tumour" & prop == 1,
         Tumour = if_else(zero.prop, "Tumour", Tumour),
         prop = if_else(zero.prop, 0, prop)) %>%
  select(-zero.prop)
prop %>%
  filter(Tumour == "Tumour") %>%
  arrange(tcs)

metadata <- prop %>%
  filter(Tumour == "Tumour") %>%
  select(-Tumour) %>%
  left_join(metadata, by = "tcs", relationship = "many-to-many") %>%
  mutate(tcs = as.numeric(as.character(tcs)),
         tcs = factor(paste0("TC", tcs + 1), levels = paste0("TC", 1:16)),
         major.tcs = case_when(prop <= 0.15 ~ "TME-rich",
                               prop >= 0.95 ~ "Tumour-rich",
                               TRUE ~ "Mixed"))

proportion <- metadata %>%
  select(major.tcs, tcs, prop) %>%
  unique() %>%
  select(tcs, major.tcs, prop) %>%
  arrange(desc(prop))
proportion

write.table(proportion, file = "tables/SuppTable2.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

bcmerged@meta.data <- metadata %>%
  select(-prop) %>%
  column_to_rownames("spots") %>%
  unique()

# UMAP before and after regression and by technology
before <- (bcClusters(bcnoregress, idents = "tcs") + ggtitle("No regressed") +
             theme_classic()) /
  (bcClusters(bcnoregress, idents = "patient") +
     scale_colour_manual(values = scales::hue_pal()(length(patient.order)),
                         breaks = patient.order) +
     theme_classic())

after <- (bcClusters(bcmerged, idents = "tcs") + ggtitle("Regressed") +
            theme_classic()) /
  (bcClusters(bcmerged, idents = "patient") +
     scale_colour_manual(values = scales::hue_pal()(length(patient.order)),
                         breaks = patient.order) +
     theme_classic())

sampletype <- ((bcClusters(bcmerged, idents = "sample.type") +
                  theme_classic()) |
                 (bcClusters(bcmerged, idents = "sample.type",
                             split.by = "sample.type") + theme_classic())) +
  plot_layout(widths = c(1, 2))

umap.regression <- before | after
umap.regression
sampletype

ggsave(plot = umap.regression, width = 18, height = 15,
       filename = paste0("figures/UMAP/UMAP_regression.png"))
ggsave(plot = sampletype, filename = paste0("figures/UMAP/UMAP_sampletype.png"),
       width = 28)

# Assign each spot to the cell type with the maximum proportion
spot.type <- metadata %>%
  select(-prop) %>%
  pivot_longer(cols = B.cells:T.cells, names_to = "cell.type",
               values_to = "cell.prop") %>%
  filter(cell.type != "Cancer.Epithelial") %>%
  group_by(spots) %>%
  filter(cell.prop == max(cell.prop)) %>%
  ungroup() %>%
  mutate(main.cell.type = case_when(Tumour == "Tumour" ~ "Cancer.Epithelial",
                                    TRUE ~ cell.type))

# Format cell types
spot.type <- spot.type %>%
  mutate(main.cell.type = str_replace(main.cell.type, pattern = "\\.",
                                      replacement = " "))
# Update metadata
new.metadata <- spot.type %>%
  select(spots, main.cell.type) %>%
  unique() %>%
  column_to_rownames("spots")

bcmerged <- bcAddMetadata(bcmerged, new.metadata)

# Save
saveRDS(bcmerged, file = paste0(out.dir, "tcs_bcs.rds"))
