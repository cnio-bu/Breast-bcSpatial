# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(eulerr) # eulerr_7.0.0
library(viridis) # viridis_0.6.4
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
set.seed(1)
out.dir <- "figures/reviewers/benchmark/deconvolution/"

# --- Functions ---
# Reads deconvolution results
readDeconv <- function(type = c("giotto", "card")) {
  # Read data
  deconvolution <- lapply(patients, FUN = function(p) {
    read.table(paste0("data/benchmarking/deconvolution/", type, "/", p,
                      "_deconvolution.tsv"), sep = "\t", header = TRUE)
  }) |>
    bind_rows() %>%
    rownames_to_column("spots")
  return(deconvolution)
}

# Classifies spots into tumour or TME.
# For TME spots, asigns the cell type with maximum deconvoluted proportion
tumourTME <- function(deconvolution, other.sources = annotation) {
  deconvolution %>%
    left_join(other.sources, by = "spots") %>%
    mutate(tumour.patho = Pathologist == "DCIS" |
             str_detect(tolower(Pathologist), pattern = "cancer") |
             str_detect(tolower(Pathologist), pattern = "carcinoma"),
           tumour.estimate = ESTIMATEScaled <= 0.4,
           tumour.deconv = Cancer.Epithelial >= 0.6,
           tumour.scevan = !is.na(subclone),
           tumour = case_when(tumour.scevan ~ "Tumour",
                              tumour.patho & tumour.estimate ~ "Tumour",
                              tumour.patho & tumour.deconv ~ "Tumour",
                              tumour.estimate & tumour.deconv ~ "Tumour",
                              TRUE ~ "TME")) %>%
    pivot_longer(cols = colnames(deconvolution)[-1], names_to = "cell.type",
                 values_to = "cell.prop") %>%
    filter(cell.type != "Cancer.Epithelial") %>%
    group_by(spots) %>%
    filter(cell.prop == max(cell.prop)) %>%
    ungroup() %>%
    add_count(spots) %>%
    mutate(main.cell.type = case_when(n > 1 ~ "Undetermined",
                                      tumour == "Tumour" ~ "Cancer.Epithelial",
                                      TRUE ~ cell.type),
           main.cell.type = str_replace_all(main.cell.type, pattern = "\\.",
                                            replacement = " ")) %>%
    select(spots, tumour, main.cell.type, tumour.patho:tumour.scevan) %>%
    unique()
}

# Computes overlap between annotation sources
getOverlap <- function(deconvolution) {
  list(Pathologist = deconvolution %>% filter(tumour.patho) %>% pull(spots),
       ESTIMATE = deconvolution %>% filter(tumour.estimate) %>% pull(spots),
       Deconvolution = deconvolution %>% filter(tumour.deconv) %>% pull(spots),
       SCEVAN = deconvolution %>% filter(tumour.scevan) %>% pull(spots))
}

# Computes cosine similarity
cosine.similarity <- function(v1, v2) {
  sum(v1 * v2) / (sqrt(sum(v1^2)) * sqrt(sum(v2^2)))
}

# --- Data ---
# Merged sensitivity BCS with TCs
bcsensitivity <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# Patient IDs
patients <- unique(bcsensitivity@meta.data$patient)

# Other deconvolution results
dwls <- readDeconv("giotto")
card <- readDeconv("card")

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")
colour.venn <-
  c("Pathologist" = "#CEBCDB", "ESTIMATE" = "darkgoldenrod1",
    "Deconvolution" = "#A3BEE3", "SCEVAN", "#AEF067")

# Extract annotation sources (pathologist, ESTIMATE and SCEVAN)
annotation <- bcsensitivity@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, Pathologist, ESTIMATEScaled, subclone)

# Classify spots into Tumour or TME (with main cell type)
dwls <- tumourTME(dwls) %>%
  mutate(deconvolution = "spatialDWLS")
card <- tumourTME(card) %>%
  mutate(deconvolution = "CARD")

# Extract major TCs and RCTD deconvolution results
metadata <- bcsensitivity@meta.data %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs =
           factor(major.tcs,
                  levels = c("Tumour-rich", "Interphase", "TME-rich")),
         patient = factor(patient, levels = patient.order)) %>%
  select(spots, patient, major.tcs)

rctd <- bcsensitivity@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, B.cells:T.cells) |>
  tumourTME() %>%
  mutate(deconvolution = "RCTD")

# Figure 1B
deconvolution <- list(CARD = card, RCTD = rctd, spatialDWLS = dwls)
fig1B <- lapply(seq_along(deconvolution), FUN = function(i) {
  plot(euler(getOverlap(deconvolution[[i]])), labels = NULL,
       quantities = list(fontsize = 15), main = names(deconvolution)[i],
       fills = c("#CEBCDB", "darkgoldenrod1", "#A3BEE3", "#AEF067"))
}) |>
  wrap_plots(ncol = 3)
fig1B
ggsave(plot = fig1B, width = 15, height = 6,
       filename = paste0(out.dir, "venn.pdf"))

# Merge deconvolution results and metadata
cols <- c("spots", "main.cell.type", "deconvolution")
metadata <- rbind(rctd[, cols], dwls[, cols], card[, cols]) %>%
  left_join(metadata, by = "spots")

# Figure 3A
id.levels <- paste(rep(sort(unique(metadata$deconvolution)),
                       each = length(patient.order)), patient.order, sep = "\n")
fig3A <- metadata %>%
  mutate(ID = paste(deconvolution, patient, sep = "\n"),
         ID = factor(ID, levels = id.levels)) %>%
  ggplot(aes(x = major.tcs, fill = main.cell.type)) +
  geom_bar(position = "fill") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
  facet_wrap(~ID, nrow = 3) +
  scale_fill_manual(values = colour.celltypes) +
  labs(y = "Proportion", fill = "Cell type") +
  pub_theme_facet +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 18))
fig3A
ggsave(plot = fig3A, height = 10, width = 20,
       filename = paste0(out.dir, "cell_proportions.pdf"))

# Summarize deconvolution results
metadata.grouped <- metadata %>%
  count(patient, major.tcs, deconvolution, main.cell.type)

# Pivot to wide format
metadata.wide <- metadata.grouped %>%
  pivot_wider(names_from = deconvolution, values_from = n, values_fill = 0)

# Cosine similarity
metadata.cosine <- metadata.wide %>%
  group_by(patient, major.tcs) %>%
  mutate(ID = paste(patient, major.tcs, sep = "_"),
         similarity.RCTD.vs.CARD = cosine.similarity(RCTD, CARD),
         similarity.RCTD.vs.DWLS = cosine.similarity(RCTD, spatialDWLS),
         similarity.CARD.vs.DWLS = cosine.similarity(CARD, spatialDWLS)) %>%
  ungroup()

metadata.cosine %>%
  summarise(mean.RCTD.CARD = mean(similarity.RCTD.vs.CARD),
            mean.RCTD.DWLS = mean(similarity.RCTD.vs.DWLS),
            mean.CARD.DWLS = mean(similarity.CARD.vs.DWLS))

# Heatmap
heat <- metadata.cosine %>%
  pivot_longer(cols = c("similarity.RCTD.vs.CARD", "similarity.RCTD.vs.DWLS",
                        "similarity.CARD.vs.DWLS"), names_to = "method",
               names_prefix = "similarity.", values_to = "similarity") %>%
  mutate(method = str_replace_all(method, pattern = "\\.", replacement = " ")) %>%
  ggplot(aes(method, ID, fill = similarity)) +
  geom_tile() +
  scale_fill_viridis() +
  labs(x = NULL, y = NULL, fill = "Cosine similarity") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 18))
heat
ggsave(plot = heat, width = 10, height = 10,
       filename = paste0(out.dir, "cosine_similarity.pdf"))
