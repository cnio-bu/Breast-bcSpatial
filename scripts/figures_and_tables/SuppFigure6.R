# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(Seurat) # Seurat_4.3.0.1
library(SeuratObject) # SeuratObject_4.1.3
library(beyondcell) # beyondcell_2.2.0
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/"

# --- Data ---
# Beyondcell object
bcobj.tnbc <-
  readRDS("data/beyondcell/sensitivity/CID44971_spatial_bcs_sensitivity.rds")

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "SuppFigure6"), recursive = TRUE,
           showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Plot TNBC subclones
figSupp6A <- bcClusters(bcobj.tnbc, idents = "subclone", spatial = TRUE,
                           image.alpha = 0.7, pt.size = 1.7) +
  scale_fill_manual(values = colour.subclones, na.value = NA) +
  ggtitle(NULL) +
  labs(fill = "Subclone") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(text = element_text(size = 18), legend.key = element_rect(fill = NA))
figSupp6A

subclones <- bcobj.tnbc@meta.data %>%
  select(region, subclone) %>%
  rownames_to_column("spots")

# Plot docetaxel response
response.docetaxel <-
  bcSignatures(bcobj.tnbc, signatures = list(values = "docetaxel_CTRP_A05821830"),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) +
  ggtitle(NULL) +
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5),
        legend.position = "bottom")

# Plot response by subclone
sensitivity.full <- bcobj.tnbc@normalized %>%
  t() %>%
  as.data.frame() %>%
  select(docetaxel_CTRP_A05821830) %>%
  rownames_to_column("spots") %>%
  left_join(subclones, by = "spots") %>%
  filter(region != "A" & subclone %in% c("1", "2", "3", "4")) %>%
  add_count(region, subclone) %>%
  filter(n >= 10)

boxplot.docetaxel <- sensitivity.full %>%
  ggplot(aes(x = subclone, y = docetaxel_CTRP_A05821830, fill = subclone)) +
  geom_boxplot() +
  scale_fill_manual(values = colour.subclones) +
  facet_wrap(~region) +
  labs(x = "Subclone", y = "Sensitivity", fill = "Subclone") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))

figSupp6B <- (response.docetaxel | boxplot.docetaxel) +
  plot_annotation(title = "Docetaxel (CTRP) sensitivity",
                  subtitle = "Microtubule agent, Taxane") &
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5))
figSupp6B

# Subset object
bcobj.tnbc.subset <- bcSubset(bcobj.tnbc, cells = sensitivity.full$spots)

# Identify specific drugs for region C
bcobj.tnbc.subset <- bcRanks(bcobj.tnbc.subset, idents = "subclone")
diff.response <- "TOP-Differential-HighSensitivity"
response <- bcobj.tnbc.subset@ranks$subclone %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("group"), names_to = "subclone",
               values_to = "response", names_prefix = "group.") %>%
  filter(response %in% diff.response) %>%
  select(IDs, subclone, response) %>%
  filter(subclone == "4")

# Plot drug response
response.oxaliplatin <-
  bcSignatures(bcobj.tnbc, signatures = list(values = "Oxaliplatin_GDSC_1806"),
               spatial = TRUE, pt.size = 1.5, image.alpha = 0.7) +
  ggtitle(NULL) +
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5),
        legend.position = "bottom")

boxplot.oxaliplatin <- bcobj.tnbc.subset@normalized %>%
  t() %>%
  as.data.frame() %>%
  select(Oxaliplatin_GDSC_1806) %>%
  rownames_to_column("spots") %>%
  left_join(subclones, by = "spots") %>%
  pivot_longer(cols = "Oxaliplatin_GDSC_1806", names_to = "IDs", values_to = "BCS") %>%
  ggplot(aes(x = subclone, y = BCS, fill = subclone)) +
  geom_boxplot() +
  scale_fill_manual(values = colour.subclones) +
  facet_wrap(~region ) +
  labs(x = "Subclone", y = "Sensitivity", fill = "Subclone") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))

figSupp6C <- (response.oxaliplatin | boxplot.oxaliplatin) +
  plot_annotation(title = "Oxaliplatin (GDSC) sensitivity",
                  subtitle = "DNA replication inhibitor, Platinum agent") &
  theme(text = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 22, hjust = 0.5))
figSupp6C

# Save figure
figSupp6 <- ((figSupp6A | figSupp6B) + plot_layout(widths = c(1, 3))) /
  ((plot_spacer() | figSupp6C) + plot_layout(widths = c(1, 3))) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = figSupp6, width = 21, height = 16,
       filename = paste0(out.dir, "SuppFigure6/SuppFigure6.pdf"))
