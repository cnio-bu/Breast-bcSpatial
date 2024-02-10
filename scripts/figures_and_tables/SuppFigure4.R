# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (1823-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
library(viridis) # viridis_0.6.4
library(rstatix) # rstatix_0.7.2
library(ggpubr) # ggpubr_0.6.0
library(patchwork) # patchwork_1.1.3
library(figpatch) # figpatch_0.2
set.seed(1)
out.dir <- "figures/"

# --- Data ---
# Tumour sensitivity BCS with TCs
bctumour <- readRDS("results/tumour/sensitivity/tcs_tumour.rds")

# Merged functional BCS
bcfunctional <- readRDS("results/all/functional/functional_bcs_regressed.rds")

# Radial distances to the tumour core
radial.dist <- read.table("results/all/radial_distance.tsv", header = TRUE,
                          sep = "\t")
# CellChat scores
cellchat.scores <-
  read.table("results/tumour/CellChat/signalling_scores_scaled.tsv",
             header = TRUE, sep = "\t")

# --- Code ---
# Create outdir
dir.create(paste0(out.dir, "SuppFigure4"), recursive = TRUE,
           showWarnings = FALSE)

# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")

# Get metadata
metadata <- bctumour@meta.data %>%
  rownames_to_column("spots") %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs),
         major.tcs = factor(major.tcs, levels = c("Tumour-rich", "Interphase",
                                                  "TME-rich")))

# Distribution of tumour TCs accross patients
figSupp4A <- metadata %>%
  count(patient, tumour.tcs) %>%
  ggplot(aes(x = patient, y = n, fill = tumour.tcs)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.tcs) +
  labs(x = "Patient ID", y = "Proportion", fill = "Tumour TC") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), axis.line = element_blank())
figSupp4A
ggsave(plot = figSupp4A,
       filename = paste0(out.dir, "SuppFigure4/patientTC.png"))

# Create tumour TC proportion table
proportion <- metadata %>%
  add_count(patient) %>%
  add_count(patient, tumour.tcs, name = "n.tumour.tcs") %>%
  mutate(prop = round(n.tumour.tcs/n, digits = 2)) %>%
  select(patient, subtype, tumour.tcs, prop) %>%
  unique() %>%
  pivot_wider(id_cols = c(patient, subtype), names_from = tumour.tcs,
              values_from = prop)
proportion[is.na(proportion)] <- 0
proportion <- proportion[match(patient.order, proportion$patient), ]
proportion <- proportion[c("patient", "subtype", "TC1.1", "TC1.2", "TC2",
                           "TC3")]
write.table(proportion, file = "tables/SuppTable9.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Distances to the tumour core
mean.radial.dist <- metadata %>%
  select(spots, tumour.tcs) %>%
  left_join(radial.dist, by = "spots") %>%
  mutate(ID = paste(patient, slide, region, sep = "_")) %>%
  filter(!ID %in% c("CID44971_slide1_A", "CID4535_slide_A")) %>%
  group_by(ID, subtype, tumour.tcs) %>%
  summarise(mean.dist = mean(r_dist_tumour),
            mean.dist.squared = sign(mean.dist) * sqrt(abs(mean.dist))) %>%
  ungroup()

figSupp4B <- mean.radial.dist %>%
  ggplot(aes(x = tumour.tcs, y = mean.dist.squared, fill = tumour.tcs)) +
  geom_boxplot() +
  scale_fill_manual(values = colour.tcs) +
  labs(x = "Tumour TC", y = "Distance to tumour core",
       fill = "Tumour TC") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
figSupp4B
ggsave(plot = figSupp4B, width = 10,
       filename = paste0(out.dir, "SuppFigure4/distanceTC.png"))

# Number of tumour spots in each region, per TC
metadata %>%
  add_count(tumour.tcs) %>%
  filter(major.tcs == "Tumour-rich") %>%
  add_count(tumour.tcs, major.tcs, name = "n.tumour.rich") %>%
  mutate(prop.tumour.rich = round(n.tumour.rich/n, digits = 2)) %>%
  select(tumour.tcs, prop.tumour.rich) %>%
  unique() %>%
  arrange(tumour.tcs)

# Number of interphase spots in each region, per TC
metadata %>%
  add_count(tumour.tcs) %>%
  filter(major.tcs == "Interphase") %>%
  add_count(tumour.tcs, major.tcs, name = "n.interphase") %>%
  mutate(prop.interphase = round(n.interphase/n, digits = 2)) %>%
  select(tumour.tcs, prop.interphase) %>%
  unique() %>%
  arrange(tumour.tcs)

figSupp4C <- metadata %>%
  count(tumour.tcs, major.tcs) %>%
  ggplot(aes(x = tumour.tcs, y = n, fill = major.tcs)) +
  geom_bar(position = position_fill(reverse = TRUE),
           stat = "identity") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values = colour.tctypes) +
  labs(x = "Major TC", y = "Proportion", fill = "Tumour TC") +
  theme_classic() +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank(), axis.line = element_blank())
figSupp4C
ggsave(plot = figSupp4C,
       filename = paste0(out.dir, "SuppFigure4/compartmentTC.png"))

# CellChat communication scores
group.order <- c("TC1.1", "TC1.2", "TC2", "T cells", "CAFs")
figSupp4D <- cellchat.scores %>%
  mutate(source.target = paste(source, target, sep = " -> "),
         source = factor(source, levels = group.order),
         target = factor(target, levels = group.order)) %>%
  arrange(as.numeric(source), as.numeric(target)) %>%
  mutate(source.target = factor(source.target, levels = unique(source.target))) %>%
  ggplot(aes(x = source.target, y = interaction_name_2, color = prob)) +
  geom_point(aes(size = prob), shape = 16) +
  scale_color_viridis(breaks = c(0, 1),labels = c("min","max")) +
  scale_size_continuous(range = c(7, 13)) +
  labs(colour = "Interaction\nstrength") +
  guides(size = "none") +
  theme_linedraw() +
  theme(legend.position = "right",
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        text = element_text(size = 18))
figSupp4D
ggsave(plot = figSupp4D, width = 10,
       filename = paste0(out.dir, "SuppFigure4/cellchat.png"))

# TIS and ECM stiffness enrichment
tumour.tcs <- bctumour@meta.data %>%
  rownames_to_column("spots") %>%
  select(spots, patient, tumour.tcs)

functional <- t(bcfunctional@normalized[c("TIS", "ECM_STIFFNESS"), ]) %>%
  as.data.frame() %>%
  rownames_to_column("spots") %>%
  left_join(tumour.tcs, by = "spots") %>%
  filter(!is.na(tumour.tcs)) %>%
  mutate(tumour.tme = factor(tumour.tcs))

# Differential enrichment of TIS
functional.tis <- functional %>%
  group_by(patient) %>%
  filter("TC1.1" %in% tumour.tcs) %>%
  ungroup()

functional.tis %>%
  group_by(patient) %>%
  kruskal_test(TIS ~ tumour.tcs) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

stat.test <- functional.tis %>%
  group_by(patient) %>%
  pairwise_wilcox_test(TIS ~ tumour.tcs) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  add_y_position()

figSupp4E <- functional.tis %>%
  ggboxplot(x = "tumour.tcs", y = "TIS", fill = "tumour.tcs",
            facet.by = "patient") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  facet_grid(~patient, scales = "free_x") +
  scale_fill_manual(values = colour.tcs) +
  labs(x = NULL, y = "ICI sensitivity", fill = "Tumour TC") +
  pub_theme_facet +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
figSupp4E
ggsave(plot = figSupp4E, width = 15, height = 5,
       filename = paste0(out.dir, "SuppFigure4/TIS_enrichment.png"))

# Differential enrichment of ECM stiffness
functional.stiffness <- functional %>%
  group_by(patient) %>%
  filter("TC1.2" %in% tumour.tcs) %>%
  ungroup()

functional.stiffness %>%
  group_by(patient) %>%
  kruskal_test(ECM_STIFFNESS ~ tumour.tcs) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

stat.test <- functional.stiffness %>%
  group_by(patient) %>%
  pairwise_wilcox_test(ECM_STIFFNESS ~ tumour.tcs) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance() %>%
  add_y_position()

figSupp4F <- functional.stiffness %>%
  ggboxplot(x = "tumour.tcs", y = "ECM_STIFFNESS", fill = "tumour.tcs",
            facet.by = "patient") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  facet_grid(~patient, scales = "free_x") +
  scale_fill_manual(values = colour.tcs) +
  labs(x = NULL, y = "ECM stifness enrichment", fill = "Tumour TC") +
  pub_theme_facet +
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1))
figSupp4F
ggsave(plot = figSupp4F, width = 15, height = 5,
       filename = paste0(out.dir, "SuppFigure4/ECMstiffness_enrichment.png"))

# Save figure
top <- (figSupp4A | figSupp4B) +
  plot_layout(widths = c(1.5, 1))
second <- (figSupp4C | figSupp4D) +
  plot_layout(widths = c(1.5, 2.5))
third <- (plot_spacer() | figSupp4E | plot_spacer()) +
  plot_layout(widths = c(0.05, 1, 0.05))
bottom <- (plot_spacer() | figSupp4F | plot_spacer()) +
  plot_layout(widths = c(0.18, 1, 0.18))

figSupp4 <- (top / second / third / bottom) +
  plot_layout(heights = c(1.5, 1.5, 1.5, 1.5)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))

ggsave(plot = figSupp4, width = 21, height = 24,
       filename = paste0(out.dir, "SuppFigure4/SuppFigure4.pdf"))
