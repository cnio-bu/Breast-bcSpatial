# conda activate pub-cellchat
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
library(CellChat) # CellChat_2.1.0
library(ggsankey) # ggsankey_0.0.99999
library(rstatix) # rstatix_0.7.2
set.seed(1)
out.dir <- "figures/reviewers/"

# --- Data ---
# Tumour sensitivity BCS with TCs
bctumour <- readRDS("results/tumour/sensitivity/tcs_tumour.rds")

# CellChat output
cellchat <- readRDS("results/tumour/CellChat/spatial_cellchat.rds")

# --- Code ---
# Aesthetics
load("scripts/figures_and_tables/aesthetics.RData")
patient.ids <- unique(bctumour@meta.data$patient)
subclone.n <- bctumour@meta.data %>%
  rownames_to_column("spots") %>%
  select(patient, subclone) %>%
  filter(subclone != "diploid") %>%
  unique() %>%
  count(patient) %>%
  column_to_rownames("patient")

colour.subclones <- c(colour.subclones, "#C6B4DB", "#7974D3")
colour.subclones <- rep(colour.subclones, times = subclone.n[patient.ids, ])
names.colour <- sapply(patient.ids, FUN = function(x) 1:subclone.n[x, "n"]) |>
  unlist()
names.colour <- paste0(rep(patient.ids, times = subclone.n[patient.ids, ]),
                       " (SC", names.colour, ")")
colour.subclones <- setNames(colour.subclones, names.colour)
colour.all <- c(colour.subclones, colour.tcs, colour.celltypes, colour.annotations,
                c("CAFs IGF+" = "#F3D000"))

# Create outdirs
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Tumour TCs
tumour.tcs <- c("TC1.1", "TC1.2", "TC2", "TC3")

# Interaction activity for all cell types
weights <- cellchat@net[["weight"]] %>%
  as.data.frame() %>%
  rownames_to_column("sender") %>%
  pivot_longer(cols = colnames(cellchat@net[["weight"]]), names_to = "receiver",
               values_to = "weight") %>%
  group_by(sender, receiver) %>%
  mutate(interaction = paste(sort(c(sender, receiver)), collapse = "-")) %>%
  group_by(interaction) %>%
  summarise(value = sum(weight)) %>%
  ungroup() %>%
  mutate(TME = str_remove(interaction, pattern = "-.*"),
         TME = if_else(TME %in% tumour.tcs, "Tumour", TME),
         tumour = str_remove(interaction, pattern = "^.*-")) %>%
  group_by(TME, tumour) %>%
  summarise(value = sum(value)) %>%
  filter(tumour %in% tumour.tcs)

# Interaction activity for IGF+ CAFs
weights.igf <- cellchat@netP$prob[, , "IGF"] %>%
  as.data.frame() %>%
  rownames_to_column("sender") %>%
  pivot_longer(cols = colnames(cellchat@net[["weight"]]), names_to = "receiver",
               values_to = "weight") %>%
  group_by(sender, receiver) %>%
  mutate(interaction = paste(sort(c(sender, receiver)), collapse = "-")) %>%
  group_by(interaction) %>%
  summarise(value = sum(weight)) %>%
  ungroup() %>%
  mutate(TME = str_remove(interaction, pattern = "-.*"),
         tumour = str_remove(interaction, pattern = "^.*-")) %>%
  filter(TME == "CAFs" & tumour %in% tumour.tcs)

# Plot
weights %>%
  ggplot(aes(x = TME, y = value, fill = tumour)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = colour.tcs)

weights %>%
  filter(TME != "CAFs") %>%
  bind_rows(weights.igf) %>%
  mutate(TME = if_else(TME == "CAFs", "CAFs IGF+", TME)) %>%
  ggplot(aes(x = TME, y = value, fill = tumour)) +
  geom_col(position = "fill") +
  scale_fill_manual(values = colour.tcs)

# Sankey diagram
weights.full <- weights %>%
  rename(tumour.tcs = tumour) %>%
  group_by(TME) %>%
  mutate(sum.value = sum(value)) %>%
  ungroup() %>%
  filter(sum.value > 0) %>%
  mutate(percentage = round((value/sum.value) * 100), digits = 0)

ocurrences <- bctumour@meta.data %>%
  select(patient, subclone, tumour.tcs) %>%
  filter(subclone != "diploid") %>%
  mutate(subclone = paste0(patient, " (SC", subclone, ")")) %>%
  select(subclone, tumour.tcs) %>%
  left_join(weights.full, by = c("tumour.tcs"), relationship = "many-to-many")

ocurrences <- lapply(ocurrences, rep, ocurrences$percentage) |>
  bind_rows()

sankey.caf <- ocurrences %>%
  make_long(subclone, tumour.tcs, TME) %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node,
             fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE) +
  geom_sankey_label(size = 3, color = "black") +
  scale_fill_manual(values = colour.all) +
  theme_void() +
  theme(text = element_text(size = 18),
        legend.position = "none")
sankey.caf
ggsave(plot = sankey.caf, height = 10.5,
       filename = paste0(out.dir, "sankey_TMEinteractions.pdf"))

weights.igf <- weights %>%
  filter(TME != "CAFs") %>%
  bind_rows(weights.igf) %>%
  mutate(TME = if_else(TME == "CAFs", "CAFs IGF+", TME)) %>%
  rename(tumour.tcs = tumour) %>%
  group_by(TME) %>%
  mutate(sum.value = sum(value)) %>%
  ungroup() %>%
  filter(sum.value > 0) %>%
  mutate(percentage = round((value/sum.value) * 100), digits = 0)

ocurrences.igf <- bctumour@meta.data %>%
  select(patient, subclone, tumour.tcs) %>%
  filter(subclone != "diploid") %>%
  mutate(subclone = paste0(patient, " (SC", subclone, ")")) %>%
  select(subclone, tumour.tcs) %>%
  left_join(weights.igf, by = c("tumour.tcs"), relationship = "many-to-many")

ocurrences.igf <- lapply(ocurrences.igf, rep, ocurrences.igf$percentage) |>
  bind_rows()

sankey.igf <- ocurrences.igf %>%
  make_long(subclone, tumour.tcs, TME) %>%
  ggplot(aes(x = x, next_x = next_x, node = node, next_node = next_node,
             fill = factor(node), label = node)) +
  geom_sankey(flow.alpha = 0.5, node.color = "black", show.legend = FALSE) +
  geom_sankey_label(size = 3, color = "black") +
  scale_fill_manual(values = colour.all) +
  theme_void() +
  theme(text = element_text(size = 18),
        legend.position = "none")
sankey.igf
ggsave(plot = sankey.igf, height = 10.5,
       filename = paste0(out.dir, "sankey_TMEinteractions_IGF.pdf"))

# Fisher's exact test
stat.list <- ocurrences.igf %>%
  count(tumour.tcs, TME) %>%
  spread(tumour.tcs, n) %>%
  mutate(TC1.1 = if_else(is.na(TC1.1), 0, TC1.1),
         TC1.2 = if_else(is.na(TC1.2), 0, TC1.2),
         TC2 = if_else(is.na(TC2), 0, TC2),
         TC3 = if_else(is.na(TC3), 0, TC3)) %>%
  column_to_rownames("TME")

stat.list %>%
  mutate(pval = fisher_test(stat.list, simulate.p.value = TRUE)$p) %>%
  adjust_pvalue(method = "fdr") %>%
  select(-pval)
