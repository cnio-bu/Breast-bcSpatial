# conda activate pub-beyondcell
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
set.seed(1)
out.dir <- "results/all/"

# --- Data ---
# Merged BCS with TCs
bcmerged <- readRDS("results/all/sensitivity/tcs_bcs.rds")

# MoAs
moas <- read.table("data/signatures/MoAs.tsv", header = TRUE, sep = "\t")

# --- Code ---
# Create outdir
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# Subset Beyondcell object to remove the interphase spots
spots.not.interphase <- bcmerged@meta.data %>%
  mutate(major.tcs = if_else(major.tcs == "Mixed", "Interphase", major.tcs)) %>%
  filter(major.tcs != "Interphase") %>%
  rownames()

bcsubseted <- bcSubset(bcmerged, cells = spots.not.interphase)

# Retrieve the most specific drugs to target the Tumour or TME compartments
bcsubseted <- bcRanks(bcsubseted, idents = "major.tcs",
                      resm.cutoff = c(0.05, 0.95), extended = FALSE)

diff.response <- c("TOP-Differential-HighSensitivity")

response <- bcsubseted@ranks$major.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("group"), names_to = "major_TC",
               values_to = "response", names_prefix = "group.") %>%
  filter(response %in% diff.response) %>%
  select(IDs, major_TC, response)

switch.point <- bcsubseted@ranks$major.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("switch.point"), names_to = "major_TC",
               values_to = "switch_point", names_prefix = "switch.point.") %>%
  select(IDs, major_TC, switch_point)

mean.bcs <- bcsubseted@ranks$major.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("mean"), names_to = "major_TC",
               values_to = "mean", names_prefix = "mean.") %>%
  select(IDs, major_TC, mean)

residuals.mean <- bcsubseted@ranks$major.tcs %>%
  rownames_to_column("IDs") %>%
  pivot_longer(cols = starts_with("residuals.mean"), names_to = "major_TC",
               values_to = "residuals_mean", names_prefix = "residuals.mean.") %>%
  select(IDs, major_TC, residuals_mean)

ranking <- response %>%
  left_join(switch.point, by = c("IDs", "major_TC")) %>%
  left_join(mean.bcs, by = c("IDs", "major_TC")) %>%
  left_join(residuals.mean, by = c("IDs", "major_TC")) %>%
  rename(drugID = IDs)

# Add drug names and MoAs
ranking <- ranking %>%
  left_join(moas, by = "drugID") %>%
  mutate(study = str_remove_all(str_extract(drugID, pattern= "_.*_"),
                                pattern = "_"),
        drug_name = toupper(paste0(drug_name, " (", study, ")"))) %>%
  select(drugID, drug_name, MoA, major_TC, response, switch_point, mean,
         residuals_mean)

write.table(ranking, file = "tables/SuppTable8.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Save
saveRDS(bcsubseted, file = paste0(out.dir, "sensitivity/ranked_bcs.rds"))
