# conda activate pub-figures
# install.packages("figpatch")
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
library(rstatix) # rstatix_0.7.2
library(ggpubr) # ggpubr_0.6.0
library(patchwork) # patchwork_1.1.3
set.seed(1)
out.dir <- "figures/reviewers/benchmark/normalisation/"

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

# Reads, merges BCS results and regress patient effect
bcRead <- function(method = c("giotto", "scanpy"), 
                   mode = c("functional", "sensitivity")) {
  # Create a list of Beyondcell objects
  bcs <- lapply(patients, FUN = function(p) {
    bcobj <- readRDS(paste0("data/benchmarking/normalisation/", method,
                            "/", mode, "/", p, "_", method, "_bcs_", mode, ".rds"))
    bcobj@meta.data <- bcobj@meta.data %>%
      rownames_to_column("spots") %>%
      select(spots) %>%
      left_join(metadata, by = "spots") %>%
      column_to_rownames("spots")
    return(bcobj)
  })
  # Merge Beyondcell objects
  bcmerged <- reduce(bcs, bcJoin)
  # Replace NAs by 0s
  bcmerged@normalized[is.na(bcmerged@normalized)] <- 0
  bcmerged <- bcRecompute(bcmerged, slot = "normalized")
  # Regress patient effect
  bcregressed <- bcRegressOut(bcmerged, vars.to.regress = "patient")
  return(bcregressed)
}

# --- Data ---
# Metadata
metadata <- readRDS("results/all/sensitivity/tcs_bcs.rds")@meta.data %>%
  rownames_to_column("spots")

# Patient IDs
patients <- unique(metadata$patient)

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Compute correlations for functional and sensitivity BCS
corr.plots <- lapply(c("functional", "sensitivity"), FUN = function(mode) {
  # BCS results
  if (mode == "functional") {
    seurat <- readRDS("results/all/functional/functional_bcs_regressed.rds")
  }
  else if (mode == "sensitivity") {
    seurat <- readRDS("results/all/sensitivity/tcs_bcs.rds")
  } 
  giotto <- bcRead("giotto", mode = mode)
  scanpy <- bcRead("scanpy", mode = mode)
  # Pivot to long format
  seurat.long <- seurat@normalized %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spots") %>%
    pivot_longer(cols = rownames(seurat@normalized),
                 names_to = "signature", values_to = "seurat.enrichment")
  giotto.long <- giotto@normalized %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spots") %>%
    pivot_longer(cols = rownames(giotto@normalized),
                 names_to = "signature", values_to = "giotto.enrichment")
  scanpy.long <- scanpy@normalized %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("spots") %>%
    pivot_longer(cols = rownames(scanpy@normalized),
                 names_to = "signature", values_to = "scanpy.enrichment")
  normalised.long <- seurat.long %>%
    left_join(giotto.long, by = c("spots", "signature")) %>%
    left_join(scanpy.long, by = c("spots", "signature"))
  # Compute Pearson correlations between enrichment scores
  comparisons <-
    list("Seurat-Giotto" = paste(c("seurat", "giotto"), "enrichment", sep = "."),
         "Seurat-Scanpy" = paste(c("seurat", "scanpy"), "enrichment", sep = "."),
         "Giotto-Scanpy" = paste(c("giotto", "scanpy"), "enrichment", sep = "."))
  corr <- lapply(comparisons, FUN = function(x) {
    normalised.long %>%
      cor_test(vars = x[1], vars2 = x[2], method = "pearson")
  }) |>
    bind_rows() %>%
    adjust_pvalue(method = "fdr")
  # Downsample spots
  normalised.downsampled <- normalised.long %>%
    sample_n(500)
   values <- normalised.downsampled %>%
    select(seurat.enrichment:scanpy.enrichment) %>%
    as.matrix()
  limits <- pretty(c(min(values), max(values)))
  limits <- c(limits[1], limits[length(limits)])
  # Scatter plots
  corr.plot <- lapply(comparisons, FUN = function(x) {
    df <- corr %>% filter(var1 == x[1] & var2 == x[2])
    sp <- ggscatter(normalised.downsampled, x = x[1], y = x[2], 
                    add = "reg.line", conf.int = TRUE,
                    add.params = list(color = "blue", fill = "lightgray"))
    sp <- ggpar(sp, xlim = limits, ylim = limits)
    test <- stat_cor(label.x = limits[1]+10, label.y = limits[2]-10,
                     aes(label = paste0("r = ", round(df$cor, 2),
                                        "; FDR = ", df$p.adj)))
    full.plot <- sp + test +
      labs(x = str_to_title(str_remove(x[1], pattern = "\\..*$")),
           y = str_to_title(str_remove(x[2], pattern = "\\..*$"))) +
      theme(text = element_text(size = 18))
    return(full.plot)
  }) |>
    wrap_plots(ncol = 3) +
    plot_annotation(title = paste(str_to_title(mode), "BCS correlation")) &
    theme(text = element_text(size = 18),
          plot.title = element_text(size = 22, hjust = 0.5, face = "bold"),
          legend.position = "bottom") 
  corr.plot
}) |>
  wrap_plots(nrow = 2) + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 30))
corr.plots
ggsave(plot = corr.plots, width = 20, height = 14, 
      filename = paste0(out.dir, "correlations.pdf"))
