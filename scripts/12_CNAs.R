# conda activate pub-cna
rm(list = ls()) # R version 4.2.3 (2023-03-15)
library(tidyverse) # tidyverse_2.0.0
library(GenomicRanges) # GenomicRanges_1.50.0
set.seed(1)
out.dir <- "results/tumour/ROI/"

# --- Data ---
# Gene segments
load("data/clones/V19L29/_count_mtx_annot.RData")

# --- Code ---
# Create outdir
dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# Get the gene CNA for each V19L29 subclone
subclones <- list.files("data/clones/V19L29/", pattern = "[0-9]{1}_CN.seg") %>%
  str_extract(pattern = "_.*_") %>%
  str_remove_all(pattern = "(_)|(subclone)")

cna <- lapply(subclones, FUN = function(i) {
  # Read segments
  segments <- read.table(paste0("data/clones/V19L29/_subclone", i, "_CN.seg"))
  # Process gene and segments data
  segments <- segments %>%
    mutate(segment_id = paste(Chr, Pos, End, sep = "_"))
  rownames(segments) <- segments$segment_id
  # Create GRanges objects
  genes.gr <-
    makeGRangesFromDataFrame(count_mtx_annot, ignore.strand = TRUE,
                             keep.extra.columns = TRUE)
  segments.gr <-
    makeGRangesFromDataFrame(segments, ignore.strand = TRUE,
                             keep.extra.columns = TRUE, seqnames.field = "Chr",
                             start.field = "Pos", end.field = "End")
  # Filter out genes that overlap in > or < than 1 region
  one.region <- countOverlaps(genes.gr, segments.gr) == 1
  genes.gr <- genes.gr[names(one.region)[one.region]]

  # Find overlap between regions
  overlap <- findOverlaps(genes.gr, segments.gr, ignore.strand = FALSE)

  # Merge data frames according to the overlap
  merged.coords <- data.frame(gene_id = names(genes.gr)[overlap@from],
                              segment_id = names(segments.gr)[overlap@to]) %>%
    left_join(count_mtx_annot, by = "gene_id") %>%
    left_join(segments, by = "segment_id") %>%
    select(gene_id, gene_name, Chr, Pos, End, CN, segm.mean) %>%
    unique() %>%
    mutate(subclone = i)
  return(merged.coords)
}) |>
  bind_rows()

# Save
write.table(cna, file = paste0(out.dir, "CNA_ROIs.tsv"), sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
