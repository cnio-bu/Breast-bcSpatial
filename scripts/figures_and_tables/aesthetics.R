# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
set.seed(1)

# --- Patient order ---
patient.order <-
  c("CID4290", "CID4535", "1142243F", "1160920F", "CID4465", "CID44971",
    "738811QB", "1168993F", "V19L29")

# --- Ggplot themes ---
theme_neighbours <- theme_bw() +
  theme(text = element_text(size = 18),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

pub_theme_facet <- theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.spacing = unit(1, "lines"))

# --- Colour palettes ---
# Histopathologicl annotations
colour.histo <-
  c("Normal + stroma +\nlymphocytes" = "#7FC97F",
    "Stroma +\nadipose tissue" = "#FDC086", "Stroma" = "#FFFF99",
    "DCIS" ="#386CB0", "Lymphocytes" = "#F0027F",
    "Invasive cancer +\nlymphocytes" = "#BF5B17")

# Tumour and TME annotations
colour.annotations <- c("Tumour" = "#FF5D5DFF", "TME" = "#45B7C4")

# Subclones
colour.subclones <-
  c("1" = "#8DD3C7", "2" = "#FFFFB3", "3" = "#FB8072", "4" = "#80B1D3",
    "5" = "#FDB462", "Diploid" = "#B3DE69", "6" = "#FCCDE5")

# Breast cancer subtypes
colour.subtype <-
  c("Luminal" = "#C70039", "HER2+" = "#7BC57D", "TNBC" = "#A986FF")

# Breast cancer intrinsic subtypes
colour.molsubtypes <-
c("LumA" = "#FF799F", "LumB" = "#DB003F", "HER2+" = "#7BC57D",
 "Basal" = "#A986FF", "Undetermined" = "#D7D7D7")

# Cell types
colour.celltypes <-
  c("B cells" = "#29AFAA", "CAFs" = "#F3D000", "Cancer Epithelial" = "#E48095",
    "Endothelial" = "#EB5528", "Myeloid" = "#8CC2F3",
    "Normal Epithelial" = "#C9ADD9", "Plasmablasts" = "#FFA42B",
    "PVL" = "#A9DB8E", "T cells" = "#6E66D4")

# Major TCs
colour.tctypes <-
  c("Tumour-rich" = "#FF5D5DFF", "TME-rich" = "#45B7C4",
    "Interphase" = "#C2228AFF", "Mixed" = "#C2228AFF")

# Initial tumour TCs
colour.initial.tcs <-
  c("TC1" = "#C7D13A", "TC2" = "#FC8D62", "TC3" = "#F5699F", "TC4" = "#2998AC",
    "TC5" = "#8DA0CB")

# Final tumour TCs
colour.tcs <-
  c("TC1.1" = "#FFD92F", "TC1.2" = "#A6D854", "TC2" = "#FC8D62",
    "TC3" = "#F5699F")

# MoAs
colour.moas <-
  c("ALK inhibitor" = "#9C0D38",
  "Anti-inflammatory / Immunosuppressant" = "#B3446C",
  "BRAF inhibitor" = "#EA95AE", "Cell cycle arrest" = "#F19C79",
  "DNA related agent" = "#FF864C", "EGFR inhibitor" = "#FB6467",
  "HER2 inhibitor" = "#FFF154", "Hormonal therapy" = "#ADD82F",
  "HSP inhibitor" = "#5DA271", "JAK-STAT signaling inhibitor" = "#BCE09B",
  "Kinase inhibitor" = "#6EB4D1", "MAPK inhibitor" = "#1B9AAA",
  "Microtubule agent" = "#4DC8BB", "Multi-kinase inhibitor" = "#B1EDE8",
  "NFkB signaling inhibitor" = "#0067A5",
  "Non-apopototic cell death" = "#4B71E5", "Other" = "#D7D7D7",
  "PI3K/AKT/mTOR inhibitor" = "#D58ACA",
  "Pro-apoptotic agent" = "#BDADEA", "Statin" = "#B47EB3",
  "SYK inhibitor" = "#4F2561", "TGFBR inhibitor" = "#A44A3F",
  "VEGFR inhibitor" = "#CA9A8C", "WNT signaling inhibitor" = "#410200")

# ROIs
colour.rois <- c("ROI1" = "#F8766D", "ROI2" = "#B79F00", "ROI3" = "#00BA38",
                 "ROI4" = "#00BFC4", "ROI5" = "#619CFF", "ROI6" = "#F564E3")

# Expression clusters
colour.expression <-
  c("#FB9A99", "#56B4E9", "#E31A1C", "#D95F02", "#E6AB02", "#E7298A", "#66A61E",
    "#1F78B4", "#7570B3", "#1B9E77")

# CNVs
colour.cna <- c("2" = "#B3DE69", "3" = "#FF799F", "4" = "#F0027F")

# Save
save(list = ls(), file = "scripts/figures_and_tables/aesthetics.RData")
