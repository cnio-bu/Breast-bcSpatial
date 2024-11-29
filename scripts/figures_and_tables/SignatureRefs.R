# conda activate pub-figures
rm(list = ls()) # R version 4.3.1 (2023-06-16)
library(tidyverse) # tidyverse_2.0.0
library(beyondcell) # beyondcell_2.2.0
set.seed(1)
out.dir <- "tables/"

# --- Data ---
# Functional signatures
reactome <- readGMT("data/signatures/c2.cp.reactome.v2023.1.Hs.symbols.gmt")
hallmarks <- readGMT("data/signatures/h.all.v2023.1.Hs.symbols.gmt")
interesting <- readGMT("data/signatures/functional_signatures.gmt")
functional.sigs <- c(reactome, hallmarks, interesting)

# --- Code ---
# Create outdir
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

# References for functional signatures
hallmarks.ref <- "Human MSigDB Hallmarks v2023.1 (https://www.gsea-msigdb.org/gsea/downloads_archive.jsp)"
reactome.ref <- "Human MSigDB Reactome v2023.1 (https://www.gsea-msigdb.org/gsea/downloads_archive.jsp)"
cancersea.ref <- "CancerSEA (http://biocc.hrbmu.edu.cn/CancerSEA/goDownload; date: 02/05/2022)"
groger.ref <- "Groger et al. 2012 (doi: 10.1371/journal.pone.0051136)"
schuetz.ref <- "Human MSigDB CGP v2023.1 (https://www.gsea-msigdb.org/gsea/downloads_archive.jsp)"
ayers.ref <- "Ayers et al. 2017 (doi: 10.1172/JCI91190)"
schroth.ref <- "Schroth et al. 2020 (doi: 10.1158/1078-0432.CCR-20-1923)"
wang.ref <- "Wang et al. 2023 (doi: 10.1016/j.celrep.2023.112338)"
buffa.ref <- "Human MSigDB CGP v2023.1 (https://www.gsea-msigdb.org/gsea/downloads_archive.jsp)"
cabrita.ref <- "Cabrita et al. 2020 (doi: 10.1038/s41586-019-1914-8)"

# Create table
df.sigs <- data.frame(sigID = names(functional.sigs)) %>%
  mutate(Signature =
    str_remove_all(sigID, pattern = "(^sig_)|(_UP)|(_DOWN)")) %>%
  add_count(Signature) %>%
  mutate(Signature = if_else(n > 1, str_remove(sigID, pattern = "^sig_"),
                             Signature),
        Signature = str_replace_all(Signature, pattern = "_", replacement = " "),
        Signature =
          str_remove_all(Signature,
                         pattern = paste0("(^HALLMARK)|(^REACTOME)|",
                                          "(CancerSEA$)|(BUFFA)")),
        Signature = gsub("\\s+", " ", str_trim(toupper(Signature))),
        Signature =
          case_when(startsWith(sigID, prefix = "HALLMARK") ~
                      paste(Signature, "(Hallmark)"),
                    sigID == "REACTOME_PD_1_SIGNALING_UP" ~
                      "PD-1 SIGNALING (Reactome)",
                    sigID ==
                      paste0("REACTOME_APC_C_CDH1_MEDIATED_DEGRADATION_OF_",
                             "CDC20_AND_OTHER_APC_C_CDH1_TARGETED_PROTEINS_",
                             "IN_LATE_MITOSIS_EARLY_G1_UP") ~
                      paste("APC C CDH1 MEDIATED DEGRADATION OF CDC20 AND",
                            "OTHER PROTEINS IN LATE MITOSIS EARLY G1",
                            "(Reactome)"),
                    sigID ==
                      paste0("REACTOME_ANTIGEN_ACTIVATES_B_CELL_RECEPTOR_BCR_",
                             "LEADING_TO_GENERATION_OF_SECOND_MESSENGERS_UP") ~
                      paste("ANTIGEN ACTIVATES BCR LEADING TO GENERATION OF",
                            "SECOND MESSENGERS (Reactome)"),
                    startsWith(sigID, prefix = "REACTOME") ~
                      paste(Signature, "(Reactome)"),
                    endsWith(sigID, suffix = "CancerSEA_UP") ~
                      paste(Signature, "(CancerSEA)"),
                    sigID == "sig_EMT_GROGER_2012_UP" ~
                      "EMT UP (Groger et al. 2012)",
                    sigID == "sig_EMT_GROGER_2012_DOWN" ~
                      "EMT DOWN (Groger et al. 2012)",
                    sigID == "SCHUETZ_BREAST_CANCER_DUCTAL_VS_INVASIVE_UP" ~
                      "DUCTAL BREAST CANCER (Schuetz et al. 2006)",
                    sigID == "SCHUETZ_BREAST_CANCER_DUCTAL_VS_INVASIVE_DOWN" ~
                      "INVASIVE BREAST CANCER (Schuetz et al. 2006)",
                    sigID == "TIS_UP" ~
                      paste(Signature, "(Ayers et al. 2017)"),
                    sigID == "BRCANESS_UP" ~
                    paste(Signature, "(Schroth et al. 2020)"),
                    startsWith(sigID, prefix = "ECM") ~
                      paste(Signature, "(Wang et al. 2023)"),
                    startsWith(sigID, prefix = "BUFFA") ~
                      paste(Signature, "(Buffa et al. 2010)"),
                    sigID == "TLS_UP" ~
                      paste(Signature, "(Cabrita et al. 2020)")),
        Reference =
          case_when(startsWith(sigID, prefix = "HALLMARK") ~ hallmarks.ref,
                    startsWith(sigID, prefix = "REACTOME") ~ reactome.ref,
                    endsWith(sigID, suffix = "CancerSEA_UP") ~ cancersea.ref,
                    startsWith(sigID, prefix = "sig_EMT_GROGER") ~ groger.ref,
                    startsWith(sigID, prefix = "SCHUETZ") ~ schuetz.ref,
                    sigID %in% c("TIS_UP", "BRCANESS_UP") ~ schroth.ref,
                    startsWith(sigID, prefix = "ECM") ~ wang.ref,
                    startsWith(sigID, prefix = "BUFFA") ~ buffa.ref,
                    sigID == "TLS_UP" ~ cabrita.ref)) %>%
  select(-n)

up.down <-
  data.frame(sigID = c("SCHUETZ_BREAST_CANCER_DUCTAL_VS_INVASIVE",
                       "sig_EMT_GROGER_2012"),
            Signature = c("DUCTAL VS INVASIVE (Schuetz et al. 2006)",
                          "EMT (Groger et al. 2012)"),
            Reference = c(schuetz.ref, groger.ref))
df.sigs <- df.sigs %>%
  rbind(up.down) %>%
  arrange(Signature)

write.table(df.sigs, file = paste0(out.dir, "SuppTable7.tsv"),
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
