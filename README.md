# Breast-bcSpatial
Code used to produce the final figures and tables of the study "Spatial Transcriptomics in Breast Cancer Reveals Tumour Microenvironment-Driven Drug Responses and Clonal Therapeutic Heterogeneity".

## Installation
Use the git clone command to create a local copy:

```
git clone https://github.com/cnio-bu/Breast-bcSpatial
```

## How to run

### Set up
You need to download additional data folders from Zenodo (DOI: 10.5281/zenodo.10638906) for the code to be functional:

* **`visium/`:** Contains processed spatial transcriptomics Seurat objects with deconvoluted spots, SCTransform-normalised counts, and clonal composition predicted with SCEVAN [1]. Please make sure to merge the contents of this folder with the `data/visium/scalefactors` folder that is provided in this repository.

* **`single-cell/`:** Contains raw and filtered merged single-cell RNA-seq Seurat objects with unnormalised counts used as a reference for spot deconvolution.

* **`beyondcell/sensitivity`:** Contains Beyondcell sensitivity objects with prediction scores for all drug response signatures in [SSc breast](https://github.com/cnio-bu/drug_susceptibility_collection/tree/breast).

* **`beyondcell/functional`:** Contains Beyondcell functional objects with enrichment scores for all functional signatures.

These objects were generated with the code available at [cnio-bu/bcSpatial](https://github.com/cnio-bu/bcSpatial).

**References**

1. De Falco A, Caruso F, Su X-D, Iavarone A, Ceccarelli M. [A variational algorithm to detect the clonal copy number substructure of tumors from scRNA-seq data.](https://www.nature.com/articles/s41467-023-36790-9) *Nat Commun.* 2023;**14**:1074.

### Execution
Once all data is downloaded, just run the code in the `scripts/` folder in order. Then, run the code in the `scripts/figures_and_tables` folder to generate the final figures and tables.

## Authors

* María José Jiménez-Santos

<!-- ## Citation -->

## Support
If you have any questions, feel free to submit an [issue](https://github.com/cnio-bu/Breast-bcSpatial/issues).
