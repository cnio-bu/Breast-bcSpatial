# breast-bcspatial
Code used to produce the final figures and tables of the study "Spatial Transcriptomics in Breast Cancer Reveals Tumour Microenvironment-Driven Drug Responses and Clonal Therapeutic Heterogeneity".

## Installation
Use the git clone command to create a local copy:

```
git clone https://github.com/cnio-bu/breast-bcspatial
```

## How to run

### Set up
You need to download additional data folders from Zenodo (DOI: 10.5281/zenodo.10638906) for the code to be functional:

* **`visium/`:** Contains processed spatial transcriptomics Seurat objects with deconvoluted spots, SCTransform-normalised counts, and clonal composition predicted with SCEVAN [1]. Please make sure to merge the contents of this folder with the `data/visium/scalefactors` folder that is provided in this repository.

* **`single-cell/`:** Contains raw and filtered merged single-cell RNA-seq Seurat objects with unnormalised counts used as a reference for spot deconvolution.

* **`beyondcell/sensitivity`:** Contains Beyondcell sensitivity objects with prediction scores for all drug response signatures in [SSc breast](https://github.com/cnio-bu/SSc-breast).

* **`beyondcell/functional`:** Contains Beyondcell functional objects with enrichment scores for all functional signatures.
  
* **`benchmarking/deconvolution`:** Spot-wise deconvolution according to CARD [2] and the spatialDWLS method [3] implemented in Giotto, two deconvolution tools that were compared to RCTD [4], our final selection.

* **`benchmarking/normalisation`:** Beyondcell sensitivity and functional objects computed using Scanpy normalisation with log-transformation [5] or Giotto normalisation with log-transformation and z-scoring [6]. These two methods were compared to Seurat SCTransform [7], our final selection.

These objects were generated with the code available at [cnio-bu/ST-preprocess](https://github.com/cnio-bu/ST-preprocess).

**References**

1. De Falco A, Caruso F, Su X-D, Iavarone A, Ceccarelli M. [A variational algorithm to detect the clonal copy number substructure of tumors from scRNA-seq data.](https://www.nature.com/articles/s41467-023-36790-9) *Nat Commun.* 2023;**14**(1):1074.
2. Ma Y, Zhou X. [Spatially informed cell-type deconvolution for spatial transcriptomics.](https://www.nature.com/articles/s41587-022-01273-7) *Nat Biotechnol.* 2022;**40**(9):1349-1359.
3. Dong R, Yuan GC. [SpatialDWLS: accurate deconvolution of spatial transcriptomic data.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02362-7) *Genome Biol.* 2021;**22**(1):145.
4. Cable DM, Murray E, Zou LS, Goeva A, Macosko EZ, Chen F, Irizarry RA. [Robust decomposition of cell type mixtures in spatial transcriptomics.](https://www.nature.com/articles/s41587-021-00830-w) *Nat Biotechnol.* 2022;**40**(4):517-526.
6. Wolf FA, Angerer P, Theis FJ. [SCANPY: large-scale single-cell gene expression data analysis.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1382-0) *Genome Biol.* 2018;**19**(1):15.
7. Dries R, Zhu Q, Dong R, Eng CL, Li H, Liu K, Fu Y, Zhao T, Sarkar A, Bao F, George RE, Pierson N, Cai L, Yuan GC. [Giotto: a toolbox for integrative analysis and visualization of spatial expression data.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02286-2) *Genome Biol.* 2021;**22**(1):78.
8. Hafemeister C, Satija R. [Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1) *Genome Biol.* 2019;**20**(1):296.

### Execution
Once all data is downloaded, just run the code in the `scripts/` folder in order. Then, run the code in the `scripts/figures_and_tables` folder to generate the final figures and tables. The reviewer's suggested analyses are stored in the `scripts/reviewers` folder.

## Authors

* María José Jiménez-Santos

<!-- ## Citation -->

## Support
If you have any questions, feel free to submit an [issue](https://github.com/cnio-bu/breast-bcspatial/issues).
