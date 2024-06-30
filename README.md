# SCIPAC
**Combining single-cell and bulk RNA-sequencing data to identify phenotype-associated cells**

## Introduction
Numerous algorithms have been proposed to identify cell types in single-cell RNA sequencing data, yet a fundamental problem remains: determining associations between cells and phenotypes such as cancer. We develop SCIPAC (Single- Cell and bulk data-based Identifier for Phenotype Associated Cells), the first algorithm that quantitatively estimates the association between each cell in single-cell data and a phenotype. SCIPAC also provides a p-value for each association and applies to data with virtually any type of phenotype, for instance, binary (e.g., cancer vs. normal), ordinal (e.g., different stages of cancer), continuous (e.g., quantitative traits), or survival. SCIPAC also requires minimum tuning and is computationally very fast.

## Updates
* Jun, 2024: minor updates in the package due to the changes in the dependencies.
* Nov, 2022: SCIPAC version 1.0.0 is launched.

## Installation
* System requirements: SCIPAC is developed under R (version >= 4.1.2).
* Latest version: the latest developmental version of SCIPAC can be downloaded from GitHub and installed from source by `devtools::install_github('RavenGan/SCIPAC')`.

## Manual
Please see [http://RavenGan.github.io/SCIPAC/vignettes/introduction.html](http://RavenGan.github.io/SCIPAC/vignettes/introduction.html) for details. In the introduction, we also include how to use TCGA databases as the source of bulk RNA-seq data for the application of SCIPAC. In the R terminal, users can also use `?SCIPAC` to access the help documents.

## Examples
In the [SCIPAC tutorial](http://RavenGan.github.io/SCIPAC/vignettes/introduction.html), we use multiple examples to show the function and performance of SCIPAC. Example datasets are also included in SCIPAC package.

## Citation
If you use SCIPAC in your research, please cite the following paper:

Gan, D., Zhu, Y., Lu, X., & Li, J. (2024). SCIPAC: quantitative estimation of cell-phenotype associations. Genome Biology, 25(1), 119.