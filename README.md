# SCIPAC
**Combining single-cell and bulk RNA-sequencing data to identify phenotype-associated cells**

## Introduction
SCIPAC (Single- Cell and bulk data-based Identifier for Phenotype Associated Cells) is a computational method that identifies cells in single-cell data that are associated with a given phenotype. This phenotype can be binary (e.g., cancer vs. normal), ordinal (e.g., different stages of cancer), continuous (e.g., quantitative traits), or survival.

## Updates
* Nov, 2022: SCIPAC version 1.0.0 is launched

## Installation
* System requirements: SCIPAC is developed under R (version >= 4.1.2).
* Latest version: the latest developmental version of SCIPAC can be downloaded from GitHub and installed from source by `devtools::install_github('RavenGan/SCIPAC')`.

## Manual
Please see [http://RavenGan.github.io/SCIPAC/vignettes/introduction.html](http://RavenGan.github.io/SCIPAC/vignettes/introduction.html) for details. In the introduction, we also include how to use TCGA databases as the source of bulk RNA-seq data for the application of SCIPAC. In the R terminal, users can also use `?SCIPAC` to access the help documents.

## Examples
In the [SCIPAC tutorial](http://RavenGan.github.io/SCIPAC/vignettes/introduction.html), we use multiple examples to show the function and performance of SCIPAC. Example datasets are also included in SCIPAC package.
