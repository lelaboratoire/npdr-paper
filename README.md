# Nearest-neighbor Projected-Distance Regression (NPDR) detects network interactions and controls for confounding and multiple testing

## Repository description

This repository provides simulation and analysis code to reproduce the results in our NPDR paper.
The repository is structured as followed:

- The `r` directory contains the `.Rmd` files to perform the analysis and produce the figures. 
Each file can be run separately in the order indicated in the file name.
`*_visualize_*.Rmd` files will have notes on which files need to be run prior to executing chunks.
`2_*_100_*.Rmd` files produces results from 100 replications.
- The `data` directory contains real and simulated data from running said `.Rmd` files.
- The `results` directory contains `.Rdata` results from running code in the `r` directory.
`*.Rdata` result files can be loaded in to `R` using the `load()` function.
`*.csv` result files should be straightforward.
- The `figs` directory contains generated figures from running the `*_visualize_*.Rmd` code in the `r` directory.
These figures are included in the main text `ms/main_npdr.pdf` or supplement `ms/sup_npdr.pdf`.
- The `ms` directory contains the text files to generate the manuscript and supplement (with `Sweave`).

## Abstract

We develop a new feature selection technique that uses the generalized linear model (GLM) to perform regression between nearest-neighbor pair distances projected onto attributes to address a broad spectrum of statistical challenges in high-dimensional data, including interactions and confounding variables.
Recently we developed STatistical Inference Relief (STIR), a pseudo t-test approach to estimate the statistical significance of Relief-based attribute scores for case-control (classification) problems, where the data may involve main effects and complex statistical interactions (network epistasis).
However, efficient statistical inference methods are needed to detect complex network effects in more complicated modeling situations, including continuous outcomes, mixtures of categorical and continuous predictors, and correcting for potential confounding variables.
Our new Nearest-neighbors Projected-Distance Regression (NPDR) encompasses STIR for case-control data and extends its capabilities to statistically correct for covariates, which previously was not feasible for Relief-based methods.
NPDR provides a novel and improved way to compute regression-based Relief scores (RRelief) for quantitative outcomes that also allows statistical corrections and various combinations of predictor data types such as genetic variants and gene expression.
In addition, we implement a penalized version of NPDR and derive the theoretical constant-$k$ approximation to the expected number of neighbors for spatially uniform radial neighborhoods.
Using realistic simulations that include main effects and gene-network interactions, we show for that NPDR improves attribute estimates compared to standard Relief while also being able to compute statistical significance of attributes and adjust for covariates.
We demonstrate that the NPDR framework has similar statistical properties to the pseudo t-test based method for case-control data with the added flexibility of regression modeling. We also compare glmnet with the penalized version of NPDR.
We use RNA-Seq data from a study of major depressive disorder to show that NPDR with covariate adjustment removes spurious associations due to confounding by sex.

## The NPDR `R` package
Installing:
```
devtools::install_github("insilico/npdr") 
library(npdr)
```

Github page:
[https://github.com/insilico/npdr](https://github.com/insilico/npdr)