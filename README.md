# Nearest-neighbor Projected-Distance Regression (NPDR) detects network interactions and controls for confounding and multiple testing

## Repository description

This repository provides simulation and analysis code to reproduce the results in our NPDR paper.
The repository is structured as followed:

- The `r` directory contains the `.Rmd` files to perform the analysis and produce the figures. Each file can be run separately in arbitrary order.
- The `data` directory contains real and simulated data from running said `.Rmd` files.
- The `results` directory contains `.Rdata` results from running code in the `r` directory.
- The `figs` directory contains generated figures from running code in the `r` directory.
- The `ms` directory contains the text files to generate the manuscript and supplement (with `Sweave`).

