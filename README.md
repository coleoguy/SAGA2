# SAGA2

<!-- badges: start -->
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
<!-- badges: end -->

**Software for the Analysis of Genetic Architecture**

SAGA2 implements an information-theoretic approach to line cross analysis (LCA) of genetic architecture. It uses AICc-based model selection to explore all possible models of composite genetic effects (CGEs) contributing to variation among generation means, producing model-averaged parameter estimates and unconditional standard errors. SAGA2 automatically constructs the appropriate C-matrix from user-supplied breeding designs, supporting XY, XO, ZW, ZO, and NSC sex chromosome systems.

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("coleoguy/SAGA2")
```

## Quick Example

```r
library(SAGA2)

# Load example data (Tribolium castaneum cross productivity)
data(BanInf)

# Run model-averaged line cross analysis
result <- LCA(BanInf, SCS = "NSC")

# View model-averaged estimates and unconditional SEs
result$estimates

# Plot results with variable importance coloring
plot(result, min.vi = 0.3)

# Visualize model space
VisModelSpace(result)

# Traditional LCA plot of observed data
plotObserved(BanInf)
```

## Key Features

- **Automatic C-matrix construction** from user-supplied breeding designs
- **Full model space exploration** with AICc-based model selection
- **Model-averaged parameter estimates** with unconditional standard errors
- **Variable importance scores** for each composite genetic effect
- **Publication-quality plotting** with S3 plot method for genarch objects
- **Support for multiple sex chromosome systems** (XY, XO, ZW, ZO, NSC)
- **Environmental and gene-by-environment effects** via the `env` argument
- **Observed parental effects** as an alternative to calculated maternal effects

## Citation

If you use SAGA2 in your research, please cite:

```r
citation("SAGA2")
```

> Blackmon, H. and Demuth, J.P. (2016). An information-theoretic approach to estimating the composite genetic effects contributing to variation among generation means: Moving beyond the joint-scaling test for line cross analysis. *Evolution*, 70(2), 420--432. doi:10.1111/evo.12867

## Issues

Please report bugs or request features at: https://github.com/coleoguy/SAGA2/issues
