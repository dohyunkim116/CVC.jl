# CVC.jl

| **Documentation** | **Build Status** | **Code Coverage**  |
|-------------------|------------------|--------------------|
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://dohyunkim116.github.io/CVC.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://dohyunkim116.github.io/CVC.jl/dev/) | [![Build Status](https://github.com/dohyunkim116/CVC.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dohyunkim116/CVC.jl/actions/workflows/CI.yml?query=branch%3Amain)  | [![codecov](https://codecov.io/gh/dohyunkim116/CVC.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dohyunkim116/CVC.jl) |  

The CVC model stands for Censored Variance Component model. It is a variant of Haseman-Elston (HE) regression for estimating partitioned heritability of censored time-to-event outcomes. It leveerages synthetic variable and randomized trace estimator for handling right-censored data and analyzing large-scale genetic data from biobank.

## Key Features

- **Time-to-event (survival) data analysis** with right-censoring via synthetic variable (Leurgans and Koul-Susarla-van Ryzin estimators)
- **Partitioned heritability estimation** across K genetic components (e.g., chromosomes and gene sets)
- **Computational efficiency** through memory-mapped PLINK files, randomized trace estimator, and column-wise block jackknife standard error estimators.
- **Optional non-negative constrained heritability estimates** via NNLS optimization

## Installation

This package requires Julia v1.5 or later. See the documentation for tutorial. 

The package is currently under development. To install, start Julia and use the `]` key to switch to the package manager REPL:

```julia
(@v1.9) pkg> add https://github.com/dohyunkim116/CVC.jl
```

Use the backspace key to return to the Julia REPL.

## Quick Example

```julia
using CVC, DelimitedFiles

# Load data
w = readdlm("covariates.txt")  # N × p covariate matrix
u = vec(readdlm("observed_times.txt"))  # observed event times
δ = Bool.(vec(readdlm("censoring_indicators.txt")))  # censoring indicators
temprary_directory = tempdir() # for working arrays

# Fit time-to-event model with K genetic components
cvcm = cvc(u, δ, w, "path/to/genotypes", temporary_directory)
fit_me!(cvcm)

# View partitioned heritability estimates
println("Total h²: ", cvcm.h2[], " (SE: ", cvcm.h2se[], ")")
println("Component h²: ", cvcm.h2k)  # K+1 components (K genetic + 1 environmental)
```

## Data Requirements

CVC.jl expects genotype data in partitioned PLINK format (`.bed`/`.bim`/`.fam` files):
- Partition genotypes across K components: `G1.bed`, `G2.bed`, ..., `GK.bed`
- Each partition should have corresponding `.bim` and `.fam` files
- Phenotypes and genotypes must be aligned (same sample ordering)

## Citation

If you use CVC.jl in your research, please cite:

*Kim D, Zhou H, Chau B, Jensen A, Shen J, Mehrotra D, Li G, Zhou J. Estimating Heritability of Survival Traits Using Censored Multiple Variance Component Model. arXiv. 2025. doi: 10.48550/ARXIV.2510.26226.*

## Acknowledgments

This work was supported by National Institutes of Health under grants P30
CA-16042, UL1TR000124-02, P50CA211015, T32 HG002536, R35 GM141798, R01 HG006139, and R01 DK142026; National Science Foundation under grants DMS 2054253 and IIS 2205441.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
