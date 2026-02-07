# CVC.jl Example Data

This directory contains a toy dataset for demonstrating CVC.jl functionality. These data are synthetically generated for tutorial purposes only.

## Directory Structure

```
data/
├── README.md (this file)
├── G1.bed/bim/fam
├── G2.bed/bim/fam
├── ...
├── G10.bed/bim/fam
├── u.txt                 # Observed event times
├── delta.txt             # Censoring indicators (1=event, 0=censored)
└── w.txt                 # Covariate matrix
```

## Example

### Description
- **Sample size (N)**: 5000
- **Number of genetic components (K)**: 10
- **Total SNPs (M)**: 10000 (1000 per component)
- **Covariates (C)**: 10
- **Censoring rate**: 0.20
- **Target heritability (h2)**: 0.50
- **LD correlation (rho)**: 0.03

### File Formats

**u.txt**: N x 1 text file with observed event times
```
0.935533873943533
-10.656942253227506
1.1283210895382028
...
```

**delta.txt**: N x 1 text file with censoring indicators
```
0  # Event observed
0  # Censored
1  # Event observed
...
```

**w.txt**: N x C text file of covariates

**Genotypes**: Partitioned PLINK files
- `G1.bed/bim/fam` through `G10.bed/bim/fam`

### Usage Example

```julia
using CVC, DelimitedFiles

# Load data
data_dir = joinpath(pwd(), "data")
u = vec(readdlm(joinpath(data_dir, "u.txt")))
delta = Bool.(vec(readdlm(joinpath(data_dir, "delta.txt"))))
w = readdlm(joinpath(data_dir, "w.txt"))

# Fit model
cvcm = cvc(u, delta, w, data_dir, tempdir())
fit_me!(cvcm)

# Results
println("h² = ", cvcm.h2[], " ± ", cvcm.h2se[])
```

## Data Generation Script

The data were generated using a CVCData-based script with these steps:
- Simulate genotypes with `simulate_geno` using K=10 components and M=10000 SNPs
- Simulate C=10 covariates with `simulate_cov`
- Build the genetic mean component and simulate right-censored outcomes
- Write `u.txt`, `delta.txt`, and `w.txt` to this folder

## Data Generation Details

### Genotype Simulation
- MAFs are drawn from Beta(2,8), then scaled to the range [0.01, 0.50]
- LD structure uses AR(1) with rho = 0.03

### Phenotype Simulation
- Right-censored outcomes are generated using `simulate_right_censored_data`
- Censoring rate is set to 0.20

### Covariate Generation
- Covariates are simulated with `simulate_cov` (C = 10)

## Notes

- **Not for publication**: These are toy datasets for learning/testing only
- **Synthetic data**: No real individual or genetic information
- **Small scale**: Real analyses typically have N>100,000
- **Simplified LD**: Real data has complex LD structure
- **No relatedness**: Assumes unrelated individuals

## Citation

If you use these example datasets in publications or presentations, please cite:

*Kim D, Zhou H, Chau B, Jensen A, Shen J, Mehrotra D, Li G, Zhou J. Estimating Heritability of Survival Traits Using Censored Multiple Variance Component Model. arXiv. 2025. doi: 10.48550/ARXIV.2510.26226.*

## Contact

For questions about the example data or to report issues:
- GitHub Issues: [github.com/dohyunkim116/CVC.jl/issues](https://github.com/dohyunkim116/CVC.jl/issues)
- Email: [dohyunkim.116@gmail.com]
