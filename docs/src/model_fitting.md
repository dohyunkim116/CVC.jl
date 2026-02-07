# Model Fitting

This guide provides instructions for fitting variance component models using CVC.jl.

### Genotype Data

CVC.jl requires genotype data in **partitioned PLINK format**. Each genetic component must be a separate set of PLINK files:

```
analysis_dir/
├── G1.bed       # First genetic component
├── G1.bim
├── G1.fam
├── G2.bed       # Second genetic component
├── G2.bim
├── G2.fam
└── ...
```

**Important requirements:**
- Files must be named `G1.bed`, `G2.bed`, ..., `GK.bed` (and corresponding `.bim`/`.fam`)
- All `.fam` files must have identical sample ordering
- Samples in `.fam` must match the order of phenotype/covariate data
- Each SNP should appear in one and only one genotype partition
- SNPs should be quality-controlled before partitioning

### Phenotype Data

For time-to-event phenotypes:
- `u`: Vector of observed times (numeric, e.g., age at event or censoring)
- `delta`: Vector of censoring indicators (`true` = event observed, `false` = censored)

```julia
using DelimitedFiles

# Load time-to-event data
u = vec(readdlm("u.txt"))  # e.g., [42.3, 51.7, 38.9, ...]
delta = Bool.(vec(readdlm("delta.txt")))  # e.g., [true, false, true, ...]
```

### Covariate Matrix

The covariate matrix should include:
- Intercept (column of ones) if desired
- Age, sex, principal components, etc.

## Basic Usage

The tutorial examples assume the example data are located in the repository's data folder.

### Initialize and Fit Model

```@example cvc_tutorial
using CVC
using DelimitedFiles
using Random
using Statistics

Random.seed!(2025)

# Paths
data_dir = normpath(joinpath(dirname(pathof(CVC)), "..", "data"))
genotype_dir = data_dir  # Directory with G1.bed, G2.bed, etc.
temp_dir = tempdir()  # Or specify custom temporary directory

# Load example data
u = vec(readdlm(joinpath(data_dir, "u.txt")))
delta = Bool.(vec(readdlm(joinpath(data_dir, "delta.txt"))))
w = readdlm(joinpath(data_dir, "w.txt"))

# Initialize CVC model (use short settings for reproducible docs output)
cvcm = cvc(
    u,    # N-vector of observed times
    delta,          # N-vector of censoring indicators
    w,                 # N × p covariate matrix
    genotype_dir,      # Path to partitioned genotypes
    temp_dir;          # Path for temporary files
    J = 100,           # Number of jackknife blocks 
    B = 10             # Number of random Gaussian vectors
)

# Fit variance components
fit_me!(cvcm)

# Extract results
println("Total heritability:")
println("  h2 = ", round(cvcm.h2[], digits=4), " ± ", round(cvcm.h2se[], digits=4))

println("\nComponent-specific heritability (first 3):")
for k in 1:min(3, length(cvcm.h2k))
    println("  Component $k: ", round(cvcm.h2k[k], digits=4), " ± ", round(cvcm.h2kse[k], digits=4))
end

# Variance components
println("\nGenetic variance components (first 3): ", round.(cvcm.ϕg[1:min(3, end)], digits=4))
println("Environmental variance: ", round(cvcm.ϕe[], digits=4))

# Check data
println("\nSample size: ", length(u))
println("Censoring rate: ", round(mean(.!delta), digits=4))
println("Number of covariates: ", size(w, 2))

# NNLS estimates
println("\n=== NNLS (Non-Negative) Estimates ===")
println("h2 (NNLS) = ", round(cvcm.h2_nn[], digits=4), " ± ", round(cvcm.h2se_nn[], digits=4))
```

## Model Parameters

### Constructor Parameters

```julia
cvc(u, delta, w, datapath, tempdir_parent; J=100, B=10)
```

**Required arguments:**
- `u::Vector{Float64}`: Observed event times
- `delta::Vector{Bool}`: Censoring indicators (true = event observed, false = censored)
- `w::Matrix{Float64}`: N × p covariate matrix (include intercept if desired)
- `datapath::String`: Directory containing partitioned PLINK files (G1.bed, G2.bed, ...)
- `tempdir_parent::String`: Parent directory for temporary files

**Keyword arguments:**
- `J::Int`: Number of jackknife blocks for standard errors (default: 100)
- `B::Int`: Number of random projection vectors for efficiency (default: 10)

## Output Interpretation

### Heritability Estimates

For time-to-event models, results are on the **observed scale**:

- `h2`: Total narrow-sense heritability
- `h2k`: Component-specific heritabilities (K genetic + 1 environmental)
- `ϕg`: Genetic variance components (length K vector)
- `ϕe`: Environmental variance
- All estimates have corresponding `*se` fields for standard errors

Interpretation: `h2` represents the proportion of variance in the survival outcome explained by genetics.

### Non-Negative Estimates

All estimates have both unconstrained and constrained versions:
- Regular: `h2`, `h2k`, `ϕg`, etc.
- Constrained: `h2_nn`, `h2k_nn`, `ϕg_nn`, etc.

Constrained estimates impose variance components to be non-negative.

### Accessing Internal Values

The `cvc` object contains many internal fields:

```julia
# Sample information
cvcm.nobs        # Number of observations
cvcm.ncens       # Number censored
cvcm.propcens    # Proportion censored

# Dimensions
cvcm.ngpred      # Number of SNPs across all components
cvcm.nfpred      # Number of fixed effect predictors
cvcm.ngcomp      # Number of genetic components (K)

# Kaplan-Meier estimator
cvcm.km          # KaplanMeier object
```