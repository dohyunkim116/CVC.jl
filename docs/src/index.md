# CVC.jl Introduction

**CVC.jl** (Censored Variance Component) is a Julia package for estimating heritability and variance components from genetic data with censored time-to-event phenotypes.

## Overview

Heritability estimation is a fundamental problem in statistical genetics, but traditional methods assume continuous, normally distributed phenotypes. Many important clinical outcomes are time-to-event (e.g., age at diagnosis, survival time). CVC.jl addresses this challenge by:

1. **Handling censored data**: Uses synthetic variable methods (Leurgans and Koul-Susarla-Van Ryzin) to transform right-censored time-to-event data into pseudo-continuous responses suitable for variance component analysis.

2. **Partitioned heritability**: Estimates variance components across K genetic components (e.g., chromosomes and gene sets) to understand the genetic architecture of complex traits.

3. **Scalability**: Handles biobank-scale data through memory-mapping, randomized trace estimator, and column-wise block jackknife for standard error estimation.

## What CVC Does

For survival/time-to-event phenotypes with right-censoring:

```math
h^2 = \frac{\sum_{k=1}^K \sigma_k^2}{\sum_{k=1}^K \sigma_k^2 + \sigma_e^2}
```

where $\sigma_k^2$ represents the variance component of the $k$th genetic partition and $\sigma_e^2$ is the environmental variance.

**Synthetic Variable Approach**: CVC transforms censored data $(T_i, \delta_i)$ where $T_i$ is observed time and $\delta_i$ indicates censoring status into synthetic continuous variables using Kaplan-Meier estimator of the censoring distribution.

## Key Features

- **Time-to-event phenotypes** with right-censoring support
- **Partitioned analysis**: Estimate $K$ genetic components + 1 environmental component
- **Standard errors**: Block jackknife with customizable number of blocks ($J=100$ default)
- **Computational efficiency**: Randomized trace estimator to reduce computation ($B = 10$ random Gaussian vectors default)
- **Flexible constraints**: Both unconstrained and constrained non-negative estimates
- **Memory efficiency**: Memory-mapped PLINK and working arrays for analyzing large-scale genetic data
