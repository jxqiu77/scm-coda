# On eigenvalues of sample covariance matrices based on high-dimensional compositional data

This repository contains the simulation code for the paper "[On eigenvalues of sample covariance matrices based on high-dimensional compositional data](https://jxqiu77.github.io/assets/files/papers/2025-arXiv-CoDA-v4.pdf)" (arXiv: 2312.14420) by Qianqian Jiang, Jiaxin Qiu, and Zeng Li (2023).

The paper studies the asymptotic spectral behavior of sample covariance matrices constructed from high-dimensional compositional data, including:

- the limiting spectral distribution (LSD),
- the extreme eigenvalues,
- the central limit theorem (CLT) for linear spectral statistics (LSS),
- and a covariance-structure test for the basis data.

---

## Repository structure

```text
scm-coda/
├── Project.toml             # Julia project dependencies
├── Manifest.toml            # Julia environment lockfile
├── codes/
│   ├── CoDA_utils.jl        # Shared Julia utilities for covariance testing
│   ├── CovTest.jl           # Empirical size/power simulations for covariance testing
│   ├── CovTest_heatmap.jl   # Power heatmap simulation for the supplementary material
│   ├── CovTest_norm.jl      # Supplementary Frobenius-norm comparison figures
│   ├── install_R_packages.R # R script to install required R packages
│   ├── LSD.R                # Simulation for limiting spectral distribution
│   ├── LSS.R                # Simulation for CLT of linear spectral statistics
│   ├── LSS_plot.R           # Plotting script for normalized LSS histograms
│   └── Mp.R                 # Supplementary simulation for M_p(z)
├── data/
│   ├── CovTest/             # CSV / XLSX / LaTeX outputs for covariance testing
│   ├── LSS/                 # Simulation outputs for LSS
│   └── Stieltjes/           # Simulation outputs for M_p(z)
└── figure/
    ├── CovTest/             # Power heatmap and related figures
    ├── LSD/                 # LSD figures
    └── LSS/                 # LSS histogram figures
```


---

## Quick start

All commands below assume you are in the repository root.

```bash
# ============================
# 1 Dependencies
# ============================
julia --project=. -e 'using Pkg; Pkg.instantiate()' # Install Julia dependencies
Rscript codes/install_R_packages.R # Install R dependencies

# ============================
# 2 Main Simulation
# ============================
Rscript codes/LSD.R # Figures 1 - 2
Rscript codes/LSS.R # Tables 1 - 2
Rscript codes/LSS_plot.R # Figures 3 - 4
julia --project=. codes/CovTest.jl # Tables 3 - 4

# ============================
# 3 Supplementary simulations
# ============================
Rscript codes/Mp.R # Tables S1 - S2
```
