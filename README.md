# Unveiling the Drivers of Annual Rhythms in Gut Microbiota of Giant Pandas: A Coupling of Microbial Interaction Networks and Environmental Factors

## Overview

This repository contains data and code for analyzing seasonal dynamics of the giant panda gut microbiome across a full annual cycle. We combine Legendre polynomial basis expansion, ordinary differential equation (ODE) modeling, and adaptive sparse regression (ADSIHT) to infer interaction networks among 19 bacterial genera and five environmental factors. The model trained on Year 1 (12 months) is externally validated on an independent Year 2 dataset (4 months).

## Repository structure

```
├── README.md
├── data_raw/                          # Raw data
│   ├── asv_abundance_year1.xls        # Year 1 ASV table (1933 ASVs, 83 samples, 12 months)
│   ├── asv_abundance_year2.tsv        # Year 2 ASV table (2755 ASVs, 28 samples, 4 months)
│   ├── metadata_year1.csv             # Year 1 metadata (temperature, bamboo nutrition)
│   ├── metadata_year2.csv             # Year 2 metadata
│   └── sample_groups.txt              # Sample-to-month grouping
├── data_processed/                    # Preprocessed data (output of joint_preprocess.py)
│   ├── all_month.csv                  # Year 1: 19 genera + 5 env × 12 months
│   ├── all_month_year2.csv            # Year 2: 19 genera + 5 env × 4 months
│   ├── bacteria_month.csv             # Year 1 bacteria-only matrix
│   ├── bacteria_month_year2.csv       # Year 2 bacteria-only matrix
│   └── metadata_year2_clean.csv
├── scripts/                           # Analysis scripts
│   ├── joint_preprocess.py            # Joint preprocessing (intersection mode)
│   ├── preprocess_year2.py            # Year 2 standalone preprocessing
│   ├── run_model.R                    # Main ODE modeling pipeline
│   ├── validate.R                     # Year 2 external validation
│   ├── plot_validation.R              # Validation visualization
│   ├── analysis.R                     # Year 1 standalone analysis (optional)
│   ├── jtk_analysis.R                 # JTK_CYCLE rhythm detection
│   └── out-in.R                       # Network in/out degree analysis
└── results/                           # Model outputs and figures
    ├── interaction_res.Rdata          # R model object (BETA_est, basis, etc.)
    ├── edge_matrix.csv                # Interaction network edges
    ├── edge_matrix_bar.pdf            # In/out degree bar plot
    ├── effect_curves.csv              # Effect decomposition curves
    ├── effect_plot.pdf                # Effect decomposition figure
    ├── PC1_dynamics.pdf               # PC1 trajectories
    ├── pc1_contributions.csv          # PC1 contribution quantification
    ├── net.pdf                        # Network visualization
    ├── panda_net.cys                  # Cytoscape network file
    ├── bardata.csv                    # Bar plot data
    ├── Figure 1.pdf / Figure 2.pdf    # Main figures
    ├── jtk_input.csv / jtk_results.csv / jtk_summary.csv  # JTK_CYCLE results
    ├── JTKresult_jtk_input.csv        # MetaCycle raw output
    ├── year2_predicted.csv / year2_observed.csv  # Validation predictions
    ├── year2_genus_metrics.csv        # Per-genus validation metrics
    └── validation_*.pdf               # Validation figures (4 panels)
```

## Data

### Year 1 (Training set, Sep 2019 – Aug 2020)

- 1,933 ASVs across 83 fecal samples from 8 giant pandas (DD, FX, GG, JN, MD, ME, ML, MM)
- 12 monthly time points
- Environmental variables: temperature, crude fat%, crude protein%, crude fiber%, carbohydrate%

### Year 2 (Validation set)

- 2,755 ASVs across 28 fecal samples from 4 pandas (CP, DP, QP, XP; same individuals as Year 1 with different identifiers)
- 4 time points: CP = March (spring), XP = August (summer), QP = October (autumn), DP = December (winter)
- Same environmental variables as Year 1

## Analysis pipeline

### 1. Joint preprocessing (`joint_preprocess.py`)

- Aggregate ASVs to genus-level relative abundance
- Filter each year independently: mean relative abundance > 0.05%, prevalence > 5%
- Take intersection of retained genera across both years → **19 shared genera**
- Average by month within each year
- Output: `data_processed/all_month.csv` (24 variables × 12 months) and `all_month_year2.csv` (24 variables × 4 months)

### 2. Dynamic system modeling (`run_model.R`)

**Legendre polynomial basis expansion:**
- Smooth each variable's time series with `smooth.spline(spar = 0.78)`
- Expand into Legendre polynomials (order = 3) via ODE path integration:

  $$\Phi_j(t) = \int_{1}^{t} P_j(x(s))\,ds$$

- Construct block-diagonal design matrix (one block per bacterium equation)
- Remove non-self leg0 columns to avoid redundancy

**Sparse regression:**
- ADSIHT (Adaptive Double Sparse Iterative Hard Thresholding) with group structure:
  - Each variable's 4 basis coefficients form a group
  - Parameters: `kappa = 0.99`, `ic.scale = 0.6`, `ic.coef = 0.6`
  - Lambda selected by minimizing information criterion (IC)
- Equations without self-effect (leg0 + leg1 both zero) are supplemented with Ridge + BFGS optimization (`alpha = 5e-3`)

**Coefficient rescaling:**
- Transform coefficients back to original scale via `(sY/sX)` ratio

### 3. Effect decomposition

For each genus *j*, the predicted trajectory is decomposed into three components:

$$\hat{y}_j(t) = \text{Self}_j(t) + \sum_{k \neq j,\ k \in \text{bacteria}} \text{Eff}_k(t) + \sum_{k \in \text{env}} \text{Eff}_k(t)$$

- Effect curves saved to `results/effect_curves.csv`

### 4. Interaction network construction

From derivative effects: for each genus *j*, compute the mean derivative effect of each variable *k* on *j*'s rate of change. Build a directed network: `source = k → target = j`, edges weighted by |mean_effect|. Self-loops excluded. Output: `results/edge_matrix.csv`.

### 5. PC1 beta-diversity analysis

Construct three effect matrices (total, bacterial interaction, environmental), compute Bray-Curtis distances, PCoA, extract PC1. Quantify interaction vs. environment contributions via three methods: variance ratio, R<sup>2</sup>, and standardized regression. Output: `results/pc1_contributions.csv`.

### 6. JTK_CYCLE rhythm analysis (`jtk_analysis.R`)

Apply JTK_CYCLE (non-parametric rhythm detection, period = 12 months) to each of the four effect curves (Total, Independent, Bacterial, External) for all 19 genera. Identify rhythmically significant curves (Bonferroni-adjusted P < 0.05) and classify each genus's rhythmic driver: self-driven, interaction-driven, or environment-driven.

### 7. External validation (`validate.R`)

Year 2 has only 4 time points. The validation procedure:

1. Interpolate each Year 2 variable to 12 months via `smooth.spline`
2. Project onto the same Legendre basis as Year 1
3. Multiply by Year 1's `BETA_est` to obtain 12-month predicted trajectories
4. Extract predictions at months 3, 8, 10, 12 and compare with observations

**Why interpolate to 12 months?** Legendre basis functions are computed via ODE path integration from *t* = 1 to 12. The same integration path must be used for prediction to align with the training basis space.

## Key results

### Model performance

| Metric | Value |
|--------|-------|
| Parameters | n_order = 3, spar = 0.78 |
| Year 1 model residual | 0.099 |
| Year 2 overall RMSE | 0.059 |
| Year 2 overall R<sup>2</sup> | **0.802** |
| Baseline R<sup>2</sup> (Year 1 direct vs. Year 2) | 0.791 |
| Genera with R<sup>2</sup> > 0 | 2 (*Streptococcus* 0.47, *Actinobacillus* 0.05) |

The model R<sup>2</sup> (0.802) exceeds the baseline (0.791), demonstrating that the interaction network learned from Year 1 is partially transferable to Year 2.

### Interaction network

- **59 directed edges** among 19 bacterial genera and 5 environmental factors
- Dominant hub genera: *Escherichia-Shigella*, *Streptococcus*, *Clostridium*, *Turicibacter* (>80% cumulative abundance)

### PC1 contribution

| Method | Bacterial interaction | Environment |
|--------|----------------------|-------------|
| Variance ratio | 28.9% | 6.6% |
| R<sup>2</sup> percentage | 31.9% | 68.1% |
| Regression coefficient | 0.923 | 0.093 |

### JTK_CYCLE rhythm

| Component | Significant curves | Key genera |
|-----------|-------------------|------------|
| Total | 9/19 | *Streptococcus*, *Terrisporobacter*, *Lactococcus*, *Escherichia-Shigella* |
| Independent | 11/19 | *Streptococcus*, *Enterobacter*, *Terrisporobacter* |
| Bacterial | 11/19 | *Lactococcus*, *Acinetobacter*, *Ligilactobacillus* |
| External | 5/19 | *Terrisporobacter*, *Weissella*, *Hafnia-Obesumbacterium* |

*Lactococcus* and *Acinetobacter* show the clearest **interaction-driven** rhythm pattern: significant total periodicity, non-significant self-effect, but highly significant bacterial interaction effect (P < 0.001).

## Dependencies

### R packages

```
vegan, ggplot2, dplyr, tidyr, ADSIHT, deSolve, Matrix, 
orthopolynom, reshape2, patchwork, purrr, MetaCycle
```

### Python packages

```
pandas, numpy, re
```

## Reproducing the analysis

1. **Preprocess data:**
   ```bash
   python scripts/joint_preprocess.py
   ```

2. **Run Year 1 model and generate all outputs:**
   ```r
   source("scripts/run_model.R")
   ```

3. **Run JTK_CYCLE rhythm analysis:**
   ```r
   source("scripts/jtk_analysis.R")
   ```

4. **Run Year 2 external validation:**
   ```r
   source("scripts/validate.R")
   ```

5. **Generate validation figures:**
   ```r
   source("scripts/plot_validation.R")
   ```

6. **Network degree analysis (optional):**
   ```r
   source("scripts/out-in.R")
   ```

> **Note:** Scripts contain hard-coded working directory paths (`setwd("D:/data/panda")`). Adjust to your local path before running. The `run_model.R` script reads from `data_processed/all_month.csv` and writes to `results/`.

## Citation

If you use this code or data, please cite the corresponding paper (under review).
