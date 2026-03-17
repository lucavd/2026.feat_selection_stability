# The Stability Illusion

**A simulation-based assessment of feature selection methods for high-dimensional metabolomics data**

## Overview

This project implements a comprehensive benchmarking framework for evaluating the **stability** of 12 feature selection methods applied to metabolomics data. The study is entirely *in silico*, using publicly available datasets for empirical parameter extraction and controlled simulations for systematic evaluation.

## Motivation

A recent meta-analysis of 244 clinical metabolomics studies found that 72% of 2,206 metabolites reported as statistically significant appeared in only one study, with 85% of these singletons being statistical noise (PMC11999569). This project systematically investigates *why* feature selection in metabolomics is unstable and *which methods* offer the best stability-accuracy trade-off.

## Methods Compared (12)

| Category | Methods |
|----------|---------|
| Filter | Wilcoxon + FDR, Fold-change, Volcano plot |
| Embedded | LASSO, Elastic Net, Knockoff filter |
| Wrapper | Boruta, Random Forest importance |
| Meta | Stability selection |
| Bayesian | Horseshoe prior, Spike-and-slab |
| ML | SHAP (XGBoost) |

## Simulation Scenarios (7 categories)

1. **p/n ratio** (5 to 50)
2. **Effect size** (FC 1.2 to 3.0)
3. **Multicollinearity** (r = 0.3 to 0.9)
4. **Non-linearity** (interaction terms)
5. **Confounders** (age, BMI)
6. **Missing data** (0% to 30% MNAR)
7. **Preprocessing** (log, pareto, auto-scaling, PQN)

## Key Metrics

- **Stability:** Nogueira index (with CI), Jaccard
- **Accuracy:** TPR, FDR (known ground truth in simulation)
- **Prediction:** AUC on held-out test set
- **Parsimony:** Number of features selected

## Project Structure

```
config/config.yaml          — Central configuration
R/00-08_*.R                 — Analysis scripts (sequential pipeline)
R/utils/                    — Helper functions
docs/01-06_*.md             — Research documentation
data/                       — Raw, processed, simulated data
results/                    — Feature selection results, metrics, figures
```

## Data Sources

- **MetaboLights** (EMBL-EBI): via `metabolighteR` R package
- **Metabolomics Workbench** (NIH): via `metabolomicsWorkbenchR` Bioconductor package

## Requirements

- R ≥ 4.3
- See `R/00_install_packages.R` for complete dependency list
- 8+ CPU cores recommended (Script 05 is computationally intensive)

## Running

```bash
# 1. Install packages
Rscript R/00_install_packages.R

# 2. Download public datasets
Rscript R/01_download_data.R

# 3. Run pipeline sequentially
Rscript R/02_preprocess.R
Rscript R/03_extract_empirical_params.R
Rscript R/04_simulate.R
Rscript R/05_feature_selection.R  # WARNING: 3-12 days
Rscript R/06_metrics.R
Rscript R/07_cross_validation.R   # Can run in parallel with 04-06
Rscript R/08_figures.R
```

## Documentation

Detailed research notes and design decisions are in `docs/`:
- [01_datasets_inventory.md](docs/01_datasets_inventory.md)
- [02_methods_inventory.md](docs/02_methods_inventory.md)
- [03_metrics_inventory.md](docs/03_metrics_inventory.md)
- [04_literature_review.md](docs/04_literature_review.md)
- [05_simulation_design.md](docs/05_simulation_design.md)
- [06_computational_plan.md](docs/06_computational_plan.md)
- [07_implementation_status.md](docs/07_implementation_status.md)
- [08_server_deployment.md](docs/08_server_deployment.md)

## License

See [LICENSE.md](LICENSE.md)
