# Figure Legends

## Figure 1 — Stability overview

Distribution of the Nogueira stability index across all 1590 simulated datasets for each of 11 feature selection methods. Methods are ordered by median stability. The dashed line marks the chance level (0.5). Box plots show the median (central line), interquartile range (box), and 1.5 IQR whiskers. Colors indicate method category: filter (orange), embedded (blue), wrapper (green), meta (yellow), Bayesian (dark blue), and ML-based (red).

## Figure 2 — Stability-accuracy trade-off

Scatter plot of mean stability (Nogueira index, x-axis) versus mean sensitivity (true positive rate, y-axis) for each method, averaged across all simulation scenarios. Point size is proportional to the average number of features selected. Dashed lines indicate the 0.5 reference for both axes. Methods in the lower-left quadrant exhibit both low stability and low sensitivity.

## Figure 3 — Scenario-specific stability heatmap

Heatmap of mean Nogueira stability index for each combination of method (columns) and scenario level (rows), faceted by simulation scenario (S1--S7 MVN, S8 semi-synthetic). Color scale ranges from 0 (dark) to 1 (bright). Empty tiles indicate insufficient convergence for that method-scenario combination.

## Figure 4 — Effect of dimensionality (p/n ratio)

Distribution of Nogueira stability index as a function of the p/n ratio (scenario S1), faceted by method. Increasing p/n ratio reflects higher dimensionality relative to sample size. Box plots show the distribution across 30 replications per level. Colors indicate method category.

## Figure 5 — Effect of feature correlation

Distribution of Nogueira stability index as a function of the inter-feature correlation level (scenario S3), faceted by method. Each panel shows one method with box plots across correlation levels. Colors indicate method category. Higher correlation among features is expected to reduce the identifiability of individual predictors and thus decrease selection stability.

## Figure 6 — FDR control

Distribution of the empirical false discovery rate (FDR) across all simulated datasets for each method. Methods are ordered by median FDR. The y-axis shows FDR as a percentage. The dashed red line marks the nominal 5% FDR threshold. Methods below this line adequately control the false discovery rate; methods above it exhibit liberal selection behavior.

## Figure 7 — Stability on real metabolomics data

Distribution of the Nogueira stability index computed from 100 bootstrap resamplings on each of 9 real metabolomics datasets, for each of the 11 methods. This panel validates the simulation results by showing that the stability patterns observed in silico are reproduced on real data. Box plots summarize the 9 per-dataset stability estimates.

## Figure 8 — Cross-database concordance

Spearman correlation of feature selection frequencies between pairs of datasets sharing the same analytical platform (GC-MS or LC-MS), for each method. Higher values indicate that the same features tend to be selected across independent datasets. This metric assesses the external reproducibility of the feature selection beyond within-dataset bootstrap stability.

## Supplementary Figure S1 — Composite: stability and FDR

Composite panel combining Figure 1 (stability overview, panel A) and Figure 6 (FDR control, panel B) for side-by-side comparison. Methods that appear stable may still exhibit uncontrolled false discovery rates, and vice versa.
