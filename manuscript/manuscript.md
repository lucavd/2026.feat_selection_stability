---
title: "The stability illusion: A simulation-based assessment of feature selection methods for high-dimensional metabolomics data"
author:
  - name: "[Author name]"
    affiliation: "[Affiliation]"
date: "2026"
bibliography: references.bib
csl: csl/numeric.csl
link-citations: true
number-sections: true
reference-section-title: References
abstract: |
  **Motivation:** Feature selection is routinely applied to metabolomics data to identify discriminant metabolites, yet the stability of the resulting feature sets --- i.e., their reproducibility under data perturbation --- is rarely assessed. Unstable selections undermine biological interpretation and have contributed to the poor reproducibility of omics biomarker studies.

  **Results:** We benchmarked 11 feature selection methods spanning six methodological categories (filter, embedded, wrapper, meta-learning, Bayesian, and machine-learning-based) across 1,590 simulated datasets covering eight scenarios that varied dimensionality (*p/n* ratio 5--50), effect size (fold change 1.2--3.0), inter-feature correlation, missing data rate, and preprocessing pipeline. We quantified stability using the chance-corrected Nogueira index. Across all methods and scenarios, the median Nogueira index was 0.19, and no method achieved a median above 0.45. Methods that attained higher sensitivity (e.g., fold-change filtering, TPR = 0.50) did so at the cost of catastrophic false discovery rates (>80%). Validation on nine public metabolomics datasets (NMR, LC-MS, GC-MS; *n* = 80--1,005; *p* = 50--1,533) confirmed these findings, with real-data stability values of 0.18--0.57 (method medians). Cross-database concordance of selection frequencies between independent datasets sharing the same analytical platform was centered on zero (mean Spearman ρ = −0.02; 64% of dataset pairs shared no selected features), indicating that feature rankings do not generalize across studies.

  **Conclusions:** Feature selection in metabolomics produces unstable, non-reproducible results regardless of the method employed. When prediction is the goal, full-dimensional penalized models should be preferred. When biological interpretation is required, dimension reduction with pathway-level analysis offers more robust insights than single-metabolite selection. If feature selection is unavoidable, stability must be quantified and reported alongside accuracy metrics, and multi-method consensus with independent validation should be mandatory.

  **Keywords:** feature selection, stability, metabolomics, biomarker discovery, reproducibility, simulation, Nogueira index
---

# Introduction

Metabolomics generates high-dimensional datasets in which the number of measured features (*p*) frequently exceeds the number of biological samples (*n*), creating a statistical environment where overfitting and spurious associations are the norm rather than the exception [@broadhurst2006statistical]. Feature selection (FS) --- the process of identifying a subset of variables that best discriminate between experimental conditions --- is widely regarded as an essential step in the metabolomics analysis pipeline, from untargeted biomarker discovery to the construction of clinical diagnostic panels [@xia2009metaboanalyst; @alseekh2021mass].

The methodological landscape of FS is vast. @guyon2003introduction and @saeys2007review established the canonical taxonomy distinguishing filter, wrapper, and embedded methods, each offering different trade-offs between computational cost, model dependency, and theoretical guarantees. More recent approaches include stability-enhanced penalized regression [@meinshausen2010stability], knockoff-based false discovery rate control [@barber2015controlling; @candes2018panning], Bayesian variable selection with sparsity-inducing priors [@mitchell1988bayesian; @george1993variable], and model-agnostic importance measures such as SHAP values [@lundberg2017unified]. Despite this methodological richness, a fundamental question is rarely addressed in practice: *are the selected features reproducible?*

The stability of FS --- the extent to which the same features are selected under perturbation of the training data --- has been recognized as a critical issue in bioinformatics for over a decade. @lustgarten2009measuring first formalized stability measurement in biomedical contexts, demonstrating that common FS methods produce highly variable feature sets. @boulesteix2009stability showed that ranked gene lists from microarray experiments exhibit poor overlap across bootstrap samples, and @haury2011influence benchmarked 32 methods on gene expression data, finding stability values between 0.1 and 0.5 for most approaches. @he2010stable characterized the stability--accuracy trade-off, arguing that high prediction accuracy does not guarantee stable feature identification.

The problem is especially acute in metabolomics. Three structural properties of metabolomics data conspire against stable feature selection. First, the typical *p/n* ratio in untargeted studies (5--50) creates an inherently ill-posed selection problem where many distinct feature subsets can yield equivalent predictive performance --- the "Rashomon effect" applied to feature space. Second, metabolites within the same biochemical pathway are highly correlated [@dunn2011procedures], so multiple features carry nearly identical discriminative information; which one is selected becomes a function of sampling noise. Third, the analytical pipeline itself introduces variability: preprocessing choices (normalization, transformation, scaling) alter the relative importance of features [@vandenberg2006centering], amplifying the "researcher degrees of freedom" that @simmons2011false identified as a major driver of irreproducible findings.

These concerns are not merely theoretical. The history of metabolomics biomarker research is marked by poor external validation and retracted claims [@ransohoff2004rules; @ioannidis2005why]. @broadhurst2006statistical warned that small sample sizes and multiple testing render most metabolomics discoveries unreliable, a critique that remains unaddressed in current practice. More recently, @hedou2024discovery demonstrated that standard FS methods fail to produce reliable biomarker signatures across omics datasets and proposed Stabl, a stability-enhanced LASSO framework, as a remedy. Their work in *Nature Biotechnology* underscores the growing recognition that FS instability is not a minor nuisance but a fundamental barrier to translational success.

Despite this recognition, no comprehensive benchmark of FS stability exists specifically for metabolomics data. The seminal benchmarks by @haury2011influence and @bommert2020benchmark focused on gene expression and general classification datasets, respectively, and may not capture the specific challenges of metabolomics (e.g., MNAR missingness, log-normal distributions, pathway-driven correlation structures). Furthermore, most stability assessments rely on within-study resampling, which may overestimate the generalizability of selected features to independent cohorts.

In this study, we provide the first systematic evaluation of FS stability tailored to metabolomics. We benchmark 11 methods from six categories --- Wilcoxon rank-sum with FDR correction (filter), fold-change filtering (filter), volcano plot (filter), LASSO (embedded), elastic net (embedded), model-X knockoff filter (embedded), Boruta (wrapper), random forest importance (wrapper), stability selection (meta), spike-and-slab regression (Bayesian), and SHAP-based XGBoost importance (ML) --- across 1,590 simulated datasets spanning eight experimental scenarios. We complement the simulation with cross-validation on nine public metabolomics datasets and, critically, with cross-database concordance analysis that tests whether feature selections generalize across independent studies of the same platform. Our primary stability metric is the chance-corrected Nogueira index [@nogueira2018stability], which accounts for the expected agreement between random selections and provides jackknife-based confidence intervals.

Our results reveal a pervasive "stability illusion": feature selections that appear reproducible within a single study (via bootstrap resampling) fail to generalize across datasets. No method achieves consistently high stability, and the stability--accuracy trade-off is unfavorable at realistic effect sizes. We conclude with practical recommendations, arguing that --- in alignment with @harrell2015regression --- the field should generally avoid single-feature selection in favor of dimension reduction and pathway-level interpretation.


# Methods

## Datasets

### Public metabolomics datasets

We assembled nine publicly available binary-classification metabolomics datasets spanning three analytical platforms (Table 1). Seven datasets were obtained from the CIMCB benchmark collection [@mendez2019comparative], which provides pre-processed Excel files with standardized class labels. One additional dataset was downloaded from MetaboLights [@haug2020metabolights] (MTBLS28, NSCLC vs. control, urine LC-MS, *n* = 1,005, *p* = 1,359) and one from Metabolomics Workbench [@sud2016metabolomics] (ST001706, renal cell carcinoma vs. control, urine NMR, *n* = 256, *p* = 50).

Inclusion criteria required a binary case--control design with at least 20 samples per group and at least 50 measured features. Features with >50% missing values were removed, followed by samples with >40% missing values. Remaining missing values were imputed with the feature-wise median. All datasets were log~2~-transformed and auto-scaled (zero mean, unit variance) prior to analysis unless otherwise specified by the simulation scenario.

**Table 1.** Public metabolomics datasets used for cross-validation and concordance analysis.

| Accession | Platform | Matrix | Comparison | *n* | *p* | *p/n* |
|-----------|----------|--------|------------|----:|----:|------:|
| ST001047 | NMR | Urine | Gastric cancer vs. control | 83 | 149 | 1.80 |
| MTBLS404 | LC-MS | Urine | Male vs. female (Sacurine) | 184 | 120 | 0.65 |
| ST001000 | LC-MS | Stool | UC vs. CD (IBD) | 107 | 1,533 | 14.33 |
| MTBLS136 | LC-MS | Serum | E-only vs. E+P (hormone use) | 668 | 787 | 1.18 |
| MTBLS92 | LC-MS | Plasma | Breast cancer pre- vs. post-chemo | 253 | 138 | 0.55 |
| ST000369 | GC-MS | Serum | Lung adenocarcinoma vs. control | 80 | 181 | 2.26 |
| ST000496 | GC-MS | Saliva | Periodontal pre- vs. post-treatment | 100 | 69 | 0.69 |
| MTBLS28 | LC-MS | Urine | NSCLC vs. control | 1,005 | 1,359 | 1.35 |
| ST001706 | NMR | Urine | Renal cell carcinoma vs. control | 256 | 50 | 0.20 |

### Empirical parameter extraction

For each dataset, we estimated the Ledoit-Wolf shrinkage covariance matrix [@ledoit2004well; @schafer2005shrinkage], feature-wise log~2~ fold changes, and the effective dimensionality (number of eigenvalues explaining 95% of variance). These empirical parameters served as input for realistic simulation design and as the basis for semi-synthetic spike-in experiments (scenario S8).

## Simulation framework

### Parametric simulation (scenarios S1--S7)

Data were generated from a multivariate normal (MVN) distribution using the empirical covariance structure of three platform-representative datasets (NMR: ST001047; LC-MS: ST001000; GC-MS: ST000369) as the base correlation matrix, scaled by a scenario-specific factor [@genz2009computation]. The MVN samples were exponentiated to produce log-normal distributions consistent with metabolomics concentration data [@vandenberg2006centering]. Ground truth was established by applying a multiplicative fold change to *p*~true~ = 20 randomly selected features in the case group.

Seven scenarios systematically varied one experimental factor while holding all others at base-level values (*n*~case~ = *n*~control~ = 50, *p* = 500, FC = 1.5, empirical correlation, 5% MNAR missingness, log + autoscaling):

- **S1** --- Dimensionality (*p/n* ratio): 5, 10, 20, 50
- **S2** --- Effect size (fold change): 1.2, 1.5, 2.0, 3.0
- **S3** --- Inter-feature correlation scale: 0.3, 0.5, 0.7, 0.9
- **S4** --- Feature interactions: 0, 5, 10 pairwise interactions
- **S5** --- Confounding variables: 0, 1, 3 confounders
- **S6** --- Missing data rate (MNAR): 0%, 5%, 15%, 30%
- **S7** --- Preprocessing: none, log + autoscale, log + Pareto, PQN + autoscale

Each scenario--level combination was replicated 30 times with independent random seeds derived from the project master seed (42), yielding 780 MVN-simulated datasets.

### Semi-synthetic simulation (scenario S8)

To complement the parametric simulations, we implemented a semi-synthetic spike-in design [@sankaran2025semisynthetic; @franklin2014plasmode]. For each of the nine real datasets, control samples were duplicated to create a synthetic case group, and a known fold change (1.2, 1.5, or 2.0) was applied to 20 randomly selected features. This preserves the original distributional properties, correlation structure, and platform-specific artifacts while maintaining known ground truth. Each spike-in was replicated 30 times, yielding 810 semi-synthetic datasets (total: 1,590 simulated datasets).

## Feature selection methods

We evaluated 11 FS methods spanning six methodological categories (Table 2). All methods were applied through a common interface that returned a binary selection vector. For methods that produce continuous importance scores (RF importance, SHAP), an elbow-point detection algorithm was applied to the sorted importance values, with a fallback to the top-20 features if the elbow did not converge.

**Table 2.** Feature selection methods evaluated in this study.

| Method | Category | Key parameter(s) | Reference |
|--------|----------|-------------------|-----------|
| Wilcoxon + FDR | Filter | BH-adjusted α = 0.05 | @wilcoxon1945individual; @benjamini1995controlling |
| Fold change | Filter | FC threshold = 1.5 | --- |
| Volcano plot | Filter | FC = 1.5 ∩ FDR = 0.05 | @li2012volcano |
| LASSO | Embedded | λ~1se~, 10-fold CV | @tibshirani1996regression; @friedman2010regularization |
| Elastic net | Embedded | α = 0.5, λ~1se~ | @zou2005regularization; @friedman2010regularization |
| Knockoff filter | Embedded | Target FDR = 0.10 | @barber2015controlling; @candes2018panning |
| Boruta | Wrapper | maxRuns = 100, *p* = 0.01 | @kursa2010boruta |
| RF importance | Wrapper | 1,000 trees, permutation | @breiman2001random; @wright2017ranger |
| Stability selection | Meta | π̂ = 0.75, PFER = 1 | @meinshausen2010stability; @hofner2015controlling |
| Spike-and-slab | Bayesian | bigp.smalln mode | @ishwaran2005spike; @ishwaran2010spikeslab |
| SHAP + XGBoost | ML | 200 rounds, η = 0.1 | @lundberg2017unified; @chen2016xgboost |

### Bootstrap stability assessment

For each simulated dataset, every applicable FS method was applied to *B* = 100 stratified bootstrap samples (63.2% subsample fraction). The resulting *B* binary selection vectors formed the input for stability computation. This bootstrap-and-evaluate protocol follows the framework established by @nogueira2018stability.

### Method applicability

Not all methods were applicable to all datasets. The knockoff filter requires *n* > *p* [@candes2018panning] and failed on all MVN scenarios (S1--S7) where *p* ≥ 500 and *n* ≤ 100, succeeding only on semi-synthetic datasets with favorable *p/n* ratios (90 of 810 S8 datasets; 11.1% success rate). Fold-change filtering was evaluated only under no-preprocessing conditions (S7/none; 30 datasets), as log-transformation alters fold-change interpretation. Volcano plot filtering required both fold-change and statistical significance criteria to be met simultaneously, resulting in convergence for 268 of the possible scenarios. These differential success rates are reported transparently and accounted for in all comparisons.

## Evaluation metrics

### Stability metrics

The primary stability metric was the Nogueira index [@nogueira2018stability], a chance-corrected measure defined as:

$$\hat{\Phi}(\mathbf{Z}) = 1 - \frac{\overline{p}_f(1 - \overline{p}_f)}{k/p \cdot (1 - k/p)} \cdot \frac{1}{B-1} \sum_{i=1}^{B}\sum_{j=1}^{p} \frac{(z_{ij} - \overline{p}_j)^2}{\overline{p}_f(1 - \overline{p}_f)}$$

where **Z** is the *B* × *p* binary selection matrix, $\overline{p}_j$ is the selection frequency of feature *j*, $\overline{p}_f$ is the mean selection frequency, and *k* is the average number of selected features. The index ranges from −1 to 1, with 0 indicating chance-level agreement and 1 indicating perfect stability. Jackknife variance estimates and 95% confidence intervals were computed following the procedure in @nogueira2018stability. The Jaccard similarity coefficient [@jaccard1901etude] was computed as a secondary metric.

An important methodological caveat must be noted: when a method consistently selects zero features across all bootstrap samples, the Nogueira index equals 1 (perfect stability of the empty set). This artifact occurs because selecting nothing is trivially reproducible. We report this phenomenon explicitly in the results and flag affected observations.

### Accuracy metrics

Because the simulation framework provides known ground truth, we computed the true positive rate (TPR = sensitivity) and false discovery rate (FDR) for each method--dataset combination. The Matthews correlation coefficient (MCC) was computed as a secondary balanced metric [@matthews1975comparison; @chicco2020advantages].

### Predictive performance

For each method--dataset combination, a logistic regression model was fit on the selected features using a 70/30 stratified train/test split, and the area under the ROC curve (AUC) was computed on the held-out test set [@robin2011proc; @fawcett2006introduction]. When no features were selected, AUC was recorded as NA.

### Parsimony

The mean number of features selected across bootstrap samples and the number of features stably selected (selection frequency ≥ 50%) were recorded.

## Cross-validation on real datasets

The bootstrap stability protocol was applied identically to the nine public datasets, with *B* = 100 bootstrap resamples per dataset per method. Because ground truth is unknown for real data, only stability, parsimony, and stably-selected feature counts are reported.

## Cross-database concordance

To assess external reproducibility, we computed the Spearman rank correlation and Jaccard similarity of selection frequency vectors between all pairs of datasets sharing the same analytical platform. For each dataset, the selection frequency vector records how often each feature was selected across the 100 bootstrap resamples. Features were matched by name across datasets, restricting the comparison to common features (range: 69--567 common features per pair). Seven cross-database pairs were evaluable (1 GC-MS, 6 LC-MS). This analysis directly tests whether the same features are consistently identified as important across independent studies, going beyond within-study bootstrap stability.

## Software and reproducibility

All analyses were conducted in R 4.5.3 [@rcoreteam2024r]. Key packages included `glmnet` [@friedman2010regularization] for LASSO/elastic net, `Boruta` [@kursa2010boruta], `ranger` [@wright2017ranger] for random forests, `stabs` [@hofner2015controlling] for stability selection, `knockoff` [@candes2018panning], `spikeslab` [@ishwaran2010spikeslab], `xgboost` [@chen2016xgboost] with `shapviz` [@lundberg2017unified] for SHAP, `stabm` [@bommert2021stabm] for stability indices, `mvtnorm` [@genz2009computation] and `corpcor` [@schafer2005shrinkage] for simulation, `pROC` [@robin2011proc] for AUC, and `ggplot2` [@wickham2016ggplot2] for visualization. Parallelization was implemented via `mclapply` with 40 workers. The total computational time for the feature selection pipeline (Script 05) was approximately 56 hours on a 128-core server. All code and configuration files are available at [repository URL].


# Results

## Overall stability is low across all methods

Across all 12,985 successful method--dataset evaluations (from 14,576 total; 10.9% failure rate, dominated by the knockoff filter), the median Nogueira stability index was 0.19 and the mean was 0.31 (Figure 1). No method achieved a median stability above 0.45 in the simulation. The best-performing method by median stability was the knockoff filter (median = 0.45), but this result is based on only 90 successful evaluations (all from semi-synthetic scenario S8) and should be interpreted with extreme caution given the method's 94.3% failure rate. Among methods evaluable across all scenarios, Boruta achieved the highest median stability (0.23), followed by RF importance (0.22) and SHAP-XGBoost (0.21).

Three filter-based methods (volcano, wilcoxon_fdr, stability_selection) exhibited strongly bimodal distributions: a large proportion of evaluations yielded Nogueira values near 0 (no features selected or unstable selections) while a substantial minority achieved values near 1.0 (Figure 1). Closer inspection revealed that 318 evaluations (2.4% of successes) achieved Nogueira = 1.0 by selecting approximately zero features across all bootstrap resamples. This "empty-set artifact" occurred predominantly at high *p/n* ratios (S1/pn50) and low effect sizes (S2/fc1.2), where methods correctly identified that no features exceeded their significance threshold --- a trivially stable but scientifically uninformative outcome. For volcano plot filtering, 108 of 268 successful evaluations (40.3%) exhibited this artifact.

## The stability--accuracy trade-off is unfavorable

Figure 2 reveals a striking pattern: all methods cluster in the lower-left quadrant, indicating simultaneously low stability and low sensitivity. The only method achieving TPR > 0.40 was fold-change filtering (TPR = 0.50), but at the cost of an empirical FDR of 84% (Figure 6), selecting on average 99 of 500 features. Fold-change filtering essentially detects every feature with an expression difference above the threshold, regardless of statistical significance, capturing true positives alongside a vastly larger number of false positives.

Among methods with controlled FDR (<5% empirical FDR), the best sensitivity was achieved by stability selection (TPR = 0.09, FDR = 0.4%), Wilcoxon + FDR (TPR = 0.20, FDR = 0.3%), and spike-and-slab (TPR = 0.16, FDR = 1.3%). However, these TPR values indicate that even the best-calibrated methods detect fewer than one in five true features at the base scenario (FC = 1.5, *p/n* = 5).

## Effect size is the dominant determinant of stability

The scenario-specific heatmap (Figure 3) reveals that effect size (S2) exerts the strongest influence on stability. At FC = 1.2 (realistic for many metabolomics studies), all methods showed mean stability below 0.16 (Table 3). At FC = 1.5, stability remained below 0.19 for most methods. Only at FC = 3.0 --- an effect size rarely observed in untargeted metabolomics outside of pharmacological interventions --- did several methods exceed 0.50 (Boruta: 0.79; Wilcoxon + FDR: 0.84; LASSO: 0.64).

**Table 3.** Mean Nogueira stability index and true positive rate by effect size (scenario S2) for each method. FDR values are empirical.

| Method | Stab. (1.2) | TPR (1.2) | Stab. (1.5) | TPR (1.5) | Stab. (2.0) | TPR (2.0) | Stab. (3.0) | TPR (3.0) |
|--------|:-----------:|:---------:|:-----------:|:---------:|:-----------:|:---------:|:-----------:|:---------:|
| Boruta | 0.11 | 0.01 | 0.18 | 0.05 | 0.41 | 0.37 | 0.79 | 0.91 |
| Elastic net | 0.05 | 0.00 | 0.14 | 0.06 | 0.42 | 0.74 | 0.54 | 0.98 |
| LASSO | 0.04 | 0.00 | 0.14 | 0.02 | 0.46 | 0.51 | 0.64 | 0.82 |
| RF importance | 0.10 | 0.00 | 0.16 | 0.03 | 0.33 | 0.11 | 0.48 | 0.19 |
| SHAP-XGB | 0.11 | 0.01 | 0.15 | 0.02 | 0.30 | 0.10 | 0.45 | 0.14 |
| Spike-slab | 0.03 | 0.00 | 0.08 | 0.00 | 0.22 | 0.09 | 0.36 | 0.75 |
| Stab. sel. | 0.02 | 0.00 | 0.08 | 0.00 | 0.27 | 0.02 | 0.54 | 0.20 |
| Wilcoxon + FDR | 0.15 | 0.00 | 0.05 | 0.00 | 0.31 | 0.06 | 0.84 | 0.86 |

The practical implication is clear: at effect sizes typical of metabolomics case-control studies (FC = 1.2--1.5), no method can reliably distinguish true features from noise. The signal-to-noise ratio is simply too low for any FS algorithm to achieve stable performance in a *p* = 500, *n* = 100 setting.

## Dimensionality exacerbates instability, with a paradoxical artifact

Increasing the *p/n* ratio from 5 to 50 (scenario S1) caused stability to decrease for most methods (Figure 4). Boruta dropped from 0.19 (pn5) to 0.07 (pn50); LASSO from 0.14 to 0.02; elastic net from 0.15 to 0.03. However, three methods --- stability selection, wilcoxon_fdr, and volcano --- showed an apparent increase in stability at pn50, reaching Nogueira = 1.0. This paradoxical result is entirely artifactual: at extreme dimensionality, these threshold-based methods consistently selected zero features across all bootstrap resamples, achieving perfect agreement on the empty set. This artifact highlights an important limitation of the Nogueira index (and indeed any stability metric): a stability value of 1.0 is meaningless when the selected set is empty.

## Correlation structure has modest effects

Contrary to intuition, increasing inter-feature correlation (scenario S3) did not uniformly decrease stability (Figure 5). For most methods, stability values were essentially flat across correlation levels (0.3--0.9), with differences of less than 0.05 between the highest and lowest correlation conditions. The exception was volcano plot filtering, which exhibited higher stability at high correlation (0.69 at cor_high vs. 0.21 at cor_low), likely because correlated features tend to share significance status, making the joint FC + p-value criterion more deterministic. The weak effect of correlation on stability suggests that, at realistic effect sizes (FC = 1.5), the dominant source of instability is signal weakness rather than the "Rashomon effect" from correlated features.

## Missing data, confounders, and interactions have minor effects

Scenarios S4 (interactions), S5 (confounders), and S6 (missing data) produced relatively modest effects on stability (Figure 3). Increasing the missing data rate from 0% to 30% reduced mean Boruta stability from 0.21 to 0.14, a decrease of 0.07. Adding 3 confounders (S5) or 10 interactions (S4) did not meaningfully alter stability for any method (mean changes < 0.03). Preprocessing pipeline choice (S7) also had limited impact, with the exception that unprocessed data slightly altered rankings without consistently improving or worsening stability.

## Semi-synthetic simulation yields more optimistic results

Scenario S8 (semi-synthetic spike-in) produced substantially higher stability values than the parametric MVN simulations (Figure 3). At FC = 2.0, several methods exceeded 0.60: Boruta (0.73), wilcoxon_fdr (0.69), stability selection (0.64), elastic net (0.63), and LASSO (0.62). Even at FC = 1.2, mean stability values ranged from 0.12 (spike-slab) to 0.42 (stability selection). The improvement relative to MVN scenarios reflects the fact that semi-synthetic data preserve the real correlation structure, which provides additional information that FS methods can exploit. This observation is consistent with the argument by @sankaran2025semisynthetic that semi-synthetic simulations produce more realistic benchmarks. However, it also cautions that MVN-based simulations may underestimate FS performance in settings where the data possess structured covariance, while semi-synthetic results may be optimistic for datasets with different correlation patterns than the source data.

## Predictive performance is independent of stability

AUC values were generally high and largely independent of stability (Figure 2). Among the nine methods with sufficient AUC data, mean AUC ranged from 0.84 (Boruta, RF importance, SHAP-XGBoost) to 0.94 (wilcoxon_fdr, spike-slab, stability selection, LASSO). The knockoff filter achieved mean AUC = 1.00, but based on only 76 evaluations with favorable *p/n* ratios. Critically, fold-change filtering --- the only method with AUC near chance (0.59) --- achieved the highest TPR, demonstrating that indiscriminate feature selection degrades predictive performance by including noise features. The decoupling of AUC from stability confirms that good prediction does not require identifying specific features; multiple distinct feature subsets can yield equivalent predictive accuracy.

## FDR control varies dramatically across methods

Empirical FDR control (Figure 6) divided the methods into three groups. The first group maintained FDR well below 5%: stability selection (0.4%), wilcoxon_fdr (0.3%), spike-and-slab (1.3%), and volcano (0%). The second group showed moderate inflation: Boruta (12.0%), LASSO (6.5%), RF importance (7.9%), SHAP-XGBoost (10.4%), and knockoff (11.5%). The third group was severely miscalibrated: elastic net (11.7% mean but with a long tail extending to 100%) and fold-change filtering (84.2%).

The knockoff filter's FDR inflation (11.5% vs. the target 10%) is noteworthy given its theoretical guarantee of exact FDR control. This discrepancy likely reflects the sensitivity of the model-X knockoff framework to violations of the assumed feature distribution [@candes2018panning], which are inevitable when knockoffs are constructed from estimated (rather than known) covariance structures.

## Validation on real metabolomics data

Bootstrap stability assessment on the nine public datasets confirmed the simulation findings (Figure 7). Method-level median stability values ranged from 0.00 (knockoff; converged on only 2 of 9 datasets) to 0.56 (fold_change). Among broadly applicable methods, wilcoxon_fdr (median = 0.53) and volcano (median = 0.52) showed the highest real-data stability, followed by Boruta (0.46) and stability selection (0.46). These values are generally higher than the simulation-derived estimates, consistent with the semi-synthetic findings and reflecting the richer correlation structure of real data.

The dataset with the highest average stability was ST001706 (*p/n* = 0.20; mean Nogueira = 0.57), while ST000369 (*p/n* = 2.26; mean = 0.18) had the lowest. The Spearman correlation between *p/n* ratio and stability across dataset--method combinations was −0.42 (*p* < 0.001), confirming the central role of dimensionality in determining FS stability. Dataset sample size also showed a positive association with stability (ρ = 0.36).

## Cross-database concordance is near zero

The most consequential finding of this study is the near-complete absence of cross-database concordance (Figure 8). Across the 67 evaluable method--pair combinations, the mean Spearman correlation of selection frequencies between datasets sharing the same platform was −0.02 (median = −0.02). The mean Jaccard similarity of selected feature sets was 0.06 (median = 0.00), with 64% of pairs sharing zero selected features.

The only method with non-negligible concordance was fold-change filtering (mean Jaccard = 0.40, mean Spearman = 0.003), but this reflects its tendency to select a large proportion of features (40--90% of common features), making overlap inevitable. Among methods with controlled FDR, all concordance metrics were centered on zero, with Spearman ρ values scattered uniformly between −0.30 and +0.27.

This result is devastating for the practice of feature-based biomarker discovery in metabolomics. Even when FS appears "stable" within a single study (Nogueira = 0.40--0.60 on real data), the selected features do not generalize to independent datasets of the same platform and biological context. Within-study bootstrap stability is therefore a necessary but grossly insufficient condition for genuine reproducibility.


# Discussion

## Summary of findings

This study provides the first systematic, metabolomics-specific benchmark of feature selection stability across 11 methods, 8 simulation scenarios, and 9 public datasets. Our three central findings are: (1) no method achieves consistently high stability, with median Nogueira indices of 0.19 across simulations and 0.18--0.57 (method medians) on real data; (2) the stability--accuracy trade-off is unfavorable, with methods that control FDR detecting fewer than 20% of true features at typical metabolomics effect sizes; and (3) cross-database concordance is essentially zero, meaning that feature selections do not reproduce across independent studies regardless of within-study stability.

## Comparison with prior work

Our results are consistent with, and extend, the genomics-focused benchmark of @haury2011influence, who reported stability values of 0.1--0.5 for most FS methods on microarray data. The metabolomics context compounds the problem: typical effect sizes are smaller (FC = 1.2--2.0 vs. 2--10 in gene expression), sample sizes are smaller (*n* = 50--300 vs. 100--1,000 in genomics), and the correlation structure is denser due to shared metabolic pathways. @bommert2020benchmark reached similar conclusions in their benchmark of 22 filter methods, observing that stability was strongly dependent on dataset characteristics and rarely exceeded 0.5. More recently, @hedou2024discovery demonstrated that standard FS in multi-omics settings produces unreliable biomarker signatures and proposed Stabl as a remedy. Our results are consistent with their findings but go further: even stability-enhanced approaches (our stability_selection, which implements the @meinshausen2010stability framework) remain unstable when *p* >> *n* and effect sizes are small. The instability problem flagged by @saeys2007review nearly two decades ago remains fundamentally unsolved.

## Why feature selection fails in metabolomics

Three structural properties of metabolomics data explain the observed instability. First, at FC = 1.5 (typical for many case-control comparisons), the signal-to-noise ratio is insufficient for any method to reliably distinguish 20 true features from 480 noise features in a *p* = 500, *n* = 100 setting. Our data show TPR < 0.05 for most methods at this effect size, meaning that selections are dominated by noise features and are therefore sensitive to sampling variation.

Second, the high inter-feature correlation characteristic of metabolomics data creates a Rashomon-like problem at the feature level. When multiple metabolites within the same pathway carry similar discriminative information, FS methods arbitrarily select different representatives depending on the bootstrap sample. However, our scenario S3 results suggest that, at realistic effect sizes, correlation per se is a secondary driver of instability --- the primary driver is signal weakness.

Third, the combination of multiple analytical and statistical choices --- preprocessing, normalization, significance thresholds, tuning parameters --- creates a "garden of forking paths" [@simmons2011false] that amplifies selection variability. Our preprocessing scenario (S7) showed that the choice of preprocessing pipeline altered stability rankings, although the magnitude of the effect was modest compared to the effect of dimensionality or effect size.

## Practical recommendations

Based on our findings and in alignment with the statistical framework articulated by @harrell2015regression, we offer the following recommendations for metabolomics researchers.

### If the goal is prediction: do not select features

Our data demonstrate that AUC is largely independent of FS stability: models built on unstable feature sets still achieve high predictive accuracy. This occurs because multiple distinct feature subsets carry equivalent predictive information. @harrell2015regression argues forcefully in Chapter 4 that variable selection degrades predictive models by introducing instability, optimism bias, and artificially narrow confidence intervals. We recommend using full-dimensional penalized models (ridge regression, PLS-DA) or ensemble methods (random forests, gradient boosting) without variable elimination when prediction is the primary objective.

### If the goal is biological interpretation: reduce, do not select

For biological insight, we recommend dimension reduction (PCA, PLS) combined with pathway-level interpretation of loadings, rather than single-metabolite selection. Loadings provide a continuous importance measure that does not force an arbitrary binary decision and is inherently more stable than binary selection. Pathway enrichment analysis (e.g., MSEA, mummichog) is preferable to single-metabolite selection because it aggregates evidence across correlated features within the same biological process, providing robustness to the specific feature subset that drives the signal [@xia2009metaboanalyst].

### If feature selection is unavoidable

When a reduced feature panel is required --- for example, for point-of-care diagnostic devices or cost-constrained targeted assays --- we recommend the following minimum standards:

1. **Report stability.** Compute and report the Nogueira index (with confidence intervals) alongside any selected feature list. A feature set with Nogueira < 0.5 should be treated as provisional. The R package `stabm` [@bommert2021stabm] provides efficient implementations.

2. **Require multi-method consensus.** Select only features identified by ≥3 independent methods. Our data show that single-method selections are not reproducible across either bootstrap resamples or independent datasets.

3. **Validate independently.** Any selected panel must be validated on a completely independent cohort. Our cross-database concordance results (mean Jaccard = 0.06) demonstrate unambiguously that within-study bootstrap stability does not predict cross-study generalizability.

4. **Use appropriate FDR control.** Avoid methods with empirical FDR > 10% (fold-change filtering, elastic net without post-selection inference). Prefer stability selection or wilcoxon + FDR when strict false positive control is needed.

## Semi-synthetic simulation as a benchmarking standard

Our semi-synthetic scenario (S8) produced more optimistic stability estimates than the parametric MVN simulations, because real data preserve distributional properties and correlation structures that provide additional information for FS methods. We recommend semi-synthetic simulation [@sankaran2025semisynthetic] as the preferred framework for future FS benchmarks in metabolomics, as it balances the need for known ground truth with realistic data properties. Purely parametric simulations risk underestimating method performance by failing to capture the structured covariance of real metabolomics data, while real-data evaluations lack ground truth for accuracy assessment.

## Limitations

Several limitations should be noted. First, our simulation framework, while extensive, cannot capture all complexities of real metabolomics data, including batch effects, analytical drift, ion suppression, and biological heterogeneity beyond case-control labels. Second, *B* = 100 bootstrap resamples may underestimate stability variance for some methods, although the jackknife confidence intervals provided by the Nogueira index partially address this concern. Third, we excluded the horseshoe prior [@heinze2018variable] for computational reasons (estimated ~60 min per fit at *p* = 500) and did not evaluate deep learning-based FS methods, which represent an increasingly important but computationally demanding approach. Fourth, cross-database concordance was limited by feature name matching across platforms; many platform pairs shared fewer than 120 features, potentially underestimating true concordance. Finally, the Nogueira index assigns equal weight to all features and does not distinguish between "nearly selected" and "never selected" features; a stability metric that incorporates selection frequency continuity might provide more nuanced insights.

## Conclusions

Feature selection in metabolomics produces unstable, non-reproducible results regardless of the method employed. The field should shift its focus from "which metabolites are important?" to "which metabolic patterns discriminate groups?" --- using dimension reduction and pathway-level analysis instead of single-feature selection. When feature selection is unavoidable, stability must be quantified and reported as a first-class metric alongside accuracy and FDR. The era of reporting a list of "top biomarkers" from a single FS run on a single dataset should end.


# Data availability

All public datasets are available from MetaboLights (<https://www.ebi.ac.uk/metabolights/>), Metabolomics Workbench (<https://www.metabolomicsworkbench.org/>), and the CIMCB GitHub repository (<https://github.com/CIMCB/MetabComparisonBinaryML>). Analysis code and configuration files are available at [repository URL].

# Funding

[To be added.]

# Conflict of interest

The authors declare no competing interests.

# References {-}

::: {#refs}
:::


# Figures {-}

![**Figure 1. Distribution of the Nogueira stability index across all simulated datasets.** Box plots show the distribution for each of 11 feature selection methods, ordered by median stability. The dashed line marks the chance level (0.5). Colors indicate method category: filter (orange), embedded (blue), wrapper (green), meta (yellow), Bayesian (dark blue), and ML-based (red). Note the strongly bimodal distributions of volcano, wilcoxon_fdr, and stability_selection, reflecting the "empty-set artifact" where methods select zero features at high *p/n* ratios.](../results/figures/fig01_stability_overview.png){#fig:fig01}

![**Figure 2. Stability--accuracy trade-off.** Scatter plot of mean stability (Nogueira index) versus mean sensitivity (TPR) for each method, averaged across all simulation scenarios. Point size is proportional to the average number of features selected. Dashed lines indicate the 0.5 reference for both axes. All methods except fold_change cluster in the lower-left quadrant, indicating simultaneously low stability and low sensitivity.](../results/figures/fig02_stability_accuracy_tradeoff.png){#fig:fig02}

![**Figure 3. Scenario-specific stability heatmap.** Mean Nogueira stability index for each combination of method (columns) and scenario level (rows), faceted by scenario (S1--S7 MVN, S8 semi-synthetic). Color scale ranges from 0 (dark purple) to 1 (yellow). The S8 panel is visibly brighter, reflecting the higher stability achievable on semi-synthetic data with preserved correlation structure.](../results/figures/fig03_scenario_heatmap.png){#fig:fig03}

![**Figure 4. Effect of dimensionality on stability.** Distribution of Nogueira stability index as a function of the *p/n* ratio (scenario S1), faceted by method. Note the paradoxical increase at pn50 for stability_selection, wilcoxon_fdr, and volcano: these methods select zero features at extreme dimensionality, achieving trivially perfect stability on the empty set.](../results/figures/fig04_pn_ratio.png){#fig:fig04}

![**Figure 5. Effect of inter-feature correlation on stability.** Distribution of Nogueira stability index as a function of the correlation scale factor (scenario S3), faceted by method. In contrast to expectations, increasing correlation does not consistently decrease stability, suggesting that signal weakness dominates over the "Rashomon effect" at realistic effect sizes.](../results/figures/fig05_correlation_effect.png){#fig:fig05}

![**Figure 6. Empirical false discovery rate control.** Distribution of FDR across all simulated datasets for each method, ordered by median FDR. The red dashed line marks the nominal 5% threshold. Fold-change filtering exhibits catastrophic FDR (>80%), while stability_selection, wilcoxon_fdr, and volcano maintain FDR well below 5%. The knockoff filter shows unexpected FDR inflation despite its theoretical guarantees, likely due to covariance mis-specification.](../results/figures/fig06_fdr_control.png){#fig:fig06}

![**Figure 7. Stability on real metabolomics data.** Distribution of the Nogueira stability index computed from 100 bootstrap resamplings on each of 9 public datasets, for each method. Method ordering and stability patterns are broadly consistent with simulation results, with generally higher absolute values reflecting the structured covariance of real data. Knockoff failed on 7 of 9 datasets (*n* < *p*).](../results/figures/fig07_real_data_stability.png){#fig:fig07}

![**Figure 8. Cross-database concordance.** Spearman correlation of feature selection frequencies between pairs of datasets sharing the same analytical platform, for each method. All distributions are centered on zero, indicating that feature rankings do not generalize across independent studies. This is the most consequential finding: within-study bootstrap stability does not predict cross-study reproducibility.](../results/figures/fig08_cross_database_concordance.png){#fig:fig08}
