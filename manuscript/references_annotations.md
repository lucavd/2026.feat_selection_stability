# References Annotations — "The Stability Illusion"

## Verification methodology

- **DOI-VERIFIED**: DOI resolved via Crossref/doi.org API and title confirmed matching (verified 2026-03-22)
- **NO-DOI**: No DOI exists (pre-DOI era, conference proceedings, JMLR, or CRAN manual); verified via alternative identifiers (PMCID, URL, CRAN page)

---

## Annotation Table

| # | BibTeX Key | Cite for | Section | Status |
|---|-----------|----------|---------|--------|
| **Cat. 1 — Stability metrics** | | | | |
| 1 | `nogueira2018stability` | Primary stability metric: chance-corrected index (Eq. 4), jackknife variance, CI; used for all 11 FS methods across 1590 datasets | Methods §Metrics; Introduction | NO-DOI (JMLR; URL verified) |
| 2 | `kuncheva2007stability` | First chance-corrected consistency index for feature subsets; secondary stability metric | Methods §Metrics; Introduction | NO-DOI (conf. proceedings; verified via citations) |
| 3 | `lustgarten2009measuring` | First to formally measure FS stability in biomedical context; historical motivation | Introduction | NO-DOI (AMIA; PMCID PMC2815476 verified) |
| 4 | `jaccard1901etude` | Jaccard similarity coefficient; pairwise set overlap metric | Methods §Metrics | DOI-VERIFIED (DOI 10.5169/seals-266450 — retro-DOI from SEALS) |
| 5 | `dice1945measures` | Dice coefficient as alternative pairwise similarity metric | Methods §Metrics | DOI-VERIFIED |
| 6 | `sorensen1948method` | Sørensen–Dice coefficient original reference | Methods §Metrics | NO-DOI (pre-DOI era) |
| 7 | `bommert2021stabm` | R package `stabm`: implementation of Nogueira, Kuncheva, and other stability indices | Methods §Implementation | DOI-VERIFIED (JOSS 6(59):3010) |
| **Cat. 2 — FS methods** | | | | |
| 8 | `wilcoxon1945individual` | Wilcoxon rank-sum test: univariate nonparametric filter | Methods §FS methods | DOI-VERIFIED |
| 9 | `mann1947test` | Mann–Whitney U test: equivalent rank-sum formulation | Methods §FS methods | DOI-VERIFIED |
| 10 | `benjamini1995controlling` | BH FDR correction for Wilcoxon p-values | Methods §FS methods | DOI-VERIFIED |
| 11 | `tibshirani1996regression` | LASSO: L1-penalized regression as embedded FS | Methods §FS methods | DOI-VERIFIED |
| 12 | `zou2005regularization` | Elastic net: L1+L2 penalization for correlated metabolites | Methods §FS methods | DOI-VERIFIED |
| 13 | `friedman2010regularization` | glmnet R package: coordinate descent for LASSO/elastic net | Methods §FS methods | DOI-VERIFIED |
| 14 | `kursa2010boruta` | Boruta: all-relevant FS via shadow features + RF | Methods §FS methods | DOI-VERIFIED |
| 15 | `breiman2001random` | Random forests: base learner for RF importance, Boruta | Methods §FS methods | DOI-VERIFIED |
| 16 | `wright2017ranger` | ranger: fast C++ RF implementation used throughout | Methods §Implementation | DOI-VERIFIED |
| 17 | `meinshausen2010stability` | Stability selection: subsampling + selection frequency framework | Methods §FS methods | DOI-VERIFIED |
| 18 | `hofner2015controlling` | stabs R package: implements stability selection with error control | Methods §FS methods | DOI-VERIFIED (BMC Bioinf 16:144) |
| 19 | `barber2015controlling` | Original knockoff filter: FDR control via synthetic variables | Methods §FS methods | DOI-VERIFIED |
| 20 | `candes2018panning` | Model-X knockoffs: extended knockoff filter; 94% failure when n < p | Methods §FS methods; Discussion | DOI-VERIFIED |
| 21 | `mitchell1988bayesian` | Spike-and-slab prior: original Bayesian variable selection | Methods §FS methods | DOI-VERIFIED |
| 22 | `george1993variable` | SSVS via Gibbs sampling: foundational Bayesian FS | Methods §FS methods | DOI-VERIFIED |
| 23 | `ishwaran2005spike` | Spike-and-slab theory: frequentist/Bayesian unification | Methods §FS methods | DOI-VERIFIED (Ann Stat 33(2):730-773) |
| 24 | `ishwaran2010spikeslab` | spikeslab R package: implementation used in our pipeline | Methods §Implementation | DOI-VERIFIED (R Journal 2(2):68-73) |
| 25 | `lundberg2017unified` | SHAP: Shapley-value feature importance from XGBoost | Methods §FS methods | NO-DOI (NeurIPS 2017; verified via proceedings) |
| 26 | `chen2016xgboost` | XGBoost: gradient-boosted tree base model for SHAP | Methods §FS methods | DOI-VERIFIED |
| 27 | `li2012volcano` | Volcano plot as combined fold-change + p-value criterion | Methods §FS methods | DOI-VERIFIED (J Bioinf Comput Biol 10(6):1231003) |
| **Cat. 3 — FS stability literature** | | | | |
| 28 | `saeys2007review` | Canonical FS taxonomy (filter/wrapper/embedded); instability as open problem | Introduction | DOI-VERIFIED |
| 29 | `he2010stable` | Stability–accuracy tradeoff in biomarker discovery | Introduction; Discussion | DOI-VERIFIED |
| 30 | `boulesteix2009stability` | FS instability in high-dimensional ranked gene lists (Brief. Bioinf.) | Introduction | DOI-VERIFIED |
| 31 | `haury2011influence` | Benchmark of 32 FS methods: stability, accuracy, interpretability | Introduction; Discussion | DOI-VERIFIED |
| 32 | `hedou2024discovery` | Stabl: recent evidence standard FS in omics is unreliable; stability-enhanced LASSO | Introduction; Discussion | DOI-VERIFIED (Nat Biotech 42(10):1581-1593) |
| 33 | `bommert2020benchmark` | Large-scale benchmark of 22 filter methods across 16 datasets | Introduction; Discussion | DOI-VERIFIED |
| **Cat. 4 — Metabolomics** | | | | |
| 34 | `mendez2019comparative` | Source of 7 CIMCB benchmark metabolomics datasets | Methods §Datasets | DOI-VERIFIED |
| 35 | `haug2020metabolights` | MetaboLights database: hosts MTBLS28, MTBLS404, etc. | Methods §Datasets | DOI-VERIFIED |
| 36 | `sud2016metabolomics` | Metabolomics Workbench: hosts ST001706 and ST-prefixed datasets | Methods §Datasets | DOI-VERIFIED |
| 37 | `sumner2007proposed` | MSI minimum reporting standards | Introduction | DOI-VERIFIED |
| 38 | `broadhurst2006statistical` | Statistical pitfalls in metabolomics; motivates stability assessment | Introduction | DOI-VERIFIED |
| 39 | `ransohoff2004rules` | Biomarker reproducibility crisis; molecular marker validation rules | Introduction | DOI-VERIFIED |
| 40 | `dunn2011procedures` | Large-scale metabolomics QC procedures | Introduction | DOI-VERIFIED |
| 41 | `alseekh2021mass` | MS-based metabolomics best-practice guidelines | Introduction | DOI-VERIFIED |
| 42 | `xia2009metaboanalyst` | MetaboAnalyst: widely-used platform with volcano plot as standard FS | Introduction; Methods | DOI-VERIFIED |
| **Cat. 5 — Simulation** | | | | |
| 43 | `genz2009computation` | mvtnorm: multivariate normal generation for scenarios S1–S7 | Methods §Simulation | DOI-VERIFIED |
| 44 | `schafer2005shrinkage` | corpcor: shrinkage covariance estimation for simulation input | Methods §Simulation | DOI-VERIFIED |
| 45 | `ledoit2004well` | Ledoit-Wolf shrinkage: theoretical basis for well-conditioned covariance | Methods §Simulation | DOI-VERIFIED |
| 46 | `vandenberg2006centering` | Metabolomics data transformations incl. log-transform; justifies MVN → exp() | Methods §Simulation | DOI-VERIFIED |
| 47 | `sankaran2025semisynthetic` | Semi-synthetic simulation framework preserving real data structure; justifies S8 | Methods §Simulation; Discussion | DOI-VERIFIED (Brief Bioinf 26(1):bbaf051) |
| 48 | `franklin2014plasmode` | Plasmode simulation: realistic benchmarks from observed data | Methods §Simulation | DOI-VERIFIED |
| **Cat. 6 — Prediction/evaluation** | | | | |
| 49 | `robin2011proc` | pROC R package: AUC computation with DeLong CIs | Methods §Evaluation | DOI-VERIFIED |
| 50 | `matthews1975comparison` | MCC: balanced accuracy metric | Methods §Evaluation | DOI-VERIFIED |
| 51 | `chicco2020advantages` | MCC superiority over F1; justifies MCC as primary metric | Methods §Evaluation | DOI-VERIFIED |
| 52 | `fawcett2006introduction` | ROC analysis: AUC interpretation, threshold-independent evaluation | Methods §Evaluation | DOI-VERIFIED |
| **Cat. 7 — Software** | | | | |
| 53 | `rcoreteam2024r` | R environment: all analyses conducted in R | Methods §Implementation | NO-DOI (standard R citation; URL verified) |
| 54 | `wickham2016ggplot2` | ggplot2: all figures | Methods §Implementation | DOI-VERIFIED |
| 55 | `dowle2024datatable` | data.table: high-performance data manipulation | Methods §Implementation | NO-DOI (CRAN manual) |
| **Cat. 8 — General context** | | | | |
| 56 | `guyon2003introduction` | Canonical FS taxonomy and problem formulation (JMLR) | Introduction | NO-DOI (JMLR; URL verified) |
| 57 | `heinze2018variable` | Variable selection pitfalls and recommendations | Introduction; Discussion | DOI-VERIFIED |
| 58 | `simmons2011false` | "Researcher degrees of freedom": analytical flexibility → false positives | Introduction; Discussion | DOI-VERIFIED |
| 59 | `ioannidis2005why` | Most research findings are false: contextualizes biomarker irreproducibility | Introduction | DOI-VERIFIED |
| 60 | `li2022benchmark` | Multi-omics FS benchmark: block vs. concurrent strategies | Discussion | DOI-VERIFIED |
| 61 | `efron2004large` | Large-scale simultaneous testing: empirical null, FDR behavior | Introduction; Discussion | DOI-VERIFIED |
| 62 | `harrell2015regression` | "Don't do variable selection": authoritative case against FS in predictive modeling; Ch. 4 on dimensionality reduction as alternative | Discussion §Recommendations | DOI-VERIFIED |

---

## Verification Summary

| Status | Count | Meaning |
|--------|-------|---------|
| **DOI-VERIFIED** | 53 | DOI resolved via Crossref/doi.org API, title confirmed matching |
| **NO-DOI** | 9 | No DOI exists; verified via PMCID, URL, or CRAN page |
| **Total** | **62** | |

## Verification details

All 53 DOIs were programmatically resolved on 2026-03-22 via `curl -H "Accept: application/json" https://doi.org/{doi}`. Every resolved title matched the corresponding BibTeX entry. The Jaccard (1901) retro-DOI (`10.5169/seals-266450`) resolves via the Swiss e-periodica.ch platform rather than Crossref, but title and metadata are confirmed. The 9 NO-DOI entries are legitimate (pre-DOI era, JMLR, NeurIPS proceedings, AMIA, CRAN manuals) and were verified via alternative identifiers.
