# Discussion Structure — "The Stability Illusion"

## 1. Summary of findings (1 paragraph)

No FS method achieves consistently high stability across simulation scenarios (Nogueira index 0.15--0.46). The stability--accuracy trade-off is unfavorable: methods with higher TPR (fold_change) have catastrophic FDR (>80%), while methods with controlled FDR (wilcoxon_fdr, stability_selection) have TPR < 0.10 at realistic effect sizes (fc = 1.5). Cross-database concordance is near zero (Fig. 8), meaning feature rankings do not reproduce across independent datasets of the same platform.

## 2. Comparison with prior work (1--2 paragraphs)

- Haury et al. (2011): our results extend their genomics findings to metabolomics; similar stability range (0.1--0.5)
- Bommert et al. (2020): their filter benchmark showed comparable instability on high-dimensional classification data
- Hedou et al. (2024, Stabl): their Nature Biotech paper reached the same conclusion — standard FS is unreliable in omics — and proposed stability-enhanced LASSO. Our results are consistent but go further: even Stabl-like approaches (our stability_selection) remain unstable when p >> n
- Saeys et al. (2007): the instability problem they flagged 19 years ago remains unsolved

## 3. Why feature selection fails in metabolomics (1--2 paragraphs)

Three structural reasons:
1. **p >> n with small effects**: At fc = 1.5 (typical in metabolomics), the signal-to-noise ratio is too low for any method to reliably distinguish true from noise features. TPR < 0.03 means methods are essentially guessing.
2. **Correlated features**: Metabolites within the same pathway are highly correlated. When multiple features carry similar information, FS methods arbitrarily pick one — which one depends on the sample. This is the "Rashomon effect" applied to feature space.
3. **Multiple analytical choices**: Preprocessing, normalization, threshold parameters all influence which features are selected (Simmons et al. 2011). Our S7 (preprocessing scenario) confirms this.

## 4. What to do instead — Practical recommendations (2--3 paragraphs)

### 4.1 If the goal is prediction: don't select

Use full-dimensional penalized models (ridge regression, PLS-DA) or ensemble methods without variable elimination. Our data show that AUC is largely independent of stability (Fig. 2) — good prediction does not require identifying specific features. This aligns with Harrell (2015, Ch. 4), who argues that variable selection degrades predictive models by introducing instability and selection bias.

### 4.2 If the goal is biological interpretation: reduce, don't select

Use PCA or PLS and interpret loadings at the pathway level, not the individual metabolite level. Loadings provide a continuous importance measure that does not force an arbitrary binary decision. Follow up with targeted hypothesis tests on specific pathways suggested by the loading structure — but these are confirmatory tests with pre-specified hypotheses, not exploratory selection.

Pathway enrichment analysis (e.g., MSEA, mummichog) is preferable to single-metabolite selection because it aggregates evidence across correlated features within the same biological process, making results more robust to the exact feature set.

### 4.3 If feature selection is unavoidable (clinical panels)

When a reduced panel is required (e.g., point-of-care diagnostics), we recommend:
- **Report stability**: Always compute and report the Nogueira index alongside any selected feature list. A feature set with Nogueira < 0.5 should not be trusted as definitive.
- **Multi-method consensus**: Select only features identified by >=3 independent methods. Our data show that single-method selections are not reproducible.
- **Independent validation**: Any selected panel must be validated on a completely independent cohort before clinical claims. Cross-database concordance near zero (Fig. 8) shows that within-study bootstrap stability is not sufficient evidence of generalizability.

## 5. Semi-synthetic simulation as a benchmarking tool (1 paragraph)

S8 scenarios (spike-in on real data) produced more optimistic results than pure MVN simulations, because the realistic correlation structure provides additional information that FS methods can exploit. We recommend semi-synthetic simulation (Sankaran et al. 2025) as the preferred benchmarking framework for future FS comparisons in metabolomics, as it preserves the distributional and correlation properties that pure parametric simulations miss.

## 6. Limitations (1 paragraph)

- Simulation may not capture all complexities of real metabolomics data (batch effects, analytical drift, biological heterogeneity)
- 100 bootstrap resamples may underestimate stability variance for some methods
- Horseshoe prior excluded for computational reasons; deep learning-based FS methods not evaluated
- Cross-database concordance limited by feature name mismatches across platforms (0 common features for many pairs)

## 7. Conclusion (1 paragraph)

Feature selection in metabolomics produces unstable, non-reproducible results regardless of the method employed. The field should shift from "which metabolites are important?" to "which metabolic patterns discriminate groups?" — using dimension reduction and pathway-level analysis instead of single-feature selection. When FS is unavoidable, stability must be quantified and reported as a first-class metric alongside accuracy and FDR.

---

## Key references for the Discussion

| Claim | Reference |
|-------|-----------|
| Don't do variable selection | `harrell2015regression` (Ch. 4) |
| FS instability in genomics | `haury2011influence`, `boulesteix2009stability` |
| FS instability in omics (recent) | `hedou2024discovery` (Stabl/Nat Biotech) |
| Filter benchmark | `bommert2020benchmark` |
| FS taxonomy, open problem | `saeys2007review` |
| Researcher degrees of freedom | `simmons2011false` |
| Semi-synthetic simulation | `sankaran2025semisynthetic` |
| Statistical pitfalls in metabolomics | `broadhurst2006statistical` |
| Biomarker irreproducibility | `ransohoff2004rules`, `ioannidis2005why` |
| Stability metric | `nogueira2018stability` |
