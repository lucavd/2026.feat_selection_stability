# Inventario Metodi di Feature Selection

> Ultimo aggiornamento: 2026-03-18 (rev. 2)
> Status: DRAFT — horseshoe rimosso, shap_xgboost aggiornato a GPU

---

## Overview

11 metodi di feature selection organizzati in 5 categorie (6 originali, horseshoe rimosso).

**Nota (2026-03-18):** Horseshoe prior rimosso per insostenibilità computazionale (~60 min/job su p=500 con 100 bootstrap MCMC). Spike-slab copre la categoria bayesiana. SHAP aggiornato da `treeshap` (CPU, >150 min) a `predict(predcontrib = TRUE)` con xgboost GPU (35s). Dettagli in `docs/07_implementation_status.md` §5.3 e §5b.

---

## 1. FILTER METHODS

### 1.1 Wilcoxon Rank-Sum Test + FDR

- **Tipo:** Filter, univariato, non parametrico
- **Pacchetto R:** `stats` (base R, built-in)
- **Funzione:** `wilcox.test()` + `p.adjust(method = "fdr")`
- **Alternativa metabolomica:** `structToolbox::wilcox_test()` con parametro `mtc = "fdr"`
- **Implementazione:**
  ```r
  # Per ogni feature j:
  pvals[j] <- wilcox.test(x[group == "case", j], x[group == "control", j])$p.value
  # Correzione multipla:
  adj_pvals <- p.adjust(pvals, method = "BH")  # Benjamini-Hochberg
  selected <- which(adj_pvals < 0.05)
  ```
- **Parametri:** soglia FDR (default 0.05)
- **Vantaggi:** Nessuna assunzione distribuzionale, robusto a outlier
- **Limitazioni:** Univariato — ignora struttura di correlazione; non cattura interazioni
- **Paper:** Benjamini & Hochberg (1995), JRSS-B 57(1):289-300

### 1.2 Fold-Change Thresholding

- **Tipo:** Filter, univariato, effect-size based
- **Pacchetto R:** Base R (nessun pacchetto dedicato necessario)
- **Implementazione:**
  ```r
  fc <- colMeans(x[group == "case", ]) / colMeans(x[group == "control", ])
  log2fc <- log2(fc)
  selected <- which(abs(log2fc) > log2(1.5))  # FC > 1.5
  ```
- **Parametri:** soglia FC (tipicamente 1.5 o 2.0)
- **Vantaggi:** Semplice, interpretabile clinicamente
- **Limitazioni:** Ignora significatività statistica; sensibile a outlier; il metodo "barbarico" che molti usano in pratica
- **Note:** Nella metabolomica la maggioranza dei FC è < 1.5 (letteratura conferma); FC > 2 sono rari

### 1.3 Volcano Plot-Based Selection

- **Tipo:** Filter, univariato, combinazione FC + p-value
- **Pacchetto R:** Base R (plot custom) o `EnhancedVolcano` (Bioconductor)
- **Implementazione:**
  ```r
  # Combina p-value e fold-change:
  selected <- which(adj_pvals < 0.05 & abs(log2fc) > log2(1.5))
  ```
- **Parametri:** soglia p-value (0.05), soglia FC (1.5)
- **Vantaggi:** Bilancia significatività statistica e rilevanza biologica
- **Limitazioni:** Ancora univariato; soglie arbitrarie; varianti robuste disponibili per dati con outlier (BMC Bioinformatics)

---

## 2. WRAPPER / EMBEDDED METHODS

### 2.1 LASSO (L1-penalized regression)

- **Tipo:** Embedded, lineare, con sparsità
- **Pacchetto R:** `glmnet` (CRAN)
- **Funzioni:** `glmnet()`, `cv.glmnet()`
- **Implementazione:**
  ```r
  library(glmnet)
  fit <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)
  coefs <- coef(fit, s = "lambda.1se")
  selected <- which(coefs[-1] != 0)
  ```
- **Parametri chiave:**
  - `alpha = 1` (LASSO puro)
  - `lambda`: scelto via CV — `lambda.min` (errore minimo) vs `lambda.1se` (più parsimonioso, default raccomandato)
  - `nfolds`: numero di fold per CV (default 10)
- **Problema instabilità:** La selezione di lambda via CV è stocastica; ripetendo la CV si ottengono set diversi. Soluzione: ripetere CV multiple volte e aggregare, o usare `escv.glmnet` per "estimation stability with cross-validation"
- **Limitazioni:** Con features altamente correlate, seleziona arbitrariamente una del gruppo; coefficienti biased (shrinkage verso zero)
- **Paper:** Tibshirani (1996), JRSS-B 58(1):267-288

### 2.2 Elastic Net (L1 + L2 penalty)

- **Tipo:** Embedded, lineare
- **Pacchetto R:** `glmnet` (CRAN)
- **Implementazione:**
  ```r
  fit <- cv.glmnet(x, y, family = "binomial", alpha = 0.5, nfolds = 10)
  ```
- **Parametri chiave:**
  - `alpha` ∈ (0, 1): bilancia L1 (sparsità) e L2 (grouping effect). Alpha=0.5 è un buon default
  - Cerca anche `alpha` ottimale via grid search
- **Vantaggi rispetto a LASSO:** Gestisce meglio feature correlate (tende a selezionare gruppi); più stabile numericamente con `alpha = 1 - ε`
- **Quando preferirlo:** Dati metabolomici con forti correlazioni (NMR in particolare)
- **Paper:** Zou & Hastie (2005), JRSS-B 67(2):301-320

### 2.3 Boruta

- **Tipo:** Wrapper, all-relevant feature selection
- **Pacchetto R:** `Boruta` v9.0.0 (CRAN)
- **Funzione:** `Boruta()`
- **Implementazione:**
  ```r
  library(Boruta)
  bor <- Boruta(x, y, maxRuns = 500, doTrace = 0)
  selected <- which(bor$finalDecision == "Confirmed")
  # Per features "Tentative":
  bor_fixed <- TentativeRoughFix(bor)
  ```
- **Meccanismo shadow features:**
  1. Crea copie permutate (shadow) di tutte le features
  2. Addestra Random Forest su features originali + shadow
  3. Confronta importance di ogni feature originale vs max importance tra le shadow
  4. Iterativamente classifica: **Confirmed** (significativamente > max shadow), **Rejected** (significativamente <), **Tentative** (indeciso)
- **Parametri chiave:**
  - `maxRuns`: iterazioni massime (default 100, raccomandato 500 per stabilità)
  - `pValue`: soglia per il test binomiale (default 0.01)
  - `mcAdj`: correzione Bonferroni (default TRUE)
  - Parametri del Random Forest sottostante (ntree, mtry)
- **Limitazioni note:**
  - Computazionalmente intensivo (O(maxRuns × RF training))
  - Sensibile agli iperparametri del RF sottostante
  - Tende a selezionare gruppi di features correlate (tutte o nessuna)
  - Non garantisce stabilità al resampling — questo è un punto chiave del nostro studio
- **Paper:** Kursa & Rudnicki (2010), J Statistical Software 36(11)

### 2.4 Random Forest Importance

- **Tipo:** Wrapper/filter, non-lineare
- **Pacchetto R:** `ranger` (CRAN, raccomandato) o `randomForest` (CRAN, classico)
- **Implementazione:**
  ```r
  library(ranger)
  rf <- ranger(y ~ ., data = cbind(x, y = y), importance = "permutation", num.trees = 1000)
  imp <- rf$variable.importance
  selected <- names(sort(imp, decreasing = TRUE))[1:k]  # top-k features
  ```
- **Tipi di importance:**
  - **Permutation importance** (`importance = "permutation"`): più affidabile, misura il decremento di accuracy quando si permuta una feature. Computazionalmente più costoso.
  - **Impurity importance (MDI)** (`importance = "impurity"`): Mean Decrease in Gini/RSS. Bias noto verso variabili continue e ad alta cardinalità; calcolata su training data.
  - **Corrected impurity** (`importance = "impurity_corrected"` in ranger): corregge il bias verso features ad alta cardinalità.
- **Raccomandazione:** Usare **permutation importance** per il nostro studio (meno bias, più comparabile)
- **Limitazioni:** Instabile con features correlate (importance distribuita tra features ridondanti); la soglia per "selezionato" è arbitraria (top-k o importance > threshold)
- **Paper:** Breiman (2001), Machine Learning 45(1):5-32

### 2.5 Stability Selection

- **Tipo:** Meta-metodo, applicabile a qualsiasi base learner
- **Pacchetto R:** `stabs` v0.7-1 (CRAN, aggiornato 2026-01-31)
- **Autori pacchetto:** Benjamin Hofner, Torsten Hothorn
- **Implementazione:**
  ```r
  library(stabs)
  stab <- stabsel(x, y, fitfun = glmnet.lasso, cutoff = 0.75,
                  PFER = 1, sampling.type = "SS")
  selected <- stab$selected
  ```
- **Meccanismo:**
  1. Ricampionamento (subsampling) ripetuto del dataset
  2. Per ogni sottocampione, applica il metodo base (es. LASSO)
  3. Calcola la frequenza di selezione di ogni feature
  4. Seleziona features con frequenza > cutoff (es. 0.75)
- **Parametri chiave:**
  - `cutoff`: soglia di frequenza (0.6-0.9, default spesso 0.75)
  - `PFER`: Per-Family Error Rate target (controlla il numero atteso di falsi positivi)
  - `sampling.type`: "SS" (Meinshausen & Bühlmann) o "MB" (complementary pairs, Shah & Samworth 2013)
  - `fitfun`: funzione di fitting base (lasso, boosting, ecc.)
- **Vantaggi teorici:**
  - Controllo errore in campione finito (PFER bound)
  - Consistenza della selezione anche quando le condizioni di consistenza del LASSO sono violate
  - Stabilità per design — il metodo è esplicitamente costruito per essere stabile
- **Limitazioni:** Più conservativo (meno features selezionate); computazionalmente più costoso (B × fit del base learner)
- **Paper:** Meinshausen & Bühlmann (2010), JRSS-B 72(4):417-473; Shah & Samworth (2013), JRSS-B 75(1):55-80

### 2.6 Knockoff Filter

- **Tipo:** Embedded, controllo FDR esatto
- **Pacchetto R:** `knockoff` (CRAN, aggiornato 2025-07-22)
- **Implementazione:**
  ```r
  library(knockoff)
  result <- knockoff.filter(x, y, fdr = 0.1, statistic = stat.glmnet_coefdiff)
  selected <- result$selected
  ```
- **Meccanismo:**
  1. Genera "knockoff" features: variabili sintetiche con stessa struttura di correlazione ma indipendenti dalla risposta
  2. Confronta l'importanza di ogni feature con il suo knockoff
  3. Seleziona features la cui importanza supera significativamente quella del knockoff
- **Varianti:**
  - **Fixed-X knockoffs:** richiedono n ≥ p (non praticabili in metabolomica tipica)
  - **Model-X knockoffs:** per p >> n, assumono distribuzione multivariata nota (tipicamente normale); stimano media e covarianza dai dati
- **Parametri chiave:**
  - `fdr`: target FDR (default 0.1)
  - `statistic`: statistica di test (glmnet_coefdiff, random_forest, ecc.)
  - Stima della covarianza per model-X
- **Limitazioni per metabolomica:**
  - Richiede stima accurata della covarianza (difficile con p >> n)
  - Assume normalità multivariata per la generazione dei knockoff — può fallire con dati non-gaussiani
  - Computazionalmente costoso per p molto grande
  - Power ridotto quando l'assunzione di normalità è violata
- **Paper:** Barber & Candès (2015), Annals of Statistics 43(5):2055-2085; Candès et al. (2018), JRSS-B 80(3):551-577

---

## 3. BAYESIAN METHODS

### 3.1 Horseshoe Prior

- **Tipo:** Bayesian shrinkage, continuous shrinkage prior
- **Pacchetti R:**
  - `horseshoe` (CRAN) — implementazione dedicata, algoritmo Bhattacharya et al. (2016)
  - `brms` (CRAN) — `horseshoe()` prior in `set_prior()`
  - `rstanarm` (CRAN) — `stan_glm()` con horseshoe prior specification

- **Implementazione con `horseshoe`:**
  ```r
  library(horseshoe)
  fit <- horseshoe(y, x, method.tau = "halfCauchy", method.sigma = "Jeffreys",
                   burn = 1000, nmc = 5000)
  # Feature selection via credible intervals:
  selected <- which(fit$LeftCI * fit$RightCI > 0)  # CI non include zero
  ```

- **Implementazione con `brms`:**
  ```r
  library(brms)
  fit <- brm(y ~ ., data = df, family = bernoulli(),
             prior = set_prior(horseshoe(df = 1, scale_global = 0.5)))
  ```

- **Meccanismo:** Spike infinitamente alto a zero con code pesanti; segnali forti rimangono non-shrunk, features nulle fortemente penalizzate verso zero
- **Parametri chiave:**
  - `scale_global` (τ): shrinkage globale. Regola empirica: `(p0/p) / sqrt(n)` dove p0 = n. atteso di features non-nulle. Piironen & Vehtari (2017) raccomandano 0.5 come weakly informative
  - `df`: gradi di libertà (1 = Cauchy, default per horseshoe classico)
  - Burn-in e MCMC iterations per convergenza
- **Vantaggi:** Adattivo — shrinkage diverso per ogni feature; ottimo per sparsità reale
- **Limitazioni:** Computazionalmente costoso (MCMC); tempo scala con p e n; sensibile alla scelta di `scale_global`
- **Paper:** Carvalho et al. (2010), Biometrika 97(2):465-480; Piironen & Vehtari (2017), Electronic Journal of Statistics 11(2):5018-5051

### 3.2 Spike-and-Slab Prior

- **Tipo:** Bayesian, variable selection prior discreto
- **Pacchetti R:**
  - `spikeslab` (CRAN) — generalized elastic net in spike-and-slab framework, **scalabile a p grande**
  - `BoomSpikeSlab` (CRAN) — MCMC-based, errori Gaussiani o Student-t
  - `spikeSlabGAM` (CRAN) — estensione per GAM con prior peNMIG
  - `BayesSUR` (CRAN) — seemingly unrelated regression con spike-and-slab MRF priors

- **Implementazione con `spikeslab`:**
  ```r
  library(spikeslab)
  fit <- spikeslab(y ~ ., data = df, bigp.smalln = TRUE, max.var = 100)
  selected <- which(fit$gnet.scale != 0)
  ```

- **Meccanismo:** Prior mistura: "spike" (distribuzione concentrata a zero) e "slab" (distribuzione diffusa per effetti non-nulli). Ogni coefficiente ha una probabilità di inclusione.
- **Parametri chiave:**
  - Prior probability of inclusion (π)
  - Varianza dello slab
  - MCMC: burn-in, iterations, thinning
- **Fattibilità per p=1000+:** `spikeslab` è esplicitamente progettato per questo regime (`bigp.smalln = TRUE`). `BoomSpikeSlab` funziona ma è più lento per p molto grande.
- **Vantaggi:** Interpretazione probabilistica naturale (posterior probability of inclusion); inferenza completa su ogni coefficiente
- **Limitazioni:** Computazionalmente intensivo (MCMC); mixing può essere lento con p molto grande; scelta dei prior iperparametri non banale
- **Paper:** Mitchell & Beauchamp (1988), JASA 83(404):1023-1032; Ishwaran & Rao (2005), Annals of Applied Statistics

---

## 4. ML-BASED METHODS

### 4.1 SHAP-Based Feature Selection

- **Tipo:** Model-agnostic feature importance, post-hoc
- **Pacchetti R:**
  - `shapviz` v0.10.3 (CRAN, aggiornato 2025-10-13) — visualizzazione SHAP
  - `treeshap` (GitHub: ModelOriented/treeshap) — calcolo efficiente TreeSHAP per ensemble tree-based
- **Pacchetti Python:** `shap` (PyPI, slundberg/shap)

- **Implementazione R:**
  ```r
  library(xgboost)
  library(shapviz)

  # Train XGBoost
  dtrain <- xgb.DMatrix(data = x, label = y)
  model <- xgb.train(params = list(objective = "binary:logistic", max_depth = 6),
                      data = dtrain, nrounds = 100)

  # Calcola SHAP values
  shp <- shapviz(model, X_pred = x)

  # Feature importance = mean |SHAP|
  imp <- colMeans(abs(shp$S))
  selected <- names(sort(imp, decreasing = TRUE))[1:k]
  ```

- **Meccanismo:** SHAP values decompongono la predizione in contributi additivi per ogni feature, basati sulla teoria dei giochi (Shapley values)
- **Parametri chiave:**
  - Modello sottostante (XGBoost, LightGBM, RF, ecc.)
  - Iperparametri del modello (influenzano le SHAP values)
  - Soglia/top-k per selezione
- **Per stabilizzare:** Ripetere training su bootstrap + aggregare SHAP values
- **Vantaggi:** Corregge il bias della Gini importance; interpretabile; computabile su validation set
- **Limitazioni:** Dipende dal modello sottostante; computazionalmente costoso per p >> n senza pre-screening; non ha garanzie teoriche di controllo errore
- **Paper:** Lundberg & Lee (2017), NeurIPS

### 4.2 Sparse Autoencoders (Opzionale — Python)

- **Tipo:** Unsupervised feature extraction/selection, non-lineare
- **Framework:** PyTorch o Keras/TensorFlow
- **Meccanismo:** Autoencoder con penalità di sparsità (L1 o KL divergence) sul bottleneck layer. Features "importanti" sono quelle con attivazione non-nulla nello strato latente.
- **Implementazione di riferimento:** Diversi repository GitHub disponibili (K-Sparse-AutoEncoder, ecc.)
- **Limitazioni:** Approccio unsupervised — potrebbe non catturare features rilevanti per il fenotipo; richiede tuning dell'architettura; meno interpretabile
- **Nota:** Consideriamo questo metodo come **opzionale** dato che aggiunge complessità (Python + GPU) senza garanzie di superiorità in questo contesto

---

## 5. ENSEMBLE / META METHODS

### 5.1 Ensemble Feature Selection

- **Concetto:** Combinare i risultati di più metodi di feature selection per aumentare la stabilità
- **Approcci:**
  - **Union/Intersection:** selezione di features che appaiono in ≥k metodi su M totali
  - **Rank aggregation:** media (o mediana) dei ranking di importanza tra metodi diversi
  - **MVFS-SHAP** (2025): majority voting con SHAP values per valutazione multi-dimensionale

- **Implementazione custom in R:**
  ```r
  # Dopo aver applicato N metodi, ciascuno con un vettore di features selezionate:
  all_selections <- list(sel_wilcox, sel_lasso, sel_boruta, sel_rf, ...)
  # Contare frequenza di selezione:
  freq <- table(unlist(all_selections))
  # Selezionare features in almeno k metodi:
  consensus <- names(freq[freq >= k])
  ```

- **Paper recente:** MVFS-SHAP: Majority Voting with SHAP for Stable Feature Selection (2025, Computers in Biology and Medicine)

---

## 6. Riepilogo Pacchetti R Necessari

```r
# CRAN packages:
install.packages(c(
  "glmnet",        # LASSO, Elastic Net
  "Boruta",        # Boruta feature selection
  "ranger",        # Random Forest (fast)
  "stabs",         # Stability selection
  "knockoff",      # Knockoff filter
  "horseshoe",     # Horseshoe prior
  "spikeslab",     # Spike-and-slab
  "brms",          # Bayesian regression (horseshoe)
  "xgboost",       # Per SHAP
  "shapviz",       # SHAP visualization
  "stabm"          # Stability metrics
))

# GitHub packages:
# remotes::install_github("ModelOriented/treeshap")

# Bioconductor (opzionale):
# BiocManager::install("OmicsMarkeR")
```

---

## 7. Confronto Proprietà dei Metodi

| Metodo | Lineare | Gestisce correlazioni | Controllo errore | Complessità | Stabilità attesa |
|--------|---------|----------------------|------------------|-------------|-----------------|
| Wilcoxon | Sì (univ.) | No | FDR (BH) | O(n·p) | Media |
| Fold-change | Sì (univ.) | No | Nessuno | O(n·p) | Bassa |
| Volcano | Sì (univ.) | No | FDR + FC | O(n·p) | Media |
| LASSO | Sì | Parziale | No formale | O(n·p²) | Bassa |
| Elastic Net | Sì | Sì (grouping) | No formale | O(n·p²) | Media |
| Boruta | No | Parziale | Test binomiale | O(B·RF) | Media |
| RF importance | No | Parziale | Nessuno | O(B·n·p·log(n)) | Media-bassa |
| Stability sel. | Dipende | Dipende dal base | PFER bound | O(B²·base) | **Alta** |
| Knockoff | Sì | Sì (model-X) | FDR esatto | O(p³) | Media-alta |
| Horseshoe | Sì | Sì (Bayes) | Posterior CI | O(MCMC·p²) | Media-alta |
| Spike-and-slab | Sì | Sì (Bayes) | PPI | O(MCMC·p²) | Media-alta |
| SHAP-based | No | No (data-driven) | Nessuno | O(model·n·p) | Dipende |
