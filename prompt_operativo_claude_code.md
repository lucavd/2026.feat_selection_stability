# PROMPT OPERATIVO — Claude Code: Feature Selection Stability in Metabolomics

> Generato il: 2026-03-17
> Basato su: ricerca sistematica documentata in `docs/01-06`
> Destinazione: Claude Code per generazione codice R completo

---

## Chi sei

Sei un programmatore R esperto con competenze in biostatistica, chemiometria e metabolomica. Scrivi codice pulito, documentato, riproducibile e modulare. Non inventi dati e non usi placeholder: ogni funzione deve essere completa e funzionante.

## Contesto del progetto

Articolo scientifico in silico:
**"The stability illusion: A simulation-based assessment of feature selection methods for high-dimensional metabolomics data"**

Target journal: Briefings in Bioinformatics / Metabolomics / Analytical Chemistry.

L'obiettivo è confrontare sistematicamente 12 metodi di feature selection su dati metabolomici simulati (con struttura di correlazione estratta da dati reali), misurando stabilità (Nogueira index), accuratezza (TPR/FDR), performance predittiva (AUC), e parsimonia.

---

## Struttura del Progetto

Crea esattamente questa struttura:

```
2026.feat_selection_stability/
├── config/
│   └── config.yaml                    # TUTTI i parametri dello studio
├── R/
│   ├── 00_install_packages.R          # Installazione dipendenze
│   ├── 01_download_data.R             # Download dataset pubblici
│   ├── 02_preprocess.R                # Preprocessing e armonizzazione
│   ├── 03_extract_empirical_params.R  # Estrazione parametri da dati reali
│   ├── 04_simulate.R                  # Framework di simulazione
│   ├── 05_feature_selection.R         # 12 metodi di feature selection
│   ├── 06_metrics.R                   # Calcolo metriche
│   ├── 07_cross_validation.R          # Validazione cross-database
│   ├── 08_figures.R                   # Figure e tabelle per il paper
│   └── utils/
│       ├── helpers.R                  # Funzioni helper condivise
│       ├── fs_methods.R               # Wrapper per ogni metodo FS
│       └── stability_metrics.R        # Wrapper per metriche stabilità
├── data/
│   ├── raw/                           # Dataset scaricati (gitignored)
│   ├── processed/                     # Dataset armonizzati
│   ├── empirical_params/              # Parametri empirici estratti
│   └── simulated/                     # Dataset simulati (gitignored)
├── results/
│   ├── feature_selection/             # Risultati FS per ogni scenario
│   ├── metrics/                       # Metriche aggregate
│   ├── checkpoints/                   # Checkpoint per ripresa
│   └── figures/                       # Figure generate
├── docs/                              # Documentazione (già presente)
├── .gitignore
├── LICENSE.md
└── README.md
```

---

## File di Configurazione Centrale

### `config/config.yaml`

```yaml
# ============================================================
# CONFIGURAZIONE CENTRALE — Feature Selection Stability Study
# ============================================================
# Tutti i parametri dello studio definiti in un unico posto.
# NON hardcodare parametri negli script — leggere sempre da qui.

project:
  name: "feat_selection_stability"
  seed: 42
  n_cores: 8  # Modificare in base all'hardware

# ----------------------------------------------------------
# DATASET PUBBLICI
# ----------------------------------------------------------
datasets:
  metabolights:
    # Dataset MetaboLights caso-controllo
    # NOTA: gli accession number devono essere verificati interrogando l'API
    # prima della run. Lo Script 01 verifica automaticamente l'accessibilità.
    accessions:
      - id: "MTBLS1"        # Diabete tipo 2 — NMR, da verificare dettagli
      - id: "MTBLS733"      # Benchmark spike-in — LC-MS, ground truth nota
      # Aggiungere 4-6 dataset dopo verifica API (vedi docs/01)
    api_base: "https://www.ebi.ac.uk/metabolights/ws"

  metabolomics_workbench:
    accessions:
      - id: "ST000001"      # Placeholder — da sostituire con studi verificati
      # Aggiungere 3-5 dataset dopo verifica API (vedi docs/01)
    api_base: "https://www.metabolomicsworkbench.org/rest"

  # Dataset benchmark con ground truth nota:
  benchmark:
    - id: "MTBLS733"
      platform: "LC-MS"
      n_true_features: 130
      description: "Piper nigrum standard mixtures, Thermo Q Exactive HF"

# ----------------------------------------------------------
# METODI DI FEATURE SELECTION
# ----------------------------------------------------------
methods:
  filter:
    - name: "wilcoxon_fdr"
      params:
        fdr_threshold: 0.05

    - name: "fold_change"
      params:
        fc_threshold: 1.5  # Fold-change minimo

    - name: "volcano"
      params:
        fdr_threshold: 0.05
        fc_threshold: 1.5

  embedded:
    - name: "lasso"
      params:
        alpha: 1.0
        lambda_rule: "lambda.1se"  # "lambda.min" o "lambda.1se"
        nfolds: 10

    - name: "elastic_net"
      params:
        alpha: 0.5
        lambda_rule: "lambda.1se"
        nfolds: 10

    - name: "knockoff"
      params:
        fdr: 0.1
        statistic: "stat.glmnet_coefdiff"

  wrapper:
    - name: "boruta"
      params:
        maxRuns: 300
        pValue: 0.01

    - name: "rf_importance"
      params:
        num_trees: 1000
        importance: "permutation"
        top_k_method: "elbow"  # "elbow", "top_n", "threshold"
        top_n: 20

  meta:
    - name: "stability_selection"
      params:
        cutoff: 0.75
        PFER: 1
        sampling_type: "SS"
        B: 100

  bayesian:
    - name: "horseshoe"
      params:
        method_tau: "halfCauchy"
        method_sigma: "Jeffreys"
        burn: 1000
        nmc: 3000

    - name: "spike_slab"
      params:
        bigp_smalln: true
        max_var: 100

  ml:
    - name: "shap_xgboost"
      params:
        max_depth: 6
        nrounds: 200
        eta: 0.1
        top_k_method: "elbow"

# ----------------------------------------------------------
# SCENARI DI SIMULAZIONE
# ----------------------------------------------------------
simulation:
  # Scenario base (default)
  base:
    n_case: 50
    n_control: 50
    p: 500
    p_true: 20
    fc: 1.5
    fc_distribution: "fixed"  # "fixed", "uniform", "decreasing"
    correlation_source: "empirical"  # "empirical", "block", "ar1"
    missing_rate: 0.05
    missing_mechanism: "MNAR"
    preprocessing: "log_auto"
    n_confounders: 0
    n_interactions: 0

  # Scenari (variazioni dal base, uno alla volta)
  scenarios:
    # S1: p/n ratio
    - name: "pn_ratio"
      vary: "p_n"
      levels:
        - {n_case: 50, n_control: 50, p: 500}   # p/n = 5
        - {n_case: 25, n_control: 25, p: 500}    # p/n = 10
        - {n_case: 25, n_control: 25, p: 1000}   # p/n = 20
        - {n_case: 10, n_control: 10, p: 1000}   # p/n = 50

    # S2: Effect size
    - name: "effect_size"
      vary: "fc"
      levels: [1.2, 1.5, 2.0, 3.0]

    # S3: Multicollinearità
    - name: "correlation"
      vary: "correlation_scale"
      levels: [0.3, 0.5, 0.7, 0.9]
      note: "Scala la matrice di correlazione empirica"

    # S4: Non-linearità
    - name: "interactions"
      vary: "n_interactions"
      levels: [0, 5, 10]

    # S5: Confounders
    - name: "confounders"
      vary: "n_confounders"
      levels: [0, 1, 3]

    # S6: Missing data
    - name: "missing"
      vary: "missing_rate"
      levels: [0.0, 0.05, 0.15, 0.30]

    # S7: Preprocessing
    - name: "preprocessing"
      vary: "preprocessing"
      levels: ["none", "log_auto", "log_pareto", "pqn_auto"]

  # Parametri resampling
  resampling:
    n_replications: 30      # Dataset simulati per scenario
    n_bootstrap: 100         # Bootstrap resamples per stabilità
    subsample_fraction: 0.632  # Subsampling senza replacement
    stratified: true

# ----------------------------------------------------------
# METRICHE
# ----------------------------------------------------------
metrics:
  stability:
    primary: ["nogueira", "jaccard"]
    secondary: ["kuncheva", "dice", "spearman_rank"]

  accuracy:
    primary: ["tpr", "fdr"]
    secondary: ["f1", "mcc"]

  prediction:
    primary: ["auc"]
    secondary: ["balanced_accuracy"]
    classifier: "logistic_regression"  # Classificatore per valutazione predittiva
    test_fraction: 0.3

  parsimony:
    - "n_selected"
    - "prop_selected"

  stably_selected_thresholds: [0.50, 0.80, 0.95]

# ----------------------------------------------------------
# FIGURE
# ----------------------------------------------------------
figures:
  format: "pdf"
  dpi: 300
  width: 7       # inches
  height: 5
  color_palette: "viridis"
  font_family: "Helvetica"
  font_size: 10
```

---

## Script 00 — Installazione Pacchetti

### `R/00_install_packages.R`

**Funzione:** Installa tutti i pacchetti R necessari con verifica versione.
**Input:** Nessuno
**Output:** Tutti i pacchetti installati e verificati
**Tempo stimato:** 5-15 minuti (prima volta)

```r
# Lista completa dei pacchetti necessari:

cran_packages <- c(
  # Core
  "yaml", "here", "fs", "cli", "glue", "checkmate",
  # Data manipulation
  "data.table", "dplyr", "tidyr", "purrr", "tibble",
  # Feature selection
  "glmnet", "Boruta", "ranger", "stabs", "knockoff",
  "horseshoe", "spikeslab", "xgboost",
  # Stability metrics
  "stabm",
  # Simulation
  "mvtnorm", "corpcor", "MASS",
  # Statistics
  "pROC", "caret",
  # Visualization
  "ggplot2", "ggpubr", "ComplexHeatmap", "UpSetR",
  "patchwork", "scales", "RColorBrewer", "viridis",
  # Imputation
  "missForest",
  # Distribution fitting
  "fitdistrplus",
  # Parallelization
  "future", "furrr", "future.apply", "progressr",
  # SHAP
  "shapviz",
  # HTTP/API
  "httr2", "jsonlite",
  # Reporting
  "knitr", "rmarkdown",
  # Misc
  "microbenchmark"
)

bioc_packages <- c(
  "metabolomicsWorkbenchR",
  "SummarizedExperiment",
  "impute"  # Per KNN imputation
)

github_packages <- c(
  "aberHRML/metabolighteR",
  "ModelOriented/treeshap"
)
```

Lo script deve:
1. Verificare se ciascun pacchetto è già installato
2. Installare solo quelli mancanti
3. Caricare tutti e stampare versioni
4. Segnalare errori chiaramente

---

## Script 01 — Download Dati

### `R/01_download_data.R`

**Funzione:** Scarica dataset metabolomici da MetaboLights e Metabolomics Workbench.
**Input:** `config/config.yaml` (accession numbers e API endpoints)
**Output:** File in `data/raw/{source}/{accession}/`
**Dipendenze:** Pacchetti `metabolighteR`, `metabolomicsWorkbenchR`, `httr2`, `jsonlite`
**Tempo stimato:** 10-30 minuti (dipende dalla rete)

**Requisiti:**
1. Per ogni accession number in config:
   a. Verificare che lo studio esista e sia pubblico (API query)
   b. Scaricare la feature table (matrice metaboliti × campioni)
   c. Scaricare i metadata (sample information con covariate)
   d. Salvare in `data/raw/{source}/{accession}/`
   e. Loggare successo/fallimento

2. **Error handling robusto:**
   - Timeout configurabile (default 120s per richiesta)
   - Retry con backoff esponenziale (3 tentativi)
   - Se un dataset non è disponibile: warning (non errore fatale), skip e continua
   - Salvare un file `data/raw/download_log.csv` con stato di ogni download

3. **Per MetaboLights** (`metabolighteR`):
   ```r
   library(metabolighteR)
   # Per ogni accession:
   # 1. Ottenere lista file dello studio
   # 2. Identificare file della feature table (tipicamente *_maf.tsv o processed data)
   # 3. Scaricare feature table e sample metadata (ISA files)
   ```

4. **Per Metabolomics Workbench** (`metabolomicsWorkbenchR`):
   ```r
   library(metabolomicsWorkbenchR)
   # Per ogni accession:
   result <- do_query(
     context = "study",
     input_item = "study_id",
     input_value = accession,
     output_item = "data"
   )
   # Anche: output_item = "factors" per metadata campioni
   # Anche: output_item = "summary" per info studio
   ```

5. **Fallback per MetaboLights** (se API non funziona):
   - Provare download diretto via FTP: `ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{MTBLS_ID}/`
   - Se anche FTP fallisce: loggare e procedere con i dataset disponibili

6. **Verifica post-download:**
   - Feature table ha almeno 30 righe (campioni) e 50 colonne (features)?
   - Metadata contiene colonna gruppo (caso/controllo)?
   - File non corrotto (può essere letto come tabella)?

---

## Script 02 — Preprocessing e Armonizzazione

### `R/02_preprocess.R`

**Funzione:** Standardizza tutti i dataset scaricati in un formato comune.
**Input:** `data/raw/` (output di Script 01)
**Output:** `data/processed/{accession}.rds` — lista con componenti:
  - `X`: matrice numerica (campioni × features), log-trasformata
  - `y`: factor con livelli "case" e "control"
  - `metadata`: data.frame con covariate (age, sex, bmi se disponibili)
  - `feature_names`: character vector
  - `source`: "MetaboLights" o "MetabolomicsWorkbench"
  - `platform`: "NMR" o "LC-MS"
  - `accession`: ID originale
**Dipendenze:** Script 01
**Tempo stimato:** 2-5 minuti

**Pipeline per ogni dataset:**

1. **Identificazione formato e parsing:**
   - MetaboLights: file TSV con formati variabili; parsare header e identificare colonne metaboliti vs metadata
   - Metabolomics Workbench: output di `do_query()` già strutturato; convertire in matrice

2. **Pulizia features:**
   - Rimuovere features con >50% di missing values
   - Rimuovere features con varianza zero o quasi-zero (sd < 1e-10)
   - Rimuovere campioni con >40% di missing values

3. **Armonizzazione gruppi:**
   - Identificare la colonna gruppo nei metadata
   - Mappare a "case" / "control" (i nomi variano tra studi)
   - Verificare: almeno 20 campioni per gruppo dopo filtering

4. **Gestione missing values (fase iniziale):**
   - Sostituire zeri con NA (in metabolomica, zero spesso = non rilevato)
   - Imputazione con half-minimum per feature: `NA → min(feature) / 2`

5. **Log-trasformazione:**
   ```r
   X_log <- log2(X + 1)
   ```

6. **Quality check:**
   - PCA per identificare outlier evidenti
   - Distribuzione dei valori per feature (dovrebbe essere approssimativamente normale dopo log)
   - Report sintetico per ogni dataset

---

## Script 03 — Estrazione Parametri Empirici

### `R/03_extract_empirical_params.R`

**Funzione:** Stima parametri empirici dai dataset reali per guidare la simulazione.
**Input:** `data/processed/*.rds` (output di Script 02)
**Output:** `data/empirical_params/{accession}_params.rds` e `data/empirical_params/summary.rds`
**Dipendenze:** Script 02
**Tempo stimato:** 5-10 minuti

**Per ogni dataset processed, estrarre:**

1. **Distribuzione marginale per feature:**
   ```r
   library(fitdistrplus)
   # Per ogni feature j (su dati log-trasformati):
   fit_norm <- fitdist(X_log[, j], "norm")
   fit_lnorm <- fitdist(X[, j], "lnorm")  # Su dati originali
   # Scegliere il fit migliore (AIC)
   # Salvare: mean, sd per ogni feature
   ```

2. **Struttura di correlazione:**
   ```r
   library(corpcor)
   # Matrice di correlazione con shrinkage (Ledoit-Wolf):
   cor_shrink <- cor.shrink(X_log, verbose = FALSE)
   # Salvare la matrice completa
   ```

3. **Statistiche riassuntive della correlazione:**
   ```r
   # Distribuzione delle correlazioni pairwise:
   cor_vals <- cor_shrink[upper.tri(cor_shrink)]
   summary_cor <- list(
     mean = mean(abs(cor_vals)),
     median = median(abs(cor_vals)),
     q75 = quantile(abs(cor_vals), 0.75),
     q90 = quantile(abs(cor_vals), 0.90),
     max = max(abs(cor_vals)),
     prop_above_0.5 = mean(abs(cor_vals) > 0.5),
     prop_above_0.7 = mean(abs(cor_vals) > 0.7)
   )
   ```

4. **Missing values pattern:**
   ```r
   missing_info <- list(
     overall_rate = mean(is.na(X_original)),  # Prima dell'imputazione
     per_feature_rate = colMeans(is.na(X_original)),
     per_sample_rate = rowMeans(is.na(X_original))
   )
   ```

5. **Effect size empirici (caso vs controllo):**
   ```r
   fc <- colMeans(X[y == "case", ]) / colMeans(X[y == "control", ])
   log2fc <- log2(fc)
   cohen_d <- apply(X_log, 2, function(x) {
     (mean(x[y == "case"]) - mean(x[y == "control"])) /
       sqrt((var(x[y == "case"]) + var(x[y == "control"])) / 2)
   })
   ```

6. **Dimensionalità:**
   ```r
   dim_info <- list(
     n = nrow(X), p = ncol(X), p_n_ratio = ncol(X) / nrow(X)
   )
   ```

7. **Summary aggregato** (`summary.rds`):
   - Mediana dei parametri attraverso tutti i dataset
   - Range (min-max) per ogni parametro
   - Differenze NMR vs LC-MS

---

## Script 04 — Simulazione

### `R/04_simulate.R`

**Funzione:** Genera dataset simulati per tutti gli scenari definiti in config.
**Input:** `config/config.yaml` + `data/empirical_params/summary.rds`
**Output:** `data/simulated/{scenario}/{rep_XXX}/` con `X.rds`, `y.rds`, `true_features.rds`, `params.rds`
**Dipendenze:** Script 03
**Tempo stimato:** 30-60 minuti

**Implementazione della funzione principale:**

```r
simulate_dataset <- function(params, empirical_cor, seed) {
  set.seed(seed)

  n <- params$n_case + params$n_control
  p <- params$p
  p_true <- params$p_true

  # 1. Preparare matrice di correlazione
  if (params$correlation_source == "empirical") {
    # Usare la matrice empirica, ridimensionata a p
    Sigma <- prepare_correlation_matrix(empirical_cor, p, params$correlation_scale)
  } else if (params$correlation_source == "block") {
    Sigma <- generate_block_correlation(p, block_size = 50, within_cor = 0.7)
  }

  # 2. Generare dati multivariati normali (in scala log)
  Z <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma)

  # 3. Definire features "vere" (con segnale)
  true_idx <- sample(1:p, p_true)

  # 4. Impiantare segnale nelle features vere per il gruppo caso
  case_idx <- 1:params$n_case
  control_idx <- (params$n_case + 1):n

  delta <- generate_effect_sizes(p_true, params$fc, params$fc_distribution)
  Z[case_idx, true_idx] <- Z[case_idx, true_idx] + matrix(
    rep(delta, each = params$n_case), nrow = params$n_case
  )

  # 5. Aggiungere interazioni (scenario S4)
  if (params$n_interactions > 0) {
    Z <- add_interactions(Z, case_idx, true_idx, params$n_interactions)
  }

  # 6. Aggiungere confounders (scenario S5)
  confounders <- NULL
  if (params$n_confounders > 0) {
    result <- add_confounders(Z, case_idx, control_idx, params$n_confounders)
    Z <- result$Z
    confounders <- result$confounders
  }

  # 7. Trasformare in scala originale (log-normale)
  X <- exp(Z)

  # 8. Aggiungere missing values (MNAR)
  if (params$missing_rate > 0) {
    X <- add_mnar_missing(X, params$missing_rate)
  }

  # 9. Imputazione
  X <- impute_half_min(X)

  # 10. Preprocessing
  X <- apply_preprocessing(X, params$preprocessing)

  # 11. Creare outcome
  y <- factor(c(rep("case", params$n_case), rep("control", params$n_control)))

  return(list(
    X = X, y = y, true_features = true_idx,
    params = params, confounders = confounders
  ))
}
```

**Funzioni helper da implementare in `R/utils/helpers.R`:**

- `prepare_correlation_matrix(empirical_cor, p, scale)`: Ridimensiona/scala la matrice empirica
- `generate_block_correlation(p, block_size, within_cor)`: Genera struttura block-diagonal
- `generate_effect_sizes(p_true, fc, distribution)`: Genera vettore di effect sizes
- `add_interactions(Z, case_idx, true_idx, n_interactions)`: Aggiunge termini di interazione
- `add_confounders(Z, case_idx, control_idx, n_confounders)`: Aggiunge effetto confounders
- `add_mnar_missing(X, rate)`: Aggiunge missing MNAR basati su LOD
- `impute_half_min(X)`: Imputazione con half-minimum
- `apply_preprocessing(X, method)`: Applica il preprocessing specificato

---

## Script 05 — Feature Selection

### `R/05_feature_selection.R`

**Funzione:** Applica tutti i 12 metodi a tutti i dataset simulati con bootstrap resampling.
**Input:** `config/config.yaml` + `data/simulated/`
**Output:** `results/feature_selection/{scenario}/{rep_XXX}/selections.rds`
  Formato: lista con nomi metodo → matrice B × p (TRUE/FALSE per ogni feature × bootstrap)
**Dipendenze:** Script 04
**Tempo stimato:** 3-12 giorni (vedi docs/06_computational_plan.md)

**Struttura principale:**

```r
run_feature_selection <- function(scenario_dir, config) {
  B <- config$simulation$resampling$n_bootstrap
  frac <- config$simulation$resampling$subsample_fraction
  methods <- get_all_methods(config)

  # Caricare dataset
  X <- readRDS(file.path(scenario_dir, "X.rds"))
  y <- readRDS(file.path(scenario_dir, "y.rds"))
  n <- nrow(X)
  p <- ncol(X)

  # Matrice risultati: metodo → matrice B × p
  selections <- setNames(
    lapply(methods, function(m) matrix(FALSE, nrow = B, ncol = p)),
    sapply(methods, function(m) m$name)
  )

  for (b in 1:B) {
    # Subsampling stratificato
    idx <- stratified_subsample(y, fraction = frac, seed = b)
    X_b <- X[idx, ]
    y_b <- y[idx]

    # Applicare ogni metodo
    for (m in methods) {
      sel <- tryCatch(
        apply_method(m, X_b, y_b),
        error = function(e) {
          warning(paste("Method", m$name, "failed on bootstrap", b, ":", e$message))
          integer(0)
        }
      )
      selections[[m$name]][b, sel] <- TRUE
    }
  }

  return(selections)
}
```

**I 12 wrapper dei metodi devono essere in `R/utils/fs_methods.R`:**

Ogni wrapper deve:
1. Accettare `(X, y, params)` come input
2. Restituire un vettore di **indici** delle features selezionate
3. Gestire errori internamente (restituire `integer(0)` se fallisce)
4. Essere deterministico dato il seed

```r
# Esempio di struttura per ogni metodo:

fs_wilcoxon_fdr <- function(X, y, params) {
  pvals <- apply(X, 2, function(x) wilcox.test(x ~ y)$p.value)
  adj_pvals <- p.adjust(pvals, method = "BH")
  which(adj_pvals < params$fdr_threshold)
}

fs_lasso <- function(X, y, params) {
  fit <- cv.glmnet(X, y, family = "binomial", alpha = params$alpha, nfolds = params$nfolds)
  coefs <- as.vector(coef(fit, s = params$lambda_rule))[-1]
  which(coefs != 0)
}

fs_boruta <- function(X, y, params) {
  bor <- Boruta(X, y, maxRuns = params$maxRuns, pValue = params$pValue, doTrace = 0)
  bor <- TentativeRoughFix(bor)
  which(bor$finalDecision == "Confirmed")
}

fs_rf_importance <- function(X, y, params) {
  rf <- ranger(y ~ ., data = data.frame(X, y = y),
               importance = params$importance, num.trees = params$num_trees)
  imp <- rf$variable.importance
  if (params$top_k_method == "elbow") {
    k <- find_elbow(sort(imp, decreasing = TRUE))
  } else {
    k <- params$top_n
  }
  order(imp, decreasing = TRUE)[1:k]
}

fs_stability_selection <- function(X, y, params) {
  stab <- stabsel(X, y, fitfun = glmnet.lasso, cutoff = params$cutoff,
                  PFER = params$PFER, sampling.type = params$sampling_type,
                  B = params$B)
  which(names(stab$selected) %in% colnames(X))
}

fs_knockoff <- function(X, y, params) {
  result <- knockoff.filter(X, as.numeric(y) - 1, fdr = params$fdr)
  result$selected
}

fs_horseshoe <- function(X, y, params) {
  fit <- horseshoe(as.numeric(y) - 1, X,
                   method.tau = params$method_tau,
                   method.sigma = params$method_sigma,
                   burn = params$burn, nmc = params$nmc)
  # Selezione via credible intervals che non includono zero:
  which(fit$LeftCI * fit$RightCI > 0)
}

fs_spike_slab <- function(X, y, params) {
  fit <- spikeslab(as.numeric(y) - 1 ~ X, bigp.smalln = params$bigp_smalln,
                   max.var = params$max_var)
  which(fit$gnet.scale != 0)
}

fs_shap_xgboost <- function(X, y, params) {
  y_num <- as.numeric(y) - 1
  dtrain <- xgb.DMatrix(data = X, label = y_num)
  model <- xgb.train(
    params = list(objective = "binary:logistic", max_depth = params$max_depth, eta = params$eta),
    data = dtrain, nrounds = params$nrounds, verbose = 0
  )
  shp <- shapviz(model, X_pred = X)
  imp <- colMeans(abs(shp$S))
  if (params$top_k_method == "elbow") {
    k <- find_elbow(sort(imp, decreasing = TRUE))
  } else {
    k <- params$top_n
  }
  order(imp, decreasing = TRUE)[1:k]
}
```

**Parallelizzazione:**

```r
library(future)
library(furrr)
plan(multisession, workers = config$project$n_cores)

# Parallelizzare su scenario × replicazione:
all_tasks <- expand.grid(
  scenario = list.dirs("data/simulated", recursive = FALSE),
  rep = list.dirs(..., recursive = FALSE)
)

results <- future_map(1:nrow(all_tasks), function(i) {
  task <- all_tasks[i, ]
  # Verificare checkpoint
  if (checkpoint_exists(task)) return(NULL)
  # Eseguire
  res <- run_feature_selection(task$path, config)
  # Salvare checkpoint
  save_checkpoint(res, task)
  return(res)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)
```

---

## Script 06 — Metriche

### `R/06_metrics.R`

**Funzione:** Calcola tutte le metriche dai risultati della feature selection.
**Input:** `results/feature_selection/` + `data/simulated/` (per ground truth)
**Output:** `results/metrics/{scenario}_metrics.rds` e `results/metrics/all_metrics.rds`
**Dipendenze:** Script 05
**Tempo stimato:** 30-60 minuti

**Le metriche di stabilità devono essere in `R/utils/stability_metrics.R`:**

```r
library(stabm)

compute_all_metrics <- function(selections, true_features, X, y, config) {
  p <- ncol(X)
  methods <- names(selections)

  results <- lapply(methods, function(m) {
    sel_matrix <- selections[[m]]  # B × p matrice TRUE/FALSE
    B <- nrow(sel_matrix)

    # Convertire in lista di indici per stabm:
    feature_lists <- lapply(1:B, function(b) which(sel_matrix[b, ]))

    # --- STABILITÀ ---
    stab_nogueira <- stabilityNogueira(feature_lists, p = p)
    stab_jaccard <- stabilityJaccard(feature_lists, p = p)
    stab_kuncheva <- tryCatch(
      stabilityKuncheva(feature_lists, p = p),
      error = function(e) NA  # Fallisce se cardinalità diverse
    )
    stab_dice <- stabilityDice(feature_lists, p = p)

    # Spearman correlation of rankings (custom):
    # Costruire ranking da frequenza di selezione
    freq <- colSums(sel_matrix)
    # ... (dettaglio in stability_metrics.R)

    # Proportion stably selected:
    freq_prop <- freq / B
    prop_50 <- sum(freq_prop > 0.50) / p
    prop_80 <- sum(freq_prop > 0.80) / p
    prop_95 <- sum(freq_prop > 0.95) / p

    # --- ACCURATEZZA (con ground truth) ---
    # Usare la selezione "aggregata" (features in >50% dei bootstrap):
    consensus <- which(freq_prop > 0.50)
    tpr <- length(intersect(consensus, true_features)) / length(true_features)
    fdr <- if (length(consensus) > 0) {
      1 - length(intersect(consensus, true_features)) / length(consensus)
    } else NA

    precision <- 1 - fdr
    f1 <- if (!is.na(precision) && !is.na(tpr) && (precision + tpr) > 0) {
      2 * precision * tpr / (precision + tpr)
    } else NA

    # --- PERFORMANCE PREDITTIVA ---
    # AUC su test set (media su bootstrap)
    aucs <- sapply(1:min(B, 50), function(b) {  # Max 50 per velocità
      sel <- which(sel_matrix[b, ])
      if (length(sel) < 2) return(NA)
      compute_auc_on_test(X, y, sel, seed = b, test_frac = config$metrics$prediction$test_fraction)
    })
    auc_mean <- mean(aucs, na.rm = TRUE)
    auc_sd <- sd(aucs, na.rm = TRUE)

    # --- PARSIMONIA ---
    n_selected_median <- median(rowSums(sel_matrix))
    n_selected_mean <- mean(rowSums(sel_matrix))

    return(data.frame(
      method = m,
      nogueira = stab_nogueira,
      jaccard = stab_jaccard,
      kuncheva = stab_kuncheva,
      dice = stab_dice,
      prop_stable_50 = prop_50,
      prop_stable_80 = prop_80,
      prop_stable_95 = prop_95,
      tpr = tpr,
      fdr = fdr,
      f1 = f1,
      auc_mean = auc_mean,
      auc_sd = auc_sd,
      n_selected_median = n_selected_median,
      n_selected_mean = n_selected_mean,
      stringsAsFactors = FALSE
    ))
  })

  do.call(rbind, results)
}

# Helper per AUC:
compute_auc_on_test <- function(X, y, selected_features, seed, test_frac) {
  set.seed(seed)
  n <- nrow(X)
  test_idx <- sample(n, round(n * test_frac))
  train_idx <- setdiff(1:n, test_idx)

  X_train <- X[train_idx, selected_features, drop = FALSE]
  X_test <- X[test_idx, selected_features, drop = FALSE]
  y_train <- y[train_idx]
  y_test <- y[test_idx]

  # Logistic regression
  df_train <- data.frame(y = y_train, X_train)
  fit <- tryCatch(
    glm(y ~ ., data = df_train, family = binomial),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(fit)) return(NA)

  df_test <- data.frame(X_test)
  pred <- predict(fit, newdata = df_test, type = "response")
  pROC::auc(pROC::roc(y_test, pred, quiet = TRUE))
}
```

---

## Script 07 — Validazione Cross-Database

### `R/07_cross_validation.R`

**Funzione:** Testa la replicabilità delle features selezionate su dataset reali indipendenti della stessa patologia.
**Input:** `data/processed/` (dataset reali) + config
**Output:** `results/cross_validation/` con metriche di replicabilità
**Dipendenze:** Script 02 (dataset processed), Script 05 (per i metodi, opzionale)
**Tempo stimato:** 1-4 ore

**Protocollo:**

1. Identificare coppie di dataset della stessa patologia (es. due studi su diabete T2)
2. Per ogni coppia (D1, D2):
   a. Applicare tutti i 12 metodi a D1 → set selezionato S1
   b. Applicare tutti i 12 metodi a D2 → set selezionato S2
   c. Calcolare overlap: Jaccard(S1, S2)
   d. Calcolare direzione: quante features hanno lo stesso segno di fold-change in entrambi i dataset?
   e. Cross-prediction: trainare su D1 con features S1, predire D2 → AUC

3. **Metriche di replicabilità:**
   - Jaccard index tra S1 e S2 (overlap diretto)
   - Concordanza direzionale: % features con stesso segno FC
   - Cross-prediction AUC
   - Per ogni metodo: media e CI su tutte le coppie disponibili

4. **Nota:** Il numero di coppie disponibili dipende dai dataset effettivamente scaricati. Se nessuna coppia è disponibile per la stessa patologia, usare un approccio "leave-one-dataset-out" su dataset di patologie diverse come proof-of-concept.

---

## Script 08 — Figure e Tabelle

### `R/08_figures.R`

**Funzione:** Genera tutte le figure per il paper scientifico.
**Input:** `results/metrics/all_metrics.rds` + `results/cross_validation/`
**Output:** `results/figures/` (PDF, 300 DPI)
**Dipendenze:** Script 06, Script 07
**Tempo stimato:** 5-15 minuti

**Figure da generare:**

### Figure 1: Heatmap di Stabilità (metodi × scenari)
```r
# Matrice: righe = 12 metodi, colonne = scenari
# Colore = Nogueira stability index (0 a 1)
# Annotazioni: valore numerico in ogni cella
# Clustering gerarchico su righe e colonne
# Palette: viridis (giallo = alta stabilità, viola = bassa)
```

### Figure 2: TPR vs FDR per Scenario
```r
# Pannello multiplo (facet_wrap per scenario)
# Ogni punto = un metodo
# X = FDR, Y = TPR
# Dimensione punto = n_selected (parsimonia)
# Colore = categoria metodo (filter, embedded, wrapper, bayesian, ML)
# Ideale: angolo in alto a sinistra (alto TPR, basso FDR)
```

### Figure 3: Radar/Spider Charts per Confronto Metodi
```r
# Un radar chart per ogni metodo (o top-6 metodi)
# Assi: Nogueira, TPR, 1-FDR, AUC, 1/n_selected (normalizzati 0-1)
# Sovrapporre metodi per confronto diretto
# Alternativa: fmsb::radarchart()
```

### Figure 4: UpSet Plot per Overlap tra Metodi
```r
library(UpSetR)
# Per lo scenario base:
# Mostrare quante features sono selezionate da 1, 2, ..., 12 metodi
# Identificare il "core" di features selezionate da molti metodi
```

### Figure 5: Stabilità vs Predizione Trade-off
```r
# Scatterplot: X = Nogueira, Y = AUC
# Ogni punto = metodo × scenario
# Colore = metodo
# Forma = tipo scenario
# Cercare: c'è un trade-off? O alcuni metodi dominano?
```

### Figure 6: Effetto p/n Ratio sulla Stabilità
```r
# Line plot: X = p/n ratio, Y = Nogueira
# Una linea per metodo
# Mostrare CI come ribbon
# Evidenziare il "breakdown point" di ogni metodo
```

### Figure 7: Cross-Database Replicabilità
```r
# Barplot: metodo vs Jaccard (overlap tra dataset indipendenti)
# Barre di errore per CI
# Ordinare per Jaccard decrescente
```

### Supplementary Figure S1: Distribuzione Empirica delle Correlazioni
```r
# Istogramma delle correlazioni pairwise dai dataset reali
# Separato per NMR e LC-MS
# Sovrapposto con la distribuzione usata nella simulazione
```

### Tabelle:

**Table 1:** Caratteristiche dei dataset pubblici usati (accession, platform, n, p, patologia)
**Table 2:** Riepilogo dei 12 metodi (nome, tipo, pacchetto R, parametri chiave)
**Table 3:** Metriche aggregate per lo scenario base (tutti i metodi)
**Table S1:** Metriche complete per tutti gli scenari (supplementary)

---

## Note Implementative Trasversali

### Error Handling

Ogni script deve:
- Usare `tryCatch` per operazioni che possono fallire
- Loggare messaggi con `cli::cli_alert_info()`, `cli::cli_alert_warning()`, `cli::cli_alert_danger()`
- Salvare checkpoint frequenti (vedi Script 05)
- Non interrompere l'esecuzione per un singolo fallimento (skip e continua)

### Riproducibilità

- Ogni operazione stocastica deve usare `set.seed()` con seed derivato dalla configurazione
- Salvare `sessionInfo()` in `results/session_info.txt`
- Tutti i percorsi relativi tramite `here::here()`
- Config letta una volta all'inizio e passata come argomento

### Logging

```r
# Usare cli per output strutturato:
cli::cli_h1("Script 05 — Feature Selection")
cli::cli_alert_info("Scenario: {scenario_name}, Replicazione: {rep_id}")
cli::cli_progress_bar("Bootstrap resampling", total = B)
```

### Struttura di ogni Script

```r
#!/usr/bin/env Rscript
# ==============================================================================
# Script XX — [Nome]
# Progetto: Feature Selection Stability in Metabolomics
# Autore: [nome]
# Data: 2026-03
# ==============================================================================

# 0. Setup
library(yaml)
library(here)
config <- yaml::read_yaml(here("config", "config.yaml"))
set.seed(config$project$seed)

# 1. [Primo step]
# ...

# 2. [Secondo step]
# ...

# N. Salvataggio risultati
# ...

cli::cli_alert_success("Script XX completato.")
sessionInfo()
```

---

## Riepilogo Dipendenze tra Script

```
00_install_packages.R    (standalone, eseguire per primo)
        ↓
01_download_data.R       (richiede connessione internet)
        ↓
02_preprocess.R          (input: data/raw/)
        ↓
03_extract_empirical_params.R  (input: data/processed/)
        ↓
04_simulate.R            (input: data/empirical_params/ + config)
        ↓
05_feature_selection.R   (input: data/simulated/ — IL PIÙ COSTOSO)
        ↓
06_metrics.R             (input: results/feature_selection/ + data/simulated/)
        ↓
08_figures.R             (input: results/metrics/)

07_cross_validation.R    (input: data/processed/ — INDIPENDENTE da 04-06)
        ↓
08_figures.R             (input anche da results/cross_validation/)
```

---

## Stima Tempo Totale

| Script | Tempo Stimato |
|--------|--------------|
| 00 | 5-15 min |
| 01 | 10-30 min |
| 02 | 2-5 min |
| 03 | 5-10 min |
| 04 | 30-60 min |
| **05** | **3-12 giorni** (dipende da hardware) |
| 06 | 30-60 min |
| 07 | 1-4 ore |
| 08 | 5-15 min |
| **Totale** | **~3-13 giorni** (dominato da Script 05) |
