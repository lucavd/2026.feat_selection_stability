# Stato Implementazione Pipeline

> Ultimo aggiornamento: 2026-03-18
> Status: IMPLEMENTATO — bugfix e hardening applicati post-deploy

---

## 1. File Implementati

| File | Status | Testato | Descrizione |
|------|--------|---------|-------------|
| `config/config.yaml` | ✅ | ✅ | Configurazione centrale unica |
| `R/00_install_packages.R` | ✅ | ✅ | 49 pacchetti (aggiunto `digest`); 3 da fonti alternative |
| `R/01_download_data.R` | ✅ | ✅ | Download da 3 fonti: CIMCB GitHub (Excel), MetaboLights (FTP/API), MW (REST). 9/10 scaricati |
| `R/02_preprocess.R` | ✅ | ✅ | Parser CIMCB (Excel), MetaboLights (MAF), MW (nested JSON). QC, imputazione mediana, normalizzazione. 9 dataset processati |
| `R/03_extract_empirical_params.R` | ✅ | ✅ | Correlazione Ledoit-Wolf, distribuzioni, eigenvalues. 9/9 dataset OK |
| `R/04_simulate.R` | ✅ | — | MVN → exp → FC + confounders + interazioni + MNAR missing. Supporta `correlation_source`: empirical, ar1, block |
| `R/05_feature_selection.R` | ✅ | — | 11 metodi × 100 bootstrap, checkpointing, parallelizzato per dataset (100 worker) |
| `R/06_metrics.R` | ✅ | — | Nogueira, Jaccard, TPR, FDR, AUC (split stratificato, feature priority), parsimonia. Output: `metrics_all.rds`, `metrics_by_scenario.rds`, `metrics_summary.rds` |
| `R/07_cross_validation.R` | ✅ | — | Bootstrap FS su dati reali + concordanza cross-database. Gestisce edge case 0/1 bootstrap convergenti |
| `R/08_figures.R` | ✅ | — | 8 figure PDF publication-quality + supplementari |
| `R/utils/helpers.R` | ✅ | ✅ | Config, logging, checkpoint, seed deterministico, parallelizzazione (cap 100 workers) |
| `R/utils/fs_methods.R` | ✅ | ✅ | 12 wrapper con interfaccia uniforme + dispatcher |
| `R/utils/stability_metrics.R` | ✅ | ✅ | Nogueira (Eq.4 + jackknife CI; B==2 → NA varianza/CI), Jaccard, Kuncheva, Dice, Spearman |

---

## 2. Risultati Esecuzione (2026-03-18)

### 2.1 Script 00 — Installazione pacchetti

- **49/49 pacchetti installati con successo** (aggiunto `digest`, `readxl`)
- Installazioni non standard: `horseshoe` (CRAN archive), `ComplexHeatmap` (Bioconductor), `treeshap` (GitHub), `knockoff` (richiede `Rdsdp` compilato manualmente)
- Tempo: ~36 secondi

### 2.2 Script 01 — Download dati

| Dataset | Source | Status | Note |
|---------|--------|--------|------|
| ST001047 | CIMCB GitHub | ✅ | 151 KB Excel — Gastric cancer NMR |
| MTBLS404 | CIMCB GitHub | ✅ | 223 KB Excel — Sacurine LC-MS |
| ST001000 | CIMCB GitHub | ✅ | 7.6 MB Excel — IBD LC-MS |
| MTBLS136 | CIMCB GitHub | ✅ | 11.6 MB Excel — Hormone use LC-MS |
| MTBLS92 | CIMCB GitHub | ✅ | 675 KB Excel — Breast cancer LC-MS |
| ST000369 | CIMCB GitHub | ✅ | 206 KB Excel — Lung cancer GC-MS |
| ST000496 | CIMCB GitHub | ✅ | 92 KB Excel — Periodontal GC-MS |
| MTBLS28 | MetaboLights | ✅ | 5 files (2 MAF NEG/POS, sample sheet, 2 assays) |
| MTBLS374 | MetaboLights | ✅ | 5 files scaricati ma MAF con solo 3 features |
| ST001706 | MW REST API | ✅ | 3 files (data JSON, factors JSON, metabolites) |

**10/10 scaricati. MTBLS374 escluso in fase di preprocessing (MAF non contiene dati numerici).**

### 2.3 Script 02 — Preprocessing

| Dataset | n | p | Classi | Missing % |
|---------|--:|---:|--------|-----------|
| ST001047 | 83 | 149 | GC(43)/HE(40) | 5.6% |
| MTBLS404 | 184 | 120 | M(101)/F(83) | 0% |
| ST001000 | 107 | 1533 | UC(58)/CD(49) | 20.2% |
| MTBLS136 | 668 | 787 | E-only(331)/E+P(337) | 6% |
| MTBLS92 | 253 | 138 | Pre(142)/Post(111) | 0% |
| ST000369 | 80 | 181 | Cancer(49)/Ctrl(31) | 0.1% |
| ST000496 | 100 | 69 | Pre(50)/Post(50) | 0% |
| MTBLS28 | 1005 | 1359 | Control(536)/Case(469) | 0% |
| ST001706 | 256 | 50 | Control(174)/RCC(82) | 0% |

**9/10 processati con successo. Tempo totale: ~5 secondi.**

### 2.4 Script 03 — Estrazione parametri empirici

| Dataset | Platform | |log2FC| med. | |cor| medio | Rank eff. | Top 10 EV | Missing raw |
|---------|----------|-------------|------------|-----------|-----------|-------------|
| ST001047 | NMR | 0.31 | 0.283 | 31 | 57.3% | 5.6% |
| MTBLS404 | LC-MS | 0.35 | 0.364 | 20 | 69.4% | 0% |
| ST001000 | LC-MS | 0.40 | 0.086 | 57 | 27.8% | 20.2% |
| MTBLS136 | LC-MS | 0.16 | 0.072 | 155 | 27.9% | 6% |
| MTBLS92 | LC-MS | 0.48 | 0.278 | 20 | 75.7% | 0% |
| ST000369 | GC-MS | 0.33 | 0.068 | 37 | 28.3% | 0.1% |
| ST000496 | GC-MS | 0.43 | 0.310 | 12 | 74.3% | 0% |
| MTBLS28 | LC-MS | 0.21 | 0.117 | 183 | 41.6% | 0% |
| ST001706 | NMR | 0.58 | 0.111 | 15 | 56.0% | 0% |

**Dataset di riferimento per simulazione (per piattaforma, il più grande):**
- NMR → ST001047 (p=149)
- LC-MS → ST001000 (p=1533)
- GC-MS → ST000369 (p=181)

**9/9 parametri estratti. Tempo totale: ~4 secondi.**

### 2.5 Script 04 — Simulazione

- **780 dataset MVN** (S1-S7): 7 scenari × 3-4 livelli × 30 rep. Tempo: ~76s con 8 worker
- **810 dataset semi-sintetici** (S8): 9 dataset reali × 3 FC (1.2, 1.5, 2.0) × 30 rep. Tempo: ~45s con 8 worker
- **Totale: 1590 dataset simulati**

### 2.6 Tempi misurati

- Script 00: ~36 secondi
- Script 01: ~21 secondi (CIMCB download veloce, MetaboLights FTP OK, MW API OK)
- Script 02: ~5 secondi
- Script 03: ~4 secondi
- Script 04: ~2 minuti (780 MVN + 810 semi-sintetici)
- Script 05: in corso (stima ~6 ore con 100 worker)

---

## 3. Scelte Implementative vs Design Originale

### 3.1 Parametri ridotti (Opzione C del piano computazionale)

| Parametro | Design originale (doc 05) | Implementato (config.yaml) | Motivazione |
|-----------|--------------------------|---------------------------|-------------|
| Replicazioni (R) | 50 | 30 | Sufficiente per CI robusti; riduce tempo 40% |
| Bootstrap (B) | 200 | 100 | Nogueira (2018) suggerisce 100 sufficiente |
| Scenari × livelli | 28 | 25 (7 scenari × 3-4 livelli) | Tutte le combinazioni in config |
| min_samples_per_group | 30 | 20 | Ampliare pool dataset disponibili |
| Horseshoe | incluso | **rimosso** | MCMC ~60 min/job a p=500; insostenibile (vedi §5.3) |
| Boruta maxRuns | 500 | 100 | Default del pacchetto (Kursa & Rudnicki 2010); convergenza verificata |
| n_cores (workers) | 8 | 100 | 128 core disponibili, 28 riservati al sistema |

**Stima computazionale finale:** 1590 dataset × 11 metodi × 100 boot = **1,749,000 fit**
(780 MVN S1-S7 + 810 semi-sintetici S8; horseshoe rimosso, da 12 a 11 metodi)

### 3.2 Bugfix e hardening (2026-03-18)

| Script | Fix | Impatto |
|--------|-----|---------|
| `00_install_packages.R` | Aggiunto pacchetto `digest` alle dipendenze CRAN | Necessario per hashing deterministico |
| `02_preprocess.R` | Fix sintassi data.table: `..data_cols` → `data_cols, with = FALSE` | Evita errore parsing MAF in certi ambienti |
| `02_preprocess.R` | Nuove helper `get_matching_colnames()` e `align_sample_info()` per matching robusto sample ↔ metadata | Gestisce ID column variabili tra dataset (Sample Name, Subject ID, etc.) |
| `02_preprocess.R` | QC filters ora aggiornano `feature_info` e `sample_info` coerentemente con `X` e `y` | Evita mismatch dimensioni dopo filtering |
| `02_preprocess.R` | Output ora include `X_raw` (pre-impute) e `X_imputed` (post-impute) separati | Permette analisi impatto imputazione |
| `04_simulate.R` | Nuovo parametro `correlation_source` ("empirical", "ar1", "block") | Supporta strutture di correlazione alternative per sensitivity analysis |
| `06_metrics.R` | `compute_prediction_metrics()`: split stratificato train/test, `feature_priority` per pruning intelligente | Evita test set senza una classe; usa importanza per selezionare features quando troppe |
| `06_metrics.R` | Guard `length(unique(y)) < 2` | Evita crash su dataset degenerati |
| `06_metrics.R` | Output aggiuntivo `metrics_by_scenario.rds` | Facilita analisi per scenario senza riaggregare |
| `07_cross_validation.R` | Gestione edge case: 0 o 1 bootstrap convergenti → metriche stabilità = NA | Evita crash di `compute_all_stability()` con input insufficiente |
| `07_cross_validation.R` | `na.rm = TRUE` nelle soglie di selezione frequenza | Evita NA propagation nel conteggio features stabili |
| `utils/helpers.R` | `setup_parallel()` cap workers a 100 | Previene esaurimento risorse su macchine con molti core |
| `utils/stability_metrics.R` | Nogueira index con B==2: restituisce NA per varianza e CI | Jackknife leave-one-out richiede B≥3; evita divisione per zero |
| `02_preprocess.R` | `pmax(X, 1)` → `log2(X + 1)` in `normalize_matrix()` per tutti i metodi (log_auto, log_pareto, pqn_auto, fallback) | `pmax(X, 1)` distruggeva feature con valori sub-unitari (es. M86 in MTBLS92: range [0.09, 0.63] → costante → NaN dopo `scale()`) |
| `03_extract_empirical_params.R` | log2FC calcolato come differenza `(mean_case - mean_control) / log(2)` anziché `log2(ratio)` | Dati su scala log: il FC è una differenza, non un rapporto. Prima dava valori ~33 (nonsense) |
| `03_extract_empirical_params.R` | Rimozione feature zero-varianza/all-NA prima della stima correlazione; tracking `cor_features_kept` e dimensione `n_features_cor` | Evita `eigen()` con NA/Inf; documenta mismatch dimensionale cor_matrix vs p per Script 04 |
| `03_extract_empirical_params.R` | `na.rm = TRUE` in `colMeans()` e `sd()` per medie/sd marginali; NA handling in skewness/kurtosis | Gestisce residui NA in X post-preprocessing |
| `04_simulate.R` | Aggiunto scenario S8 semi-sintetico: spike-in su dati reali con `simulate_semisynthetic()` | 810 dataset aggiuntivi; correlazione e distribuzioni biologicamente autentiche |
| `05_feature_selection.R` | Parallelizzazione ristrutturata: da chunk dataset×metodo a per-dataset (ogni worker = 1 dataset × 11 metodi sequenziali) | Elimina idle time; worker sempre occupati |
| `utils/fs_methods.R` | `num.threads = 1` in Boruta e rf_importance (ranger) | Evita contention tra worker paralleli |
| `utils/fs_methods.R` | `fs_shap_xgboost`: xgboost GPU (`device = "cuda"`) + `predict(predcontrib = TRUE)` al posto di `treeshap::treeshap()` | 250x speedup (da >150 min a 35s per 100 boot) |
| `utils/fs_methods.R` | Cache device GPU in `.xgb_device` (global env) | Evita overhead test GPU ripetuto per ogni chiamata |
| `config.yaml` | horseshoe rimosso, maxRuns Boruta 300→100, n_cores 8→100 | Fattibilità computazionale |

### 3.3 Imputazione

- **Design originale:** half-minimum, KNN, missForest, MinProb
- **Implementato:** median imputation come default
- **Motivazione:** L'imputation method non è il focus; median è conservativo e riproducibile. Sensitivity analysis con altri metodi può essere aggiunta come supplementary material.

### 3.4 Preprocessing

Implementate tutte e 4 le opzioni previste:
| Codice config | Operazione |
|---------------|------------|
| `none` | Nessun preprocessing |
| `log_auto` | log2(x+1) + autoscaling (mean-center, unit variance) — **DEFAULT** |
| `log_pareto` | log2(x+1) + Pareto scaling (÷ √sd) |
| `pqn_auto` | PQN + log2(x+1) + autoscaling |

### 3.5 Segnale nelle simulazioni

- Dati generati come MVN → exp() (distribuzione log-normale)
- FC applicato moltiplicativamente post-exp ai campioni caso
- Direzioni up/down randomizzate per realismo
- **Nota equivalenza:** FC moltiplicativo su scala originale = shift additivo in scala log. La struttura di correlazione è preservata perché il segnale è impiantato dopo la trasformazione.

### 3.7 Sorgenti dati (2026-03-18)

- **Design originale:** 6 MetaboLights + 3 MW, download diretto da API
- **Implementato:** 7 CIMCB benchmark (Excel pre-processati) + 1 MetaboLights (MAF) + 1 MW (REST JSON)
- **Motivazione:** Le API MetaboLights e MW restituiscono formati variabili e spesso incompleti. Il repository CIMCB di Mendez et al. (2019) fornisce dati validati in benchmark pubblicati. MTBLS28 (il più grande) scaricato direttamente da MetaboLights MAF. ST001706 scaricato via MW REST API per coprire NMR.
- **Riferimento:** Mendez, Reinke & Broadhurst (2019), *Metabolomics* 15:150

### 3.8 Metriche di stabilità — Nogueira index

- Implementazione diretta Equazione 4 di Nogueira et al. (2018, JMLR 18(174):1-54)
- Varianza: jackknife leave-one-out (Teorema 7)
- CI: φ̂ ± z_{α/2} × √(Var_hat), default α=0.05

---

## 4. Architettura del Codice

### 4.1 Flusso dati completo

```
config/config.yaml ─── [letto da ogni script all'avvio]
          │
          ▼
R/00 ─── Verifica pacchetti
          │
          ▼
R/01 ─── data/raw/<accession>/ ─── download_manifest.rds
          │
          ▼
R/02 ─── data/processed/<acc>_processed.rds ─── datasets_summary.rds
          │                                │
          │                                └───► R/07 (path indipendente)
          ▼                                          │
R/03 ─── data/empirical_params/empirical_params.rds  │
          │                                          │
          ▼                                          ▼
R/04 ─── data/simulated/sim_<scen>_<lev>_rep<r>.rds  results/cross_validation/
          │                                              cv_<acc>_<method>.rds
          ▼                                              cv_summary.rds
R/05 ─── results/feature_selection/                      concordance_results.rds
          fs_<task>_<method>.rds                          │
          │                                               │
          ▼                                               │
R/06 ─── results/metrics/                                 │
          metrics_all.rds (master table)                  │
          metrics_by_scenario.rds                         │
          metrics_summary.rds                             │
          │                                               │
          ▼                                               │
R/08 ─── results/figures/ ◄───────────────────────────────┘
          fig_01..08.pdf + supplementari
```

### 4.2 Convenzioni codice

| Aspetto | Convenzione |
|---------|-------------|
| Path | `here::here()` per tutti i percorsi |
| Config | `yaml::read_yaml(here("config", "config.yaml"))` |
| Logging | `cli::cli_alert_info/success/warning/danger()` |
| Errors | `tryCatch()` — mai interrompere per singolo fallimento |
| Seed | `derive_seed(base_seed, id)` — hash xxhash32, deterministico |
| Checkpoint | `.rds` file; script verificano esistenza prima di ricalcolare |
| Parallelo | `future::multisession` + `furrr::future_map()` |
| Session | `sessionInfo()` salvata in `results/session_info_XX.txt` |

### 4.3 Interfaccia fs_methods.R

Ogni wrapper restituisce una lista con struttura fissa:
```r
list(
  selected   = logical(p),     # TRUE = feature selezionata
  importance = numeric(p),     # Score importanza (NA se non applicabile)
  time       = numeric(1),     # Secondi
  converged  = logical(1),     # TRUE se successo
  message    = character(1)    # Errore/warning
)
```

Il dispatcher `run_fs_method(name, X, y, params)` mappa il nome (stringa) al wrapper corretto.

---

## 5. Problemi Noti e Workaround

### 5.1 Dataset download

| Problema | Workaround |
|----------|------------|
| MetaboLights FTP lento/down | Fallback su ISA-Tab bundle via API REST |
| MTBLS97 restituisce 403 | Potrebbe essere studio ad accesso limitato; sostituire se necessario |
| MW REST API ha rate limit | Retry con backoff esponenziale (3 tentativi, 5s × attempt) |
| File MAF con formato variabile | Parsing robusto: identifica colonne dati vs metadata per tipo numerico |

### 5.2 Knockoff filter (n > p richiesto)

- Fallisce silenziosamente negli scenari S1c (p/n=20) e S1d (p/n=50)
- `converged = FALSE`, `message = "n <= p: knockoff requires n > p"`
- Nelle metriche aggregate, questi run sono esclusi (non falsano i risultati)

### 5.3 Horseshoe — RIMOSSO

- **Problema:** ~55-60 min per job (100 bootstrap × MCMC burn=1000 + nmc=3000) con p=500, n=50. Con p=1000 stimato ~4h/job.
- **Pilot misurato:** su scenario S1_pn10 (p=500, n=50), horseshoe impiegava 55 min/job vs ~8 min per boruta (il secondo più lento). ~85% del tempo totale.
- **Proiezione:** 1590 dataset × 55 min = ~60 giorni di runtime solo per horseshoe (8 core). Insostenibile.
- **Decisione:** Rimosso dalla pipeline. Spike-slab (~1-2 min/job) copre la categoria bayesiana. Horseshoe citato nel paper come "tested in pilot, excluded for computational infeasibility at scale."
- **Pipeline ora ha 11 metodi** (era 12)

### 5.4 Boruta — maxRuns ridotto a 100

- maxRuns=100 è il default del pacchetto (Kursa & Rudnicki 2010, JMLR 11:271-282)
- Early stopping: se tutte le feature sono Confirmed/Rejected, termina prima (~60-80 iter su p=500)
- `num.threads = 1` in ranger per evitare contention tra worker paralleli
- TentativeRoughFix per features ancora indecise a convergenza

### 5.5 Stability selection — doppio sampling

- `stabs` ha un resampling interno (B=100 subsample del LASSO)
- Il nostro bootstrap esterno (B=100) crea doppio sampling
- Questo è intenzionale: vogliamo misurare la stabilità di stability selection, non usarla come stability selection intende

---

## 5b. Ottimizzazione Computazionale — Cronologia (2026-03-18)

### Problema originale

Il primo run di Script 05 (con 12 metodi, chunk_size=16, parallelizzazione per dataset×metodo) produceva ~0.48 job/min. Proiezione: **27.6 giorni** per 19.080 job.

Analisi dei timestamp sui file completati ha rivelato la distribuzione dei tempi per singolo job (100 bootstrap, n=100, p=500):

| Metodo | Tempo/job | % del totale |
|--------|----------|-------------|
| fold_change, knockoff | <5s | ~0% |
| elastic_net, rf_imp, wilcoxon, volcano, lasso | 13-21s | ~3% |
| spike_slab | 5 min | ~12% |
| stability_selection | 10 min | ~24% |
| boruta (maxRuns=300) | 22 min | ~52% |
| **horseshoe** | **55-60 min** | **~85%** (quando presente nel chunk) |
| **shap_xgboost (treeshap CPU)** | **>150 min** | killer |

### Intervento 1: Rimozione horseshoe

- Horseshoe MCMC (burn=1000, nmc=3000) è intrinsecamente lento e non parallelizzabile su GPU
- Rimosso dalla pipeline; spike-slab copre la categoria bayesiana
- **Impatto:** da 12 a 11 metodi. Proiezione scesa a ~17 giorni

### Intervento 2: XGBoost GPU + SHAP nativo (shap_xgboost)

- **Problema:** `treeshap::treeshap()` su CPU era il vero killer (>150 min per 100 bootstrap)
- **Soluzione:** Installato xgboost 3.2.0 precompilato con supporto CUDA (RTX 4090) + sostituzione di treeshap con `predict(fit, predcontrib = TRUE)` (SHAP nativo)
- **Risultato misurato:** da >150 min a **35 secondi** per 100 bootstrap (~250x speedup)
- File: xgboost 3.2.0 GPU binary da `s3://xgboost-nightly-builds/release_3.2.0/`
- Device detection cached in `.xgb_device` (global env) per evitare overhead ripetuto

### Intervento 3: Tentativi per Boruta (non implementati)

Tre approcci testati per accelerare Boruta, nessuno dei quali ha funzionato:

1. **Boruta XGBoost GPU (reimplementazione algoritmo in R):** L'algoritmo Boruta è stato reimplementato usando XGBoost GPU come base learner al posto di RF. Risultato: più lento (72s vs 24s singola chiamata) e troppo conservativo (3 feature selezionate vs 19). L'importanza SHAP si distribuisce diversamente dalla permutation importance di RF, rendendo la soglia shadow inefficace. **Scartato.**

2. **Boruta Python (sklearn via system()):** Wrapper R→Python con BorutaPy + sklearn RandomForestClassifier. sklearn è 2.5x più veloce di R ranger per singola chiamata, ma l'overhead I/O (npy save/load + process spawn) per ogni bootstrap annullava il vantaggio: 17 min vs 22 min. **Marginale.**

3. **Boruta Python batch (tutto il bootstrap in un processo Python):** Script `boruta_batch.py` che esegue tutti i 100 bootstrap in un unico processo Python. Problema: 27 min (peggio di R). Motivo: `n_jobs=1` nel batch sequenziale; con `n_jobs=-1` non parallelizzava (confermato dal monitoring CPU). Il loop di Boruta è intrinsecamente sequenziale (ogni iterazione dipende dalla precedente). **Scartato.**

4. **cuML GPU RandomForest via BorutaPy:** Installato RAPIDS cuML (venv `/home/user/gpu_env`). Risultato: cuML RF su GPU era **più lento** di sklearn CPU (44.5s vs 9.5s) per n=100, p=500. I dati sono troppo piccoli per saturare la GPU — l'overhead di trasferimento CPU↔GPU domina. **Scartato.**

**Conclusione Boruta:** R Boruta con `num.threads=1` e `maxRuns=100` (default pacchetto) è l'opzione più veloce: ~5.6 min per 100 bootstrap su n=100, p=500. La riduzione di maxRuns da 300 a 100 è il fattore più significativo (~4x speedup).

### Intervento 4: Ristrutturazione parallelizzazione Script 05

- **Prima:** parallelizzazione per dataset×metodo. Chunk di 16 task misti → worker veloci idle aspettando boruta/stab_sel
- **Dopo:** parallelizzazione per dataset. Ogni worker prende 1 dataset ed esegue tutti 11 metodi sequenzialmente. Nessun idle time, nessuna contention
- `num.threads = 1` in Boruta e rf_importance per evitare contention tra worker
- n_cores aumentato da 8 a 100 (128 core disponibili, 28 riservati)

### Intervento 5: Scenario semi-sintetico (S8)

- Aggiunto scenario S8_semisynthetic: spike-in su 9 dataset reali × 3 FC × 30 rep = 810 dataset
- Usa campioni controllo reali divisi in pseudo-caso/pseudo-controllo
- Segnale (FC) impiantato su p_true feature random
- Correlazione, distribuzioni, missing autenticamente biologici
- Totale dataset simulati: 780 MVN + 810 semi-sintetici = 1590

### Benchmark finale (n=100, p=500, 100 bootstrap stimati)

| Metodo | Tempo est. 100 boot | Note |
|--------|---------------------|------|
| fold_change | 0s | |
| knockoff | 1s | |
| elastic_net | 12s | |
| rf_importance | 16s | num.threads=1 |
| wilcoxon_fdr | 18s | |
| volcano | 18s | |
| lasso | 18s | |
| shap_xgboost | 35s | GPU + SHAP nativo |
| spike_slab | 4.1 min | |
| boruta | 5.6 min | maxRuns=100, num.threads=1 |
| stability_selection | 10.3 min | |
| **Totale per dataset** | **~22 min** | |

**Stima runtime finale:** 1590 dataset / 100 worker × 22 min = **~5.8 ore**

---

## 6. Prossimi Passi — Checklist Operativa

### Completati
- [x] Tutti gli script scritti e verificati sintatticamente
- [x] Config.yaml completo con tutti i parametri
- [x] Pacchetti installati (49/49 OK su R 4.5.2)
- [x] Script 00 eseguito — 49 pacchetti OK
- [x] Script 01 eseguito — 10/10 dataset scaricati (7 CIMCB + 1 MetaboLights + 1 MW + 1 MTBLS374 inutilizzabile)
- [x] Script 02 eseguito — 9/9 dataset processati con successo
- [x] Script 03 eseguito — 9/9 parametri estratti; 3 piattaforme (NMR, LC-MS, GC-MS)
- [x] Script 04 eseguito — 1590 dataset simulati (780 MVN + 810 semi-sintetici)
- [x] Ottimizzazione Script 05 — horseshoe rimosso, xgboost GPU, Boruta maxRuns=100, parallelizzazione per dataset, 100 worker

### Prossimi passi
1. [ ] Eseguire Script 05 (`bash run_fs_with_telegram.sh`)
3. [ ] **Pilot run script 05:** 1 scenario × 3 rep × 20 bootstrap
4. [ ] Verificare tempi reali e risultati del pilot
5. [ ] **Full run:** eseguire 05 → 06 (sequenziali)
6. [ ] In parallelo: eseguire `Rscript R/07_cross_validation.R`
7. [ ] Eseguire `Rscript R/08_figures.R`

### Post-esecuzione
- [ ] Verificare `metrics_all.rds`, `metrics_by_scenario.rds`, `metrics_summary.rds`
- [ ] Verificare figure in `results/figures/`
- [ ] Analisi esplorativa dei risultati
- [ ] Stesura paper

---

## 7. Dataset Alternativi (in caso di fallimento download)

Se alcuni accession nel config non sono scaricabili, sostituire con:

### MetaboLights — NMR
- `MTBLS19` — Breast cancer serum NMR (n~100)
- `MTBLS124` — Inflammatory bowel disease NMR (n~80)

### MetaboLights — LC-MS
- `MTBLS146` — Ovarian cancer plasma LC-MS (n~200)
- `MTBLS263` — Hepatocellular carcinoma serum LC-MS (n~100)

### Metabolomics Workbench
- `ST000284` — Colorectal cancer serum (n~80)
- `ST001000` — Breast cancer plasma (n~60)

**NOTA:** Questi accession sono suggerimenti da verificare interrogando le API sul server. Non sono stati validati.
