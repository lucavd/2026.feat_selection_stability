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
| `R/03_extract_empirical_params.R` | ✅ | — | Correlazione Ledoit-Wolf, distribuzioni, eigenvalues |
| `R/04_simulate.R` | ✅ | — | MVN → exp → FC + confounders + interazioni + MNAR missing. Supporta `correlation_source`: empirical, ar1, block |
| `R/05_feature_selection.R` | ✅ | — | 12 metodi × 100 bootstrap, checkpointing, parallelizzato |
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

### 2.4 Tempi misurati

- Script 00: ~36 secondi
- Script 01: ~21 secondi (CIMCB download veloce, MetaboLights FTP OK, MW API OK)
- Script 02: ~5 secondi

---

## 3. Scelte Implementative vs Design Originale

### 3.1 Parametri ridotti (Opzione C del piano computazionale)

| Parametro | Design originale (doc 05) | Implementato (config.yaml) | Motivazione |
|-----------|--------------------------|---------------------------|-------------|
| Replicazioni (R) | 50 | 30 | Sufficiente per CI robusti; riduce tempo 40% |
| Bootstrap (B) | 200 | 100 | Nogueira (2018) suggerisce 100 sufficiente |
| Scenari × livelli | 28 | 25 (7 scenari × 3-4 livelli) | Tutte le combinazioni in config |
| min_samples_per_group | 30 | 20 | Ampliare pool dataset disponibili |
| Horseshoe MCMC | burn=1000, nmc=5000 | burn=1000, nmc=3000 | Riduce tempo ~40% per fit |
| Boruta maxRuns | 500 | 300 | Bilancia accuratezza/tempo |

**Stima computazionale finale:** 25 livelli × 30 rep × 100 boot × 12 metodi = **900,000 fit**

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

### 3.3 Imputazione

- **Design originale:** half-minimum, KNN, missForest, MinProb
- **Implementato:** median imputation come default
- **Motivazione:** L'imputation method non è il focus; median è conservativo e riproducibile. Sensitivity analysis con altri metodi può essere aggiunta come supplementary material.

### 3.4 Preprocessing

Implementate tutte e 4 le opzioni previste:
| Codice config | Operazione |
|---------------|------------|
| `none` | Nessun preprocessing |
| `log_auto` | log2 + autoscaling (mean-center, unit variance) — **DEFAULT** |
| `log_pareto` | log2 + Pareto scaling (÷ √sd) |
| `pqn_auto` | PQN + log2 + autoscaling |

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

### 5.3 Horseshoe — performance

- **~60s/fit** con p=500, n=100; **~240s/fit** con p=1000
- Horseshoe da solo: ~900,000/12 × 60s = ~75,000 fit × 60s = 4,500,000s ≈ 52 giorni sequenziali
- Con 8 core: ~6.5 giorni solo per horseshoe
- **Mitigazione:** MCMC ridotto (burn=1000, nmc=3000); checkpointing per ripresa

### 5.4 Boruta — timeout possibile

- maxRuns=300 dovrebbe completare in ~10-15s per fit con p=500
- Su dataset con p=1000: ~25s/fit
- TentativeRoughFix per features indecise

### 5.5 Stability selection — doppio sampling

- `stabs` ha un resampling interno (B=100 subsample del LASSO)
- Il nostro bootstrap esterno (B=100) crea doppio sampling
- Questo è intenzionale: vogliamo misurare la stabilità di stability selection, non usarla come stability selection intende

---

## 6. Prossimi Passi — Checklist Operativa

### Completati
- [x] Tutti gli script scritti e verificati sintatticamente
- [x] Config.yaml completo con tutti i parametri
- [x] Pacchetti installati (49/49 OK su R 4.5.2)
- [x] Script 00 eseguito — 49 pacchetti OK
- [x] Script 01 eseguito — 10/10 dataset scaricati (7 CIMCB + 1 MetaboLights + 1 MW + 1 MTBLS374 inutilizzabile)
- [x] Script 02 eseguito — 9/9 dataset processati con successo

### Prossimi passi
1. [ ] Eseguire `Rscript R/03_extract_empirical_params.R`
2. [ ] Eseguire `Rscript R/04_simulate.R`
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
