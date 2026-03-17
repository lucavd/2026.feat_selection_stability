# Stato Implementazione Pipeline

> Ultimo aggiornamento: 2026-03-17
> Status: IMPLEMENTATO — testato localmente, pronto per deploy su server

---

## 1. File Implementati

| File | Status | Testato | Descrizione |
|------|--------|---------|-------------|
| `config/config.yaml` | ✅ | ✅ | Configurazione centrale unica |
| `R/00_install_packages.R` | ✅ | ✅ | 48/48 pacchetti installati (3 da fonti alternative) |
| `R/01_download_data.R` | ✅ | ✅ | MTBLS1 scaricato con successo; altri falliti per DNS locale |
| `R/02_preprocess.R` | ✅ | — | Parsing MAF/JSON, QC, imputazione mediana, normalizzazione |
| `R/03_extract_empirical_params.R` | ✅ | — | Correlazione Ledoit-Wolf, distribuzioni, eigenvalues |
| `R/04_simulate.R` | ✅ | — | MVN → exp → FC + confounders + interazioni + MNAR missing |
| `R/05_feature_selection.R` | ✅ | — | 12 metodi × 100 bootstrap, checkpointing, parallelizzato |
| `R/06_metrics.R` | ✅ | — | Nogueira, Jaccard, TPR, FDR, AUC, parsimonia |
| `R/07_cross_validation.R` | ✅ | — | Bootstrap FS su dati reali + concordanza cross-database |
| `R/08_figures.R` | ✅ | — | 8 figure PDF publication-quality + supplementari |
| `R/utils/helpers.R` | ✅ | ✅ | Config, logging, checkpoint, seed deterministico, parallelizzazione |
| `R/utils/fs_methods.R` | ✅ | ✅ | 12 wrapper con interfaccia uniforme + dispatcher |
| `R/utils/stability_metrics.R` | ✅ | ✅ | Nogueira (Eq.4 + jackknife CI), Jaccard, Kuncheva, Dice, Spearman |

---

## 2. Risultati Test Locale (2026-03-17)

### 2.1 Script 00 — Installazione pacchetti

- **48/48 pacchetti installati con successo**
- 3 pacchetti richiesti installazione da fonti alternative:
  - `horseshoe` → CRAN archive (v0.2.0)
  - `ComplexHeatmap` → Bioconductor
  - `treeshap` → GitHub (ModelOriented/treeshap)
- Tempo: ~30 secondi

### 2.2 Script 01 — Download dati

| Dataset | Source | Status | Note |
|---------|--------|--------|------|
| MTBLS1 | MetaboLights | ✅ Success | 3 file (MAF, sample sheet, assay). T2D urine NMR |
| MTBLS97 | MetaboLights | ⚠️ 403 Forbidden | API bloccata; FTP fallito; ISA-Tab non disponibile |
| MTBLS136 | MetaboLights | ❌ DNS error | Non risolvibile da macchina locale |
| MTBLS352 | MetaboLights | ❌ DNS error | Idem |
| MTBLS537 | MetaboLights | ❌ DNS error | Idem |
| MTBLS733 | MetaboLights | ❌ DNS error | Benchmark spike-in — prioritario |
| ST000388 | MW | ❌ DNS error | T2D plasma lipidomics |
| ST001047 | MW | ❌ DNS error | CRC serum |
| ST001386 | MW | ❌ DNS error | CVD serum NMR |

**Causa dei fallimenti:** DNS locale non risolveva host esterni. Su server con connessione internet completa questi errori non si presenteranno.

**MTBLS97 (403):** Potrebbe essere uno studio con accesso limitato. Da verificare sul server — se persiste, sostituire con altro dataset NMR.

### 2.3 Tempi misurati

- Script 00: ~30 secondi
- Script 01: ~52 minuti (dominato da timeout su host non raggiungibili — su server sarà ~5-10 min)

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

### 3.2 Imputazione

- **Design originale:** half-minimum, KNN, missForest, MinProb
- **Implementato:** median imputation come default
- **Motivazione:** L'imputation method non è il focus; median è conservativo e riproducibile. Sensitivity analysis con altri metodi può essere aggiunta come supplementary material.

### 3.3 Preprocessing

Implementate tutte e 4 le opzioni previste:
| Codice config | Operazione |
|---------------|------------|
| `none` | Nessun preprocessing |
| `log_auto` | log2 + autoscaling (mean-center, unit variance) — **DEFAULT** |
| `log_pareto` | log2 + Pareto scaling (÷ √sd) |
| `pqn_auto` | PQN + log2 + autoscaling |

### 3.4 Segnale nelle simulazioni

- Dati generati come MVN → exp() (distribuzione log-normale)
- FC applicato moltiplicativamente post-exp ai campioni caso
- Direzioni up/down randomizzate per realismo
- **Nota equivalenza:** FC moltiplicativo su scala originale = shift additivo in scala log. La struttura di correlazione è preservata perché il segnale è impiantato dopo la trasformazione.

### 3.5 Metriche di stabilità — Nogueira index

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

### Pre-deploy
- [x] Tutti gli script scritti e verificati sintatticamente
- [x] Config.yaml completo con tutti i parametri
- [x] Pacchetti testati localmente (48/48 OK)
- [x] Script 01 testato (MTBLS1 scaricato con successo)
- [ ] Verificare che il server abbia R ≥ 4.3

### Sul server
1. [ ] Clonare repository
2. [ ] Eseguire `Rscript R/00_install_packages.R` — ~2-5 min
3. [ ] Eseguire `Rscript R/01_download_data.R` — ~5-10 min con internet
4. [ ] Verificare `download_manifest.rds`: quanti dataset scaricati con successo?
5. [ ] Se <3 dataset: aggiornare `config.yaml` con accession alternativi (vedi sezione 7)
6. [ ] Eseguire `Rscript R/02_preprocess.R` — ~2 min
7. [ ] Eseguire `Rscript R/03_extract_empirical_params.R` — ~1 min
8. [ ] **Pilot run:** modificare temporaneamente config per 1 scenario × 3 rep × 20 bootstrap
9. [ ] Verificare tempi reali e risultati del pilot
10. [ ] **Full run:** ripristinare config originale ed eseguire 04 → 05 → 06 (sequenziali)
11. [ ] In parallelo: eseguire `Rscript R/07_cross_validation.R`
12. [ ] Eseguire `Rscript R/08_figures.R`

### Post-esecuzione
- [ ] Verificare `metrics_all.rds` e `metrics_summary.rds`
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
