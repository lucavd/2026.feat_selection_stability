# CLAUDE.md — Feature Selection Stability in Metabolomics

## Progetto

Articolo scientifico in silico: **"The stability illusion: A simulation-based assessment of feature selection methods for high-dimensional metabolomics data"**

Target journal: Briefings in Bioinformatics / Metabolomics / Analytical Chemistry.

## Regole fondamentali

1. **NON inventare dati, DOI, accession numbers, o risultati.** Ogni riferimento deve essere verificato.
2. **NON hardcodare parametri negli script.** Tutto viene da `config/config.yaml`.
3. **Documentare tutto.** Questo è un progetto scientifico — ogni decisione deve essere tracciabile.
4. **Aggiornare la documentazione in `docs/` ad ogni cambiamento significativo.**
5. **Seed riproducibili.** Ogni operazione stocastica usa seed derivato da `config$project$seed`.

## Struttura del progetto

```
config/config.yaml          — UNICO punto di configurazione
R/00_install_packages.R     — Installazione dipendenze
R/01_download_data.R        — Download dataset CIMCB + MetaboLights + MW
R/02_preprocess.R           — Preprocessing e armonizzazione
R/03_extract_empirical_params.R — Parametri empirici da dati reali
R/04_simulate.R             — Framework simulazione
R/05_feature_selection.R    — 12 metodi × bootstrap resampling (BOTTLENECK)
R/06_metrics.R              — Nogueira, TPR, FDR, AUC, parsimonia
R/07_cross_validation.R     — Validazione cross-database (dati reali)
R/08_figures.R              — Figure per il paper
R/utils/helpers.R           — Funzioni condivise
R/utils/fs_methods.R        — Wrapper per ogni metodo FS
R/utils/stability_metrics.R — Wrapper metriche stabilità
docs/01-06_*.md             — Documentazione ricerca dettagliata
```

## Dipendenze tra script

```
00 → 01 → 02 → 03 → 04 → 05 → 06 → 08
                 ↓                    ↑
                 07 ──────────────────┘
```

Script 07 (cross-validation su dati reali) è indipendente da 04-06 (simulazione).

## Convenzioni di codice R

- Usare `here::here()` per tutti i percorsi
- Leggere config con `yaml::read_yaml(here("config", "config.yaml"))`
- Logging con `cli::cli_alert_*()` (info, warning, danger, success)
- Error handling con `tryCatch()` — mai interrompere per un singolo fallimento
- Checkpoint dopo ogni unità di lavoro completata (scenario × replicazione)
- `sessionInfo()` salvata alla fine di ogni script

## Pacchetti chiave

| Scopo | Pacchetto | Fonte |
|-------|-----------|-------|
| LASSO/Elastic Net | `glmnet` | CRAN |
| Boruta | `Boruta` | CRAN |
| Random Forest | `ranger` | CRAN |
| Stability selection | `stabs` | CRAN |
| Knockoff filter | `knockoff` | CRAN |
| Horseshoe prior | `horseshoe` | CRAN |
| Spike-and-slab | `spikeslab` | CRAN |
| SHAP | `shapviz` + `treeshap` | CRAN/GitHub |
| Stabilità metriche | `stabm` | CRAN |
| Simulazione | `mvtnorm` + `corpcor` | CRAN |
| Parallelizzazione | `future` + `furrr` | CRAN |
| MetaboLights API | `metabolighteR` | CRAN |
| MW API | `metabolomicsWorkbenchR` | Bioconductor |
| CIMCB Excel | `readxl` | CRAN |

## Stime computazionali

- **Script 05 è il bottleneck:** ~3-12 giorni con 8 core
- Configurazione raccomandata: 20 scenari × 30 rep × 100 bootstrap × 12 metodi = 720,000 fit
- Usa checkpointing per ripresa dopo interruzioni
- Parallelizzazione su livello scenario × replicazione via `future::multisession`

## Metriche primarie

- **Stabilità:** Nogueira index con CI (Nogueira et al. 2018, JMLR)
- **Accuratezza:** TPR + FDR (ground truth nota in simulazione)
- **Predizione:** AUC su test set holdout
- **Parsimonia:** n features selezionate

## Note per l'autore

- L'autore è un biostatistico con esperienza in metabolomica clinica (NMR + LC-MS)
- Usa regolarmente Boruta
- Preferisce: precisione, step-by-step, fonti verificate, nessuna improvvisazione
- Comunicazione in italiano, codice e contenuto scientifico in inglese
