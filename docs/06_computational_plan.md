# Piano Computazionale

> Ultimo aggiornamento: 2026-03-17
> Status: DRAFT

---

## 1. Stime di Tempo per Metodo (singolo fit, p=500, n=100)

Stime conservative basate su benchmark tipici. Tempo per singola applicazione del metodo.

| Metodo | Tempo stimato | Note |
|--------|--------------|------|
| Wilcoxon + FDR | ~0.01s | Molto veloce, p test indipendenti |
| Fold-change | ~0.001s | Operazione vettoriale |
| Volcano | ~0.01s | Combina FC + p-value |
| LASSO (cv.glmnet) | ~0.5s | Con 10-fold CV |
| Elastic Net (cv.glmnet) | ~0.5s | Simile a LASSO |
| Boruta | ~10s | 500 iterazioni × RF |
| RF importance (ranger) | ~1s | 1000 alberi, ranger è veloce |
| Stability selection (stabs) | ~30s | B=100 subsample × LASSO |
| Knockoff filter | ~2s | Stima covarianza + knockoff generation |
| Horseshoe (`horseshoe`) | ~60s | 5000 MCMC iterations, burn-in 1000 |
| Spike-and-slab (`spikeslab`) | ~30s | Con bigp.smalln=TRUE |
| SHAP (XGBoost + shapviz) | ~5s | Train XGBoost + calcolo SHAP |

**Tempo totale per singolo dataset (12 metodi):** ~140s ≈ 2.3 minuti

---

## 2. Stima Totale

### Configurazione completa (dal simulation design)

| Parametro | Valore |
|-----------|--------|
| Scenari | 28 condizioni |
| Replicazioni/scenario | 50 |
| Bootstrap/replicazione | 200 |
| Metodi | 12 |

### Calcolo

```
Totale fit = 28 × 50 × 200 × 12 = 3,360,000 fit

Tempo medio per fit ≈ 140s / 12 ≈ 11.7s (media pesata)

Tempo sequenziale totale ≈ 3,360,000 × 11.7s ≈ 39,312,000s ≈ 455 giorni
```

**→ Sequenziale: impossibile. Serve parallelizzazione aggressiva e/o riduzione dello scope.**

---

## 3. Strategia di Riduzione

### Opzione A: Ridurre parametri

| Parametro | Originale | Ridotto | Fattore |
|-----------|-----------|---------|---------|
| Replicazioni | 50 | 20 | 2.5× |
| Bootstrap | 200 | 100 | 2× |
| Scenari | 28 | 15 (focus sui più informativi) | 1.87× |

**Totale ridotto:** 15 × 20 × 100 × 12 = 360,000 fit
**Tempo sequenziale:** ~4,212,000s ≈ 48.7 giorni

### Opzione B: Parallelizzazione

Con 8 core (laptop moderno):
```
48.7 giorni / 8 core ≈ 6.1 giorni
```

Con 32 core (workstation o cloud):
```
48.7 giorni / 32 core ≈ 1.5 giorni
```

Con 64 core (HPC cluster):
```
48.7 giorni / 64 core ≈ 0.76 giorni ≈ 18 ore
```

### Opzione C (Raccomandata): Riduzione + Parallelizzazione

**Parametri finali proposti:**

| Parametro | Valore | Giustificazione |
|-----------|--------|-----------------|
| Scenari | 20 | Eliminare combinazioni ridondanti |
| Replicazioni | 30 | Sufficiente per CI robusti |
| Bootstrap | 100 | Nogueira (2018) suggerisce 100 come sufficiente |
| Metodi | 12 | Tutti |

**Totale:** 20 × 30 × 100 × 12 = 720,000 fit

**Con 8 core (laptop M1/M2):** ~720,000 × 11.7s / 8 ≈ 1,053,000s ≈ **12.2 giorni**
**Con 16 core:** ≈ **6.1 giorni**
**Con 32 core (cloud):** ≈ **3 giorni**

---

## 4. Parallelizzazione in R

### Framework raccomandato: `future` + `furrr`

```r
library(future)
library(furrr)

# Setup parallelismo:
plan(multisession, workers = parallel::detectCores() - 1)

# Parallelizzare il loop esterno (scenario × replicazione):
results <- future_map(scenario_rep_grid, function(row) {
  scenario <- row$scenario
  rep_id <- row$rep_id
  seed <- row$seed

  # Genera dataset
  data <- simulate_dataset(scenario, seed)

  # Per ogni bootstrap:
  boot_results <- map(1:B, function(b) {
    boot_data <- subsample(data, seed = seed + b)
    # Applica tutti i metodi:
    selections <- apply_all_methods(boot_data)
    return(selections)
  })

  # Calcola metriche
  metrics <- compute_metrics(boot_results, data$true_features)
  return(metrics)
}, .options = furrr_options(seed = TRUE))
```

### Alternativa: `parallel::mclapply` (più semplice, solo Unix)

```r
library(parallel)
n_cores <- detectCores() - 1

results <- mclapply(1:nrow(grid), function(i) {
  # ... come sopra
}, mc.cores = n_cores)
```

### Approccio gerarchico

1. **Livello esterno** (parallelizzato): scenario × replicazione
2. **Livello medio** (sequenziale per semplicità): bootstrap runs
3. **Livello interno** (sequenziale): metodi (alcuni usano RNG e devono essere riproducibili)

Il livello esterno ha 20 × 30 = 600 unità → buon granularity per 8-32 core.

---

## 5. Metodi Computazionalmente Critici

### 5.1 Horseshoe Prior (MCMC)

- **Tempo stimato:** ~60s per fit con p=500, n=100
- **Con p=1000:** ~240s (scala ~p²)
- **Totale per horseshoe solo:** 720,000/12 × 60s = 3,600,000s ≈ 41.7 giorni sequenziali
- **GPU:** Non necessaria per `horseshoe` package (no Stan/GPU); `brms` può usare CmdStan con threading
- **Mitigazione:** Ridurre MCMC iterations a 2000 (burn=500, nmc=1500); usare `horseshoe` package (più veloce di brms per questo task specifico)

### 5.2 Spike-and-Slab

- **Tempo stimato:** ~30s per fit con p=500
- **`spikeslab` con bigp.smalln=TRUE:** Ottimizzato per p >> n, usa approccio gnet (non MCMC completo)
- **Mitigazione:** OK per le dimensioni proposte

### 5.3 Boruta

- **Tempo stimato:** ~10s per fit con p=500
- **Con p=1000:** ~20s
- **Mitigazione:** Ridurre maxRuns a 200 (da 500 default); usare `ranger` come backend (più veloce di `randomForest`)

### 5.4 Stability Selection (stabs)

- **Tempo stimato:** ~30s per fit (include B=100 subsample del LASSO interno)
- **Nota:** Stability selection CON bootstrap resampling esterno = B_esterno × B_interno = doppio sampling
- **Mitigazione:** Per stability selection, il bootstrap esterno potrebbe essere ridotto (il metodo ha il suo resampling interno)

---

## 6. Stima Risorse

### Hardware minimo (laptop)

- CPU: 8 core (M1/M2 Mac o equivalente)
- RAM: 16 GB (sufficiente per p=1000, n=200 per singolo fit)
- Storage: ~50 GB per risultati completi (RDS compressi)
- Tempo: ~12 giorni

### Hardware raccomandato (workstation/cloud)

- CPU: 32 core
- RAM: 64 GB
- Storage: ~50 GB
- Tempo: ~3 giorni
- Costo cloud (es. AWS c5.8xlarge 32 vCPU): ~$1.36/ora × 72 ore ≈ $100

### Hardware ideale (HPC)

- CPU: 64+ core
- RAM: 128 GB
- Tempo: ~1.5 giorni

---

## 7. Strategia di Checkpointing

```r
# Dopo ogni scenario × replicazione completata, salvare risultato:
save_checkpoint <- function(result, scenario, rep_id) {
  dir_path <- file.path("results", "checkpoints", scenario)
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(dir_path, paste0("rep_", sprintf("%03d", rep_id), ".rds")))
}

# All'inizio di ogni run, verificare se il checkpoint esiste:
checkpoint_exists <- function(scenario, rep_id) {
  file.exists(file.path("results", "checkpoints", scenario,
                         paste0("rep_", sprintf("%03d", rep_id), ".rds")))
}

# → Permette di interrompere e riprendere la computazione
```

---

## 8. Memory Footprint

| Oggetto | Dimensione stimata |
|---------|-------------------|
| Matrice X (n=100, p=1000) | ~0.8 MB |
| 200 bootstrap di X | ~160 MB |
| Risultati 12 metodi × 200 bootstrap | ~50 MB (solo indici selezionati) |
| Risultati aggregati per 1 scenario × 30 rep | ~1.5 GB |
| **Totale risultati (20 scenari)** | **~30 GB** |

RAM necessaria per singolo worker: ~1 GB (dataset + modelli temporanei). Con 8 workers: ~8 GB.

---

## 9. Profiling e Ottimizzazione

### Step 1: Micro-benchmark (prima della run completa)

```r
library(microbenchmark)

# Testare ogni metodo su un singolo dataset:
mb <- microbenchmark(
  wilcoxon = apply_wilcoxon(X, y),
  lasso = apply_lasso(X, y),
  boruta = apply_boruta(X, y),
  # ... etc
  times = 5
)
print(mb)
```

### Step 2: Pilot run

Eseguire 1 scenario × 5 replicazioni × 50 bootstrap per stimare tempo reale e verificare correttezza.

### Step 3: Full run

Solo dopo validazione del pilot.

---

## 10. Dipendenze Computazionali tra Script

```
Script 01 (Download) → Script 02 (Preprocessing) → Script 03 (Parametri empirici)
                                                          ↓
                                                    Script 04 (Simulazione)
                                                          ↓
                                                    Script 05 (Feature selection)
                                                          ↓
                                                    Script 06 (Metriche)
                                                          ↓
Script 07 (Cross-validation, indipendente da 04-06) ←→ Script 08 (Figure)
```

Script 07 usa dati REALI (non simulati) → può essere eseguito in parallelo con 04-05-06.
