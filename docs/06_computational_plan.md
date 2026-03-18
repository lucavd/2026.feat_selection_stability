# Piano Computazionale

> Ultimo aggiornamento: 2026-03-18 (rev. 2)
> Status: AGGIORNATO con benchmark reali e ottimizzazioni GPU

---

## 1. Tempi Misurati per Metodo (100 bootstrap, p=500, n=100)

Benchmark reali misurati su workstation con 128 core + NVIDIA RTX 4090.

| Metodo | Tempo 100 boot | Categoria | Note |
|--------|---------------|-----------|------|
| fold_change | <1s | filter | Operazione vettoriale |
| knockoff | 1s | embedded | Fallisce se n ≤ p (`converged = FALSE`) |
| elastic_net | 12s | embedded | cv.glmnet, 10-fold |
| rf_importance | 16s | wrapper | ranger, `num.threads = 1` |
| wilcoxon_fdr | 18s | filter | |
| volcano | 18s | filter | |
| lasso | 18s | embedded | cv.glmnet, lambda.1se |
| shap_xgboost | 35s | ml | **GPU (CUDA) + SHAP nativo** |
| spike_slab | 4.1 min | bayesian | `bigp.smalln = TRUE` |
| boruta | 5.6 min | wrapper | maxRuns=100, `num.threads = 1` |
| stability_selection | 10.3 min | meta | B=100 internal subsample × LASSO |
| **Totale per dataset** | **~22 min** | | |

### Metodi rimossi

| Metodo | Tempo misurato | Motivo rimozione |
|--------|---------------|------------------|
| horseshoe | ~55-60 min/job (p=500) | MCMC intrinsecamente lento; proiezione 60+ giorni per pipeline completa |
| shap_xgboost (treeshap CPU) | >150 min/job | Sostituito da GPU + SHAP nativo (~250x speedup) |

---

## 2. Configurazione Finale

| Parametro | Valore | Giustificazione |
|-----------|--------|-----------------|
| Dataset MVN (S1-S7) | 780 | 26 livelli × 30 rep |
| Dataset semi-sintetici (S8) | 810 | 9 dataset × 3 FC × 30 rep |
| **Totale dataset** | **1590** | |
| Bootstrap (B) | 100 | Nogueira (2018) suggerisce 100 sufficiente |
| Metodi | 11 | horseshoe rimosso |
| **Totale fit** | **1,749,000** | |

---

## 3. Parallelizzazione

### Strategia: per-dataset

Ogni worker prende 1 dataset ed esegue tutti 11 metodi sequenzialmente. Questo evita idle time (nessun worker aspetta boruta mentre altri finiscono metodi veloci).

```
Worker 1: dataset_001 → [fc, wilcox, volcano, lasso, enet, rf, knockoff, shap, spike, boruta, stabsel] → save 11 files
Worker 2: dataset_002 → [fc, wilcox, ...] → save 11 files
...
Worker 100: dataset_100 → [...]
```

### Configurazione

- **Workers:** 100 (`future::multisession`)
- **CPU:** 128 core disponibili, 28 riservati al sistema
- **GPU:** 1× RTX 4090 condivisa (usata da shap_xgboost, ~0.4s/fit → contention trascurabile)
- **Thread interni:** `num.threads = 1` in ranger/Boruta per evitare contention tra worker

### Stima runtime

```
1590 dataset / 100 worker = 16 batch
16 batch × 22 min/dataset = ~5.8 ore
```

Con overhead I/O e varianza tra scenari (p=1000 più lento): **stima conservativa ~6-8 ore**.

---

## 4. Checkpointing

Ogni dataset×metodo produce un file indipendente:
```
results/feature_selection/fs_<task_id>_<method>.rds
```

Lo script verifica `file.exists()` prima di ogni job → ripresa automatica dopo interruzione. I 1590 × 11 = 17.490 file attesi consentono monitoraggio granulare del progresso.

---

## 5. Risorse

### Hardware attuale

| Componente | Specifiche |
|-----------|-----------|
| CPU | 128 core |
| RAM | (disponibile) |
| GPU | NVIDIA RTX 4090 (24 GB VRAM) |
| Storage | ~50 GB per risultati completi |

### RAM per worker

~5-10 MB per dataset (X matrix + modelli temporanei). Con 100 worker: ~1 GB totale. Non è un collo di bottiglia.

### Dipendenze software non-standard

- `xgboost` 3.2.0 con supporto CUDA (prebuilt binary da GitHub releases)
- Python venv `/home/user/gpu_env` con cuML + sklearn (usato solo per test, non in produzione)

---

## 6. Dipendenze tra Script

```
Script 01 (Download) → Script 02 (Preprocessing) → Script 03 (Parametri empirici)
                                ↓                          ↓
                          Script 07 (CV reali)      Script 04 (Simulazione)
                                ↓                          ↓
                                │                    Script 05 (Feature selection) ← BOTTLENECK
                                │                          ↓
                                │                    Script 06 (Metriche)
                                │                          ↓
                                └──────────────────→ Script 08 (Figure)
```

Script 07 (cross-validation su dati reali) è indipendente da 04-06 e può essere eseguito in parallelo.
