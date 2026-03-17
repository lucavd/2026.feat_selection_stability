# Guida Deploy su Server

> Ultimo aggiornamento: 2026-03-17
> Status: PRONTO PER DEPLOY

---

## 1. Requisiti Server

### Hardware minimo (laptop/workstation)
| Risorsa | Minimo | Raccomandato |
|---------|--------|--------------|
| CPU | 8 core | 16-32 core |
| RAM | 16 GB | 64 GB |
| Storage | 50 GB liberi | 100 GB |
| Internet | Necessaria per Script 01 | — |

### Hardware ideale (HPC/cloud)
| Risorsa | Valore | Costo stimato (AWS) |
|---------|--------|---------------------|
| CPU | 32-64 core | c5.8xlarge ($1.36/h) o c5.16xlarge ($2.72/h) |
| RAM | 64-128 GB | Inclusa nell'istanza |
| Storage | 100 GB EBS | ~$10/mese |
| Tempo stimato | 3-6 giorni (32 core) | ~$100-200 totale |

### Software
| Requisito | Versione |
|-----------|----------|
| R | ≥ 4.3 |
| OS | Linux (Ubuntu 22.04+ / CentOS 8+) o macOS |
| git | Per clonare il repo |
| Internet | Solo per Script 00 (pacchetti) e Script 01 (dati) |

---

## 2. Setup Rapido

```bash
# 1. Clonare il repository
git clone <repo-url> feat_selection_stability
cd feat_selection_stability

# 2. Verificare R
R --version  # Deve essere ≥ 4.3

# 3. Installare pacchetti
Rscript R/00_install_packages.R 2>&1 | tee logs/00_install.log

# 4. Se horseshoe/ComplexHeatmap/treeshap falliscono:
Rscript -e '
install.packages("https://cran.r-project.org/src/contrib/Archive/horseshoe/horseshoe_0.2.0.tar.gz",
                 repos = NULL, type = "source")
BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
remotes::install_github("ModelOriented/treeshap", upgrade = "never")
'
```

---

## 3. Esecuzione Pipeline

### 3.1 Sequenza completa

```bash
# Creare directory per i log
mkdir -p logs

# Script 01-03: Download + Preprocessing + Parametri
# Tempo: ~15 minuti totale
Rscript R/01_download_data.R 2>&1 | tee logs/01_download.log
Rscript R/02_preprocess.R 2>&1 | tee logs/02_preprocess.log
Rscript R/03_extract_empirical_params.R 2>&1 | tee logs/03_params.log

# Script 04: Simulazione
# Tempo: ~30 min - 2 ore (dipende da n_cores)
Rscript R/04_simulate.R 2>&1 | tee logs/04_simulate.log

# Script 05: Feature Selection (BOTTLENECK)
# Tempo: 3-12 GIORNI con 8 core
# Usare nohup o screen/tmux!
nohup Rscript R/05_feature_selection.R > logs/05_fs.log 2>&1 &

# Script 07: Cross-validation (INDIPENDENTE da 04-06)
# Può girare in parallelo con 05
nohup Rscript R/07_cross_validation.R > logs/07_cv.log 2>&1 &

# Dopo completamento di 05:
Rscript R/06_metrics.R 2>&1 | tee logs/06_metrics.log

# Dopo completamento di 06 E 07:
Rscript R/08_figures.R 2>&1 | tee logs/08_figures.log
```

### 3.2 Pilot run (raccomandato prima della full run)

Modificare temporaneamente `config/config.yaml`:

```yaml
simulation:
  resampling:
    n_replications: 3    # Era 30
    n_bootstrap: 20      # Era 100
```

Eseguire Script 04 → 05 → 06, verificare risultati, poi ripristinare valori originali.

**Tempo pilot stimato:** ~30 min con 8 core.

### 3.3 Monitoraggio progresso

```bash
# Contare file checkpoint generati
ls data/simulated/sim_*.rds | wc -l              # Script 04
ls results/feature_selection/fs_*.rds | wc -l     # Script 05

# Confronto con atteso
# Script 04: 25 livelli × 30 rep = 750 file
# Script 05: 750 × 12 metodi = 9,000 file

# Tail del log
tail -f logs/05_fs.log
```

---

## 4. Configurazione n_cores

Modificare in `config/config.yaml`:

```yaml
project:
  n_cores: 8    # ← Cambiare in base al server
```

| Server | n_cores consigliato | Tempo Script 05 stimato |
|--------|--------------------|-----------------------|
| Laptop 8 core | 6-7 | ~12 giorni |
| Workstation 16 core | 14 | ~6 giorni |
| Cloud 32 core | 30 | ~3 giorni |
| HPC 64 core | 60 | ~1.5 giorni |

**Regola:** usare `n_cores - 2` per lasciare headroom al sistema.

---

## 5. Checkpointing e Ripresa

Tutti gli script supportano checkpointing automatico. Se un'esecuzione viene interrotta:

```bash
# Semplicemente rieseguire lo script — ripartirà da dove si era fermato
Rscript R/05_feature_selection.R 2>&1 | tee -a logs/05_fs.log
```

**Meccanismo:**
- Script 04: verifica `data/simulated/sim_<task_id>.rds` prima di simulare
- Script 05: verifica `results/feature_selection/fs_<task_id>_<method>.rds` prima di calcolare
- Script 07: verifica `results/cross_validation/cv_<acc>_<method>.rds`

**Per forzare ricalcolo:** cancellare i file .rds corrispondenti.

---

## 6. Stima Spazio Disco

| Directory | Stima | Note |
|-----------|-------|------|
| `data/raw/` | 1-5 GB | Dipende da quanti dataset scaricati |
| `data/processed/` | 100-500 MB | Dataset processati |
| `data/empirical_params/` | 50-200 MB | Matrici correlazione |
| `data/simulated/` | 5-15 GB | 750 dataset simulati |
| `results/feature_selection/` | 10-30 GB | 9,000 file con selection matrices |
| `results/metrics/` | 50-200 MB | Tabelle aggregate |
| `results/cross_validation/` | 500 MB - 2 GB | FS su dati reali |
| `results/figures/` | 10-50 MB | PDF vettoriali |
| **Totale** | **~20-50 GB** | |

---

## 7. Troubleshooting

### Pacchetto non si installa

```bash
# Controllare versione R
R --version

# Per pacchetti CRAN, provare mirror diverso
Rscript -e 'install.packages("pkg", repos = "https://cran.rstudio.com")'

# Per Bioconductor
Rscript -e 'BiocManager::install("pkg", force = TRUE)'

# Per GitHub con problemi SSL
Rscript -e 'options(download.file.method = "wget"); remotes::install_github("repo/pkg")'
```

### Script 01 fallisce (download)

1. Verificare connessione: `ping www.ebi.ac.uk`
2. Verificare firewall: porte 80, 443, 21 (FTP) aperte
3. Se FTP bloccato: i download usano anche REST API come fallback
4. Se persiste: scaricare manualmente da browser e mettere in `data/raw/<accession>/`

### Script 05 va in OOM (Out of Memory)

1. Ridurre `n_cores` in config
2. Lo script processa in chunk di `n_cores × 2` task alla volta con `gc()` tra i chunk
3. Se non basta: ridurre `n_bootstrap` a 50

### Script 05 è troppo lento

1. Aumentare `n_cores`
2. Ridurre scope: diminuire `n_replications` (30 → 20) o `n_bootstrap` (100 → 50)
3. Escludere metodi lenti: commentare horseshoe e/o spike_slab dal config
4. Usare cloud con più core per la run principale

---

## 8. Dopo l'Esecuzione

### File di output chiave

```
results/metrics/metrics_all.rds          # Master table: tutte le metriche per ogni scenario × metodo × rep
results/metrics/metrics_summary.rds      # Aggregato: media ± sd per scenario × metodo
results/cross_validation/cv_summary.rds  # Stabilità su dati reali
results/cross_validation/concordance_results.rds  # Concordanza cross-database
results/figures/fig_01..08.pdf           # Figure per il paper
results/session_info_*.txt               # Ambiente R per riproducibilità
```

### Trasferire risultati

```bash
# Comprimere risultati (escludendo dati grezzi e simulati)
tar czf results_export.tar.gz \
  results/metrics/ \
  results/cross_validation/ \
  results/figures/ \
  results/session_info_*.txt \
  config/config.yaml

# Dimensione: ~500 MB - 1 GB
```
