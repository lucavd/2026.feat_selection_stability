# PROMPT — Deep Research Phase: Feature Selection Stability in Metabolomics

## Chi sei
Sei un consulente senior di biostatistica e chemiometria con esperienza decennale in metabolomica e high-dimensional data analysis. Il tuo compito è fare una ricerca approfondita sul web e produrre come output un **prompt operativo estremamente dettagliato** destinato a Claude Code, che dovrà scrivere tutto il codice R (e eventualmente Python) per un articolo scientifico di simulazione.

## Contesto del progetto
Stiamo scrivendo un articolo scientifico intitolato provvisoriamente:

**"The stability illusion: A simulation-based assessment of feature selection methods for high-dimensional metabolomics data"**

L'articolo è interamente in silico. Non produciamo dati nuovi. Usiamo dati pubblici scaricati programmaticamente e simulazioni basate su strutture di correlazione e distribuzioni empiriche estratte da dati reali. Il target journal è Briefings in Bioinformatics, Metabolomics, o Analytical Chemistry.

L'autore principale è un biostatistico con esperienza diretta in metabolomica clinica (studi su trapianto polmonare, CLAD, dati NMR e LC-MS), quindi il paper deve riflettere consapevolezza delle problematiche reali del campo.

## Cosa devi ricercare (sii esaustivo)

### 1. Dati pubblici — accesso programmatico

#### MetaboLights (EMBL-EBI)
- Cerca la **REST API** di MetaboLights: endpoint esatti, documentazione corrente (URL), come interrogare gli studi, come scaricare i dati metabolomici processed (feature tables, non raw spectra).
- Identifica **5-8 dataset caso-controllo** specifici su MetaboLights che soddisfino TUTTI questi criteri: (a) design caso-controllo o two-group comparison, (b) feature table processed disponibile (matrice metaboliti × campioni), (c) almeno 30 campioni per gruppo, (d) metadata con covariate cliniche (età, sesso, BMI come minimo), (e) mix di piattaforme: almeno 2 NMR e 2-3 LC-MS.
- Per ciascun dataset: accession number (MTBLS...), piattaforma, patologia studiata, numero campioni, numero features, URL diretto per download, formato file.
- **Verifica che i dati siano effettivamente scaricabili** — MetaboLights ha avuto problemi di accesso in passato.

#### Metabolomics Workbench (NIH)
- Cerca la **REST API**: endpoint esatti (dovrebbe essere `https://www.metabolomicsworkbench.org/rest/` — **verifica**), documentazione, come scaricare study data programmaticamente.
- Identifica **3-5 dataset** con gli stessi criteri di sopra. Accession numbers (ST...), dettagli, URL.
- Cerca se esiste un pacchetto R per accedere a Metabolomics Workbench (`metabolomicsWorkbenchR` su Bioconductor? **Verifica se esiste e se è mantenuto**).

#### Alternative
- Cerca se **MassIVE** o **GNPS** hanno dataset metabolomici processed scaricabili via API.
- Cerca se su **Zenodo** ci sono dataset metabolomici benchmark con ground truth nota (mixture experiments, spike-in metabolomici).

### 2. Metodi di feature selection da confrontare
Per ciascuno, cerca: pacchetto R esatto, paper originale (DOI), come si usa, iperparametri chiave, assunzioni, limitazioni note.

**Filter methods:**
- **Wilcoxon rank-sum test** con correzione FDR (Benjamini-Hochberg)
- **Fold-change thresholding** (il metodo "barbarico" che molti usano)
- **Volcano plot-based selection** (combinazione FC + p-value)

**Wrapper/Embedded methods:**
- **LASSO** (`glmnet`): come funziona la selezione, come si sceglie lambda, problema della instabilità al resampling
- **Elastic net** (`glmnet`): differenze con LASSO, quando è meglio
- **Boruta** (`Boruta`): come funziona (shadow features), pacchetto R, parametri chiave, **limitazioni note** (cerca critiche pubblicate). Questo è particolarmente rilevante perché l'autore lo usa regolarmente.
- **Random forest importance** (`randomForest` o `ranger`): permutation importance vs impurity importance, bias noti (verso variabili con più categorie/valori unici)
- **Stability selection** (Meinshausen & Bühlmann, 2010): paper originale, implementazione in R (`stabs` package? **verifica**), come si parametrizza, vantaggi teorici sulla stabilità
- **Knockoff filter** (`knockoff` package in R): paper originale (Barber & Candès), come funziona, limitazioni pratiche per p >> n, se è applicabile a dati metabolomici

**Bayesian methods:**
- **Horseshoe prior** per regressione Bayesiana: implementazione in R (`brms`? `rstanarm`? pacchetto dedicato?), come si usa per feature selection, complessità computazionale
- **Spike-and-slab prior**: implementazione in R, praticabilità per p=1000+

**ML-based methods:**
- **SHAP-based feature selection**: calcolo SHAP values per un modello (es. XGBoost), selezione basata su mean |SHAP|. Pacchetti: `shapviz`, `treeshap`, o Python `shap`. Come si stabilizza?
- **Sparse autoencoders** per feature selection: architettura, implementazione in Python (PyTorch/Keras), come si estraggono le features "importanti" dal bottleneck layer
- **Variational autoencoders con sparsità**: differenze dal sopra, implementazione

**Ensemble/Meta methods:**
- Cerca se esistono metodi **ensemble di feature selection** che combinano i risultati di più metodi. Paper recenti?

### 3. Struttura di correlazione nei dati metabolomici
Questo è CRUCIALE per simulazioni realistiche.

- Cerca paper che documentino la **struttura di correlazione tipica** nei dati metabolomici: quanto sono correlati i metaboliti? Che tipo di struttura (block diagonal? hub structure? random?)?
- Cerca differenze nella struttura di correlazione tra **NMR vs LC-MS** (NMR ha correlazioni molto più forti per overlap spettrale)
- Come simulare matrici di covarianza realistiche? Cerca approcci:
  - Stima diretta da dati reali (sample covariance → simulazione multivariata)
  - Factor model (pochi fattori latenti + rumore)
  - Graphical LASSO per stimare la precision matrix
- Pacchetti R per simulare dati correlati ad alta dimensionalità: `mvtnorm`, `MASS::mvrnorm`, ma come gestire p >> n per la matrice di covarianza? Cerca soluzioni (shrinkage estimators, `corpcor`, `rags2ridges`).

### 4. Framework di simulazione — specifico per metabolomica
Cerca:

- **`MetSizeR`** o simili: pacchetti R per simulazione di dati metabolomici. Esistono?
- Come simulare **dati metabolomici realistici**: distribuzione tipica dei metaboliti (log-normale? gamma?), proporzione di missing values, pattern di missing (LOD-based = MNAR), trasformazioni tipiche (log, pareto scaling, auto-scaling).
- Come impiantare **segnali veri controllati**: modificare un sottoinsieme di metaboliti con effect size noto, mantenendo la struttura di correlazione. Questo è non-triviale — cerca paper che affrontino il problema.
- Cerca paper di simulazione specifici per metabolomica (2020-2025).

### 5. Scenari di simulazione da implementare
Per ciascuno, cerca evidenza empirica che sia un problema reale:

- **Rapporto p/n**: da p/n=5 (500 metaboliti, 100 campioni) a p/n=50 (1000 metaboliti, 20 campioni). Qual è il range realistico negli studi pubblicati?
- **Effect size**: quali fold-change sono realistici in metabolomica? Cerca paper che riportino distribuzioni di effect size in studi caso-controllo. (Nella mia esperienza: la maggior parte sono FC < 1.5, pochi arrivano a FC > 2)
- **Multicollinearità**: vari gradi, da correlazione moderata (r=0.3-0.5) a forte (r=0.7-0.9). Come interagisce con la feature selection?
- **Non-linearità**: interazioni tra metaboliti che influenzano l'outcome. Come simulare? I metodi lineari (LASSO) le perdono, quelli tree-based no?
- **Confounders**: età, BMI, farmaci che influenzano il metaboloma. Quanto è grande il loro effetto? Cerca paper.
- **Missing data**: prevalenza tipica nei dati metabolomici (cerca numeri), pattern (MNAR per LOD), metodi di imputazione comuni (KNN, minimum value, probabilistic minimum imputation). Come l'imputazione interagisce con la feature selection?
- **Preprocessing variability**: la scelta di normalizzazione (PQN, TIC, median) e scaling (auto, pareto, log) cambia i risultati della feature selection? Cerca evidenza.

### 6. Metriche di valutazione — specifiche per stabilità
Cerca definizione formale e implementazione in R per:

- **Jaccard index** su set di features selezionate da sottocampioni bootstrap
- **Kuncheva index** (versione corretta per chance del Jaccard per feature selection)
- **Average Tanimoto index**
- **Spearman correlation of feature rankings** tra resampling runs
- **Nogueira stability measure** (2018) — cerca questo paper, sembra essere lo stato dell'arte per misurare stabilità
- **Proportion of Stably Selected features** (appare in >50%, >80%, >95% dei resampling runs)
- **True Positive Rate e FDR** (quando conosci la ground truth in simulazione)
- **Prediction performance** su test set (AUC, balanced accuracy) delle features selezionate
- **Parsimonia**: numero di features selezionate (un metodo che seleziona 200 features su 1000 non è utile clinicamente)

### 7. Validazione cross-database
Cerca:
- È realistico aspettarsi che features selezionate in un dataset metabolomico si replichino in un dataset indipendente della stessa patologia? Cerca paper che abbiano provato.
- Quali patologie hanno multipli dataset metabolomici indipendenti su database pubblici? (diabete tipo 2? cancro? malattie cardiovascolari? cerca)
- Come valutare la replicabilità: stesse features? stessa direzione? stessa classificazione su dati nuovi?

### 8. Paper correlati — gap analysis
Cerca i seguenti e paper simili recenti (2022-2025):
- Paper di benchmark di feature selection per metabolomica (o omics in generale)
- **Boulesteix & Slaughter, 2009** — sulla stabilità della feature selection in high-dimensional data. DOI, cosa hanno fatto.
- **He & Yu, 2010** — stability of feature selection. Cerca.
- **Nogueira et al., 2018** — stability measures. Paper chiave.
- Paper recenti che critichino la riproducibilità dei biomarker metabolomici.
- Paper sul "vibration of effects" in omics (il concetto che piccole variazioni analitiche portano a risultati diversi).

### 9. Aspetti computazionali
- Per i metodi Bayesian: tempo stimato per feature selection su 1000 metaboliti con `brms` o `rstanarm`. È fattibile? Servono GPU?
- Per sparse autoencoders: PyTorch, GPU memory requirements per dataset metabolomici tipici
- Come parallelizzare efficientemente bootstrap resampling × metodi × scenari in R
- Stima complessiva: se facciamo 500 bootstrap resamples × 12 metodi × 8 dataset simulati × 7 scenari — quante ore/core?

## Output richiesto

Sulla base di tutta la ricerca fatta, produci un **singolo prompt dettagliatissimo** destinato a Claude Code che contenga:

1. **Struttura del progetto**: directory structure esatta
2. **Lista completa dei pacchetti** R/Python necessari con versioni
3. **Script 1 — Download dati**: codice per scaricare tutti i dataset metabolomici da MetaboLights e Metabolomics Workbench con error handling robusto, inclusi fallback se un dataset non è disponibile
4. **Script 2 — Preprocessing e armonizzazione**: standardizzazione dei dataset scaricati in un formato comune (matrice features × samples + metadata), gestione missing values, quality checks
5. **Script 3 — Estrazione parametri empirici**: stima delle distribuzioni marginali, struttura di correlazione, proporzione di zeri/missing, effect size tipici da ciascun dataset reale
6. **Script 4 — Simulazione**: framework di simulazione con TUTTI gli scenari parametrizzati, impianto di segnali noti, preservazione della struttura di correlazione, seed riproducibili
7. **Script 5 — Feature selection**: applicazione di TUTTI i metodi a tutti i dataset simulati, con bootstrap resampling per stabilità
8. **Script 6 — Metriche**: calcolo di tutte le metriche (TPR, FDR, stabilità, predizione, parsimonia)
9. **Script 7 — Validazione cross-database**: test di replicabilità su dataset reali indipendenti
10. **Script 8 — Figure e tabelle**: codice per generare tutte le figure per il paper (heatmap di stabilità, TPR/FDR curves per scenario, radar charts per confronto metodi, upset plots per overlap tra metodi)
11. **Un file di configurazione** centrale dove tutti i parametri dello studio sono definiti in un unico posto

Per ogni script, specifica:
- Funzione e ruolo esatto
- Input/output con path relativi
- Dipendenze da altri script
- Parametri configurabili
- Stima del tempo di esecuzione

Includi riferimenti concreti: URL di documentazione, DOI dei paper, accession number dei dataset, endpoint API esatti. **Non inventare nulla — ogni riferimento deve essere verificato nella tua ricerca.**

Il prompt deve essere sufficientemente dettagliato che Claude Code possa eseguirlo senza ambiguità e senza dover fare ricerche aggiuntive.
