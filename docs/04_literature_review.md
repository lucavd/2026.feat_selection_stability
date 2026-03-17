# Literature Review — Feature Selection Stability in Metabolomics

> Ultimo aggiornamento: 2026-03-17
> Status: DRAFT

---

## 1. Paper Fondamentali sulla Stabilità della Feature Selection

### 1.1 Meinshausen & Bühlmann (2010) — Stability Selection

- **Titolo:** Stability Selection
- **Journal:** Journal of the Royal Statistical Society Series B, 72(4):417-473
- **DOI:** 10.1111/j.1467-9868.2010.00740.x
- **URL:** https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9868.2010.00740.x
- **Contributo:**
  - Introduce il framework "stability selection": subsampling ripetuto + selezione ad alta dimensionalità
  - Controllo PFER (Per-Family Error Rate) in campione finito
  - Dimostra che stability selection raggiunge consistenza anche quando le condizioni di consistenza del LASSO sono violate
  - Applicazione a variable selection e modelli grafici Gaussiani
- **Rilevanza per il nostro studio:** Uno dei 12 metodi da confrontare; il paper più citato sull'argomento; il nostro studio testa empiricamente le promesse teoriche

### 1.2 Nogueira, Sechidis & Brown (2018) — Stability Measures

- **Titolo:** On the Stability of Feature Selection Algorithms
- **Journal:** Journal of Machine Learning Research, 18(174):1-54
- **URL:** https://jmlr.org/papers/v18/17-514.html
- **Contributo:**
  - Identifica almeno 15 diverse misure di stabilità usate nella letteratura (2002-2018)
  - Framework teorico basato su proprietà statistiche chiave
  - **Primo lavoro a permettere CI e test di ipotesi sulla stabilità** → confronto rigoroso tra algoritmi
  - Codice open-source per calcolo stabilità, varianza, CI, hypothesis tests
- **Rilevanza:** La metrica di stabilità primaria del nostro studio; il framework teorico su cui basiamo le analisi

### 1.3 He & Yu (2010) — Stable Feature Selection for Biomarker Discovery

- **Titolo:** Stable Feature Selection for Biomarker Discovery
- **Journal:** Computational Biology and Chemistry, 34(4):215-225
- **Autori:** Zengyou He, Weichuan Yu
- **URL:** https://pubmed.ncbi.nlm.nih.gov/20702140/
- **Contributo:**
  - Review completa con framework gerarchico generico
  - Identifica tre fonti di instabilità: (1) algoritmi non progettati per stabilità, (2) set multipli di veri marker, (3) maledizione della dimensionalità
  - Categorizzazione dei metodi esistenti
- **Rilevanza:** Framework concettuale; giustifica perché servono studi come il nostro

### 1.4 Boulesteix & Slawski (2009) — Stability and Aggregation of Ranked Gene Lists

- **Titolo:** Stability and Aggregation of Ranked Gene Lists
- **Journal:** Briefings in Bioinformatics
- **Contributo:**
  - Uno dei primi lavori sulla stabilità della selezione in dati high-dimensional
  - Identifica il problema del publication bias: ricercatori scelgono indici di stabilità che fanno apparire il loro algoritmo più stabile
- **Rilevanza:** Antecedente storico; motiva la necessità di misure standardizzate (poi fornite da Nogueira 2018)

### 1.5 Barber & Candès — Knockoff Filter

- **Paper 1:** "Controlling the false discovery rate via knockoffs"
  - **ArXiv:** 1404.5609
  - **Journal:** Annals of Statistics, 43(5):2055-2085, 2015
  - **Contributo:** Introduce il knockoff filter per controllo FDR nella selezione di variabili, condizione n ≥ p

- **Paper 2:** "A knockoff filter for high-dimensional selective inference"
  - **ArXiv:** 1602.03574
  - **Contributo:** Estensione al caso p >> n tramite sample splitting

- **Paper 3:** Candès, Fan, Janson & Lv (2018), "Panning for gold: model-X knockoffs for high dimensional controlled variable selection"
  - **Journal:** JRSS-B 80(3):551-577
  - **Contributo:** Model-X knockoffs — non richiedono modello lineare; assumono distribuzione nota delle features

---

## 2. Crisis della Riproducibilità in Metabolomica

### 2.1 Meta-Analisi sulla Riproducibilità (PMC11999569)

- **Titolo (probabile):** A Reproducibility Crisis for Clinical Metabolomics Studies
- **URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11999569/
- **Dati chiave:**
  - Meta-analisi di **244 studi metabolomici clinici**
  - **2206 metaboliti unici** identificati come statisticamente significativi
  - **72% riportati in un solo studio** (non replicati)
  - Di questi singleton, **85% sono rumore statistico**
  - Soglie p-value aggressive nonostante modesto miglioramento con correzione FDR
- **Rilevanza:** MOTIVAZIONE PRINCIPALE del nostro studio. Dimostra quantitativamente il problema che affrontiamo.

### 2.2 Sfide nella Validazione di Biomarker Metabolomici

- **URL:** https://www.mdpi.com/2218-1989/14/4/200
- **Problemi identificati:**
  - Protocolli non uniformi per raccolta e processamento campioni
  - Alta variabilità biologica (età, stile di vita, dieta, microbioma)
  - Mancanza di standardizzazione nell'analisi dati
- **Soluzione proposta:** Cross-validazione su coorti indipendenti multiple e protocolli standardizzati

### 2.3 Riproducibilità dei Biomarker

- **URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC10979063/
- **Finding:** La quantificazione relativa non è facilmente riproducibile; biomarker fluttuano tra batch e tra laboratori
- **Requisito:** Replicazione in altre popolazioni e test di riproducibilità nel tempo

---

## 3. Benchmark Recenti di Feature Selection (2020-2025)

### 3.1 Benchmarking Feature Selection for Metabolomics (2024)

- **PMID:** 38560281
- **URL:** https://pubmed.ncbi.nlm.nih.gov/38560281/
- **Risultato:** Feature extraction lineare e non-lineare accoppiata con feature selection supervisionata migliora significativamente la classificazione in dataset metabolomici
- **Rilevanza:** Benchmark diretto, ma NON valuta stabilità → gap che il nostro studio colma

### 3.2 Multi-Omics Feature Selection Benchmark (2022)

- **Journal:** BMC Bioinformatics
- **URL:** https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04962-x
- **Design:** 4 metodi filter + 2 embedded + 2 wrapper su 15 dataset TCGA
- **Rilevanza:** Grande scala ma su genomica, non metabolomica; non valuta stabilità

### 3.3 Microbiome-Metabolome Integration (2025)

- **Journal:** Communications Biology
- **URL:** https://www.nature.com/articles/s42003-025-08515-9
- **Finding:** Random Forest + NNLS migliore performance; tree-based methods con feature selection consistente
- **Rilevanza:** Approccio ensemble su dati multi-omics

### 3.4 ML per Microbiome-Metabolomics con Stable Feature Selection (2025)

- **URL:** https://www.biorxiv.org/content/10.1101/2025.06.21.660858v2
- **Rilevanza:** Uno dei pochi studi recenti che affronta esplicitamente la stabilità in contesto metabolomico

---

## 4. Effect Size e Fold-Change in Metabolomica

### 4.1 Large-Scale Analysis of Targeted Metabolomics Data

- **Journal:** Metabolomics
- **URL:** https://link.springer.com/article/10.1007/s11306-019-1564-8
- **Finding chiave:**
  - Fold-change medio per la maggior parte dei metaboliti: **< 2** in studi caso-controllo
  - Eccezioni: omocisteina (>20-fold), lattato/gliceraldeide (>9-fold)
  - Soglie tipiche: FC 1.5-2.0 per caso-controllo; 1.2-1.5 per studi meno contrastati
  - Alta variabilità da fonti multiple di incertezza complica la rilevazione di cambiamenti piccoli
- **Rilevanza:** Giustifica i nostri scenari di simulazione con FC < 1.5 come caso dominante

### 4.2 Uncertainty Budgeting in Fold Change

- **Journal:** The Analyst (RSC)
- **URL:** https://pubs.rsc.org/en/content/articlehtml/2017/an/c6an01342b
- **Punto chiave:** Il test di significatività non contiene informazione sulla magnitudine del cambiamento; la rilevanza metabolica richiede una "worthwhile amount of change"

---

## 5. Preprocessing e Impatto sulla Feature Selection

### 5.1 Comprehensive Evaluation (2021-2023)

- **Journal:** Metabolites
- **URL:** https://www.mdpi.com/2218-1989/12/3/202
- **Risultati:**
  - La scelta del preprocessing è determinata da ipotesi biologica, caratteristiche dataset, e metodo statistico
  - **Critical finding:** preprocessing e strategie di analisi dati sono determinanti critici della qualità
  - Random Forest: 27% misclassificazione più alta su dati originali vs feature-selected
  - Step inclusi: handling zero/missing, outlier detection, normalizzazione, centering, scaling, trasformazione
  - **Servono sensitivity analyses per valutare robustezza**
- **Rilevanza:** Giustifica lo scenario "preprocessing variability" nella nostra simulazione

### 5.2 Statistical Workflow for Feature Selection

- **Journal:** Metabolites
- **URL:** https://www.mdpi.com/2218-1989/9/7/143
- **Enfasi:** Selezione e applicazione corretta degli algoritmi cruciale per buoni modelli

---

## 6. Missing Data in Metabolomica

### 6.1 Characterization of Missing Values in Untargeted MS

- **Journal:** Metabolomics
- **URL:** https://link.springer.com/article/10.1007/s11306-018-1420-2
- **Classificazione:**
  - **MCAR:** errori casuali
  - **MAR:** determinati da altre variabili osservate
  - **MNAR:** valori censurati — **il più comune in metabolomica** (concentrazione sotto LOD)
- **Handling:** Sostituzione con valore LOD, zero, o minimo rilevato

### 6.2 Mechanism-Aware Imputation

- **URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9109373/
- **Approccio two-step:** Identificare tipo di missingness → applicare imputazione appropriata
- **Tool:** BayesMetab per modellazione Bayesiana dei missing values

### 6.3 imputomics (2024)

- **Journal:** Bioinformatics
- **URL:** https://academic.oup.com/bioinformatics/article/40/3/btae098/7611648
- **Tool:** Web server e pacchetto R per imputazione missing values in metabolomica

---

## 7. Struttura di Correlazione NMR vs LC-MS

### 7.1 Coverage e Complementarità

- **Journal:** Scientific Reports
- **URL:** https://www.nature.com/articles/s41598-023-43056-3
- **Copertura metaboliti:**
  - 46% esclusivi GC-MS
  - 34% esclusivi LC-MS
  - 15% esclusivi NMR
  - Solo 6% rilevati da più di una tecnica
- **NMR:** Non distruttivo, altamente riproducibile, strutturale — ma **correlazioni più forti** (overlap spettrale)
- **LC-MS:** Alta sensibilità, distruttivo, massima copertura metaboliti

### 7.2 Implicazioni per Simulazione

- NMR: simulare con correlazioni a blocchi forti (r = 0.7-0.9) — metaboliti dello stesso pathway o con overlap spettrale
- LC-MS: correlazioni più moderate (r = 0.3-0.6) — pathway metabolici condivisi ma meno overlap tecnico
- Struttura tipica: block-diagonal con hub (metaboliti centrali altamente connessi)

---

## 8. Ensemble Feature Selection

### 8.1 MVFS-SHAP (2025)

- **Titolo:** Majority Voting with SHAP for Stable Feature Selection
- **Journal:** Computers in Biology and Medicine
- **URL:** https://www.sciencedirect.com/science/article/abs/pii/S0169260725005863
- **Innovazione:** Combina frequenza di occorrenza con SHAP values per valutazione multi-dimensionale
- **Target:** Dati metabolomici high-dimensional, small-sample, con collinearità
- **Rilevanza:** Metodo recente che potremmo confrontare o usare come ispirazione per il nostro ensemble

### 8.2 Synergy Mechanistic + Ensemble (2024)

- **URL:** https://link.springer.com/chapter/10.1007/978-3-031-90714-2_9
- **Applicazione:** Eterogeneità metabolica cancro colorettale
- **Approccio:** Ensemble multi-fase con modelli meccanicistici

---

## 9. Validazione Cross-Database

### 9.1 Cross-Biobank Metabolomic Prediction (2024)

- **Journal:** Nature Communications
- **DOI:** 10.1038/s41467-024-54357-0
- **URL:** https://www.nature.com/articles/s41467-024-54357-0
- **Scala:** 700,217 partecipanti, tre biobank nazionali
- **Metodo:** NMR metabolomico su sangue
- **Risultato:** Score metabolomici identificano gruppi ad alto rischio con **replicazione cross-biobank consistente**
- **Rilevanza:** Dimostra fattibilità della validazione cross-database; è il nostro benchmark per lo Script 07

### 9.2 Replication in Epidemiological Studies

- **URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC11642057/
- **Importanza:** Step essenziale nell'identificazione di biomarker metabolomici
- **Integrazione:** ML migliora discovery identificando pattern per classificazione e predizione del rischio

---

## 10. Paper di Simulazione per Metabolomica

### 10.1 SimOmics (2025, emergente)

- **ArXiv:** 2507.09967
- **Capabilities:** Generazione dataset sintetici multi-omics multivariati, fattori latenti, strutture di sparsità, covarianza a blocchi, noise biologicamente ispirato
- **Status:** Emergente, da verificare disponibilità su CRAN/Bioconductor

### 10.2 mzrtsim (2023)

- **URL:** https://www.biorxiv.org/content/10.1101/2023.11.14.567024v3
- **Capabilities:** Simulazione raw data LC-MS/GC-MS con picchi cromatografici
- **Uso:** Benchmarking software (XCMS, mzMine, OpenMS) vs ground truth
- **Nota:** Simula raw data, non feature tables — utile come riferimento ma non direttamente per il nostro studio

### 10.3 Simulated LC-MS Datasets (2023, Analytical Chemistry)

- **URL:** https://pubs.acs.org/doi/10.1021/acs.analchem.3c04979
- **Features:** Localizzazione picchi nota, differenze metaboliti definite tra gruppi, rumore Gaussiano variabile, features mancanti
- **Rilevanza:** Approccio più vicino al nostro — simulazione a livello di feature table

---

## 11. Gap Analysis — Cosa Manca nella Letteratura

| Gap | Evidenza | Il nostro contributo |
|-----|----------|---------------------|
| Nessuno studio confronta stabilità di >10 metodi su dati metabolomici | Benchmark esistenti confrontano 3-4 metodi al massimo | Confronto sistematico di 12 metodi |
| Stabilità misurata con indici non standardizzati | Nogueira (2018) ha proposto il framework ma pochi lo applicano | Uso sistematico di Nogueira con CI |
| Simulazioni non preservano struttura di correlazione reale | La maggior parte usa dati iid o correlazioni semplici | Struttura di correlazione estratta da dati reali |
| Interazione preprocessing × feature selection non studiata | Solo studi osservazionali, non sistematici | Scenario dedicato nella simulazione |
| Missing data impact sulla stabilità non quantificato | Noto il problema, non quantificato sistematicamente | Scenario MNAR dedicato |
| Metodi Bayesiani mai confrontati con metodi classici in metabolomica | Horseshoe/spike-and-slab usati in genomica, quasi mai in metabolomica | Inclusione di 2 metodi Bayesiani |
| Cross-database validation raramente tentata | La meta-analisi mostra 72% non-replicazione | Script dedicato (Script 07) |

---

## 12. Riferimenti Pan-Repository

- Pan-repository universal identifiers (Nature Communications 2025): DOI 10.1038/s41467-025-60067-y
- Comprehensive computational metabolomics review (2025): https://www.sciencedirect.com/science/article/pii/S2001037025002806
