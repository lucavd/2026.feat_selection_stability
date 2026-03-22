# PROMPT — Deep Research Phase: DE Analysis Robustness Under Assumption Violations

## Chi sei
Sei un consulente senior di biostatistica e bioinformatica con esperienza decennale in RNA-seq analysis. Il tuo compito è fare una ricerca approfondita sul web e produrre come output un **prompt operativo estremamente dettagliato** destinato a Claude Code, che dovrà scrivere tutto il codice R (e eventualmente Python) per un articolo scientifico di simulazione.

## Contesto del progetto
Stiamo scrivendo un articolo scientifico intitolato provvisoriamente:

**"How robust is your differential expression analysis? A systematic simulation study under realistic assumption violations in bulk RNA-seq"**

L'articolo è interamente in silico. Non produciamo dati nuovi. Usiamo dati pubblici scaricati programmaticamente e simulazioni basate su distribuzioni empiriche estratte da dati reali. Il target journal è Genome Biology o Nature Methods.

## Cosa devi ricercare (sii esaustivo)

### 1. Dati pubblici — accesso programmatico
Cerca e verifica che i seguenti siano ancora accessibili e documenta esattamente come accedervi:

- **GTEx**: come scaricare count matrices gene-level via `recount3` (pacchetto Bioconductor). Verifica la versione corrente del database, quali tessuti sono disponibili, il formato dei dati, e il codice R minimo per scaricare i dati di 3-4 tessuti (es. whole blood, lung, liver, brain cortex). Cerca la documentazione ufficiale di `recount3` e eventuali tutorial recenti (2023-2025).

- **SEQC/MAQC-III spike-in dataset**: questo è il gold standard per validazione di metodi DE perché ha ground truth nota. Cerca: (a) dove si scarica oggi (GEO accession number, dovrebbe essere GSE49712 o simile — **verifica**), (b) come accedervi via `GEOquery` o altro metodo programmatico in R, (c) cosa contiene esattamente (quanti campioni, quali spike-in, quale piattaforma), (d) se ci sono paper recenti (2022-2025) che lo hanno usato per benchmark simili al nostro.

- **ERCC spike-in controls**: cerca se esistono altri dataset con ERCC spike-ins su GEO che possiamo usare come validazione indipendente. Documenta accession numbers e come scaricarli.

- Cerca se esistono **altri dataset con ground truth nota** per DE analysis (es. mixture experiments, titration series) pubblicati dopo il 2020.

### 2. Metodi DE da confrontare
Cerca la documentazione aggiornata, il pacchetto R/Bioconductor, e le assunzioni esplicite di ciascuno:

- **DESeq2**: versione corrente, assunzioni (NB distribution, shrinkage priors specifici), parametri chiave, come gestisce outliers, come gestisce low counts. Cerca se ci sono stati aggiornamenti importanti dal 2023 in poi.
- **edgeR**: versione corrente, differenze tra exact test, GLM, e quasi-likelihood approach. Assunzioni specifiche.
- **limma-voom**: versione corrente, come funziona la trasformazione voom, assunzioni sul modello lineare sottostante.
- **glmGamPoi**: cerca questo pacchetto — è relativamente nuovo, progettato per large-scale data. Come si differenzia da DESeq2? Quando è meglio?
- **dream** (dal pacchetto `variancePartition`): specifico per design con misure ripetute/campioni correlati. Documentazione, come si usa, cosa assume.
- **MAST**: originariamente per single-cell ma usato anche per bulk con zero-inflation. Documentazione corrente.
- **NOISeq**: metodo non-parametrico. Ancora mantenuto? Come performa?
- **Wilcoxon rank-sum test applicato a RNA-seq**: cerca paper recenti che lo propongono come alternativa robusta (c'è stato un paper controverso su questo — trovalo).

Per ciascun metodo cerca:
- Pacchetto R esatto e come installarlo
- Paper originale (DOI)
- Eventuali paper di benchmark recenti (2022-2025) che lo includono
- Limitazioni note documentate

### 3. Framework di simulazione
Cerca in dettaglio:

- **`splatter` (Bioconductor)**: documentazione corrente, come funziona `splatSimulate()`, quali parametri controlla, come si stimano i parametri da dati reali (`splatEstimate()`), limitazioni note.
- **`powsimR`**: ancora mantenuto? Documentazione, come si usa per simulare scenari DE, cosa può e non può simulare.
- **`compcodeR`**: altro framework di simulazione per RNA-seq. Stato corrente, funzionalità.
- **`muscat`**: per simulazione di dati multi-sample multi-group da single-cell, ma può essere usato per bulk? Cerca.
- **Simulazione semi-parametrica custom**: cerca paper recenti (2023-2025) che propongono approcci di simulazione semi-parametrica per RNA-seq dove si campiona dalle distribuzioni empiriche reali perturbandole in modo controllato. Questo è probabilmente l'approccio più realistico.

### 4. Scenari di violazione delle assunzioni
Per ciascuno dei seguenti scenari, cerca:
- Evidenza empirica che il problema esista nei dati reali (paper con dati)
- Come simularlo realisticamente
- Se qualcuno ha già studiato l'impatto sulla DE analysis

Scenari:
- **Zero-inflation oltre il modello NB**: quanto è prevalente nel bulk RNA-seq? In quali condizioni?
- **Correlazione tra campioni**: campioni dallo stesso paziente (misure ripetute), stessa famiglia, stesso batch. Come si struttura questa correlazione?
- **Batch effects non corretti**: cerca quanto sono grandi tipicamente (quale proporzione della varianza spiegano) nei dataset reali.
- **Confounders nascosti**: età, sesso, composizione cellulare (cell type deconvolution). Cerca `RUVSeq` e surrogate variable analysis (`sva`) come metodi per gestirli e come simulare la loro presenza.
- **Outlier samples**: frequenza tipica nei dataset reali, come simularli realisticamente.
- **Asimmetria dei fold-change**: in realtà la maggior parte dei geni ha fold-change vicino a zero e pochi hanno FC grandi. Cerca distribuzioni empiriche dei log2FC da studi reali.
- **Sample size piccoli**: n=3 per gruppo (comune in molti studi) vs n=5, n=10, n=20. Cerca raccomandazioni esistenti e come il sample size interagisce con le altre violazioni.

### 5. Metriche di valutazione
Cerca come vengono calcolate e quali pacchetti usare per:
- **FDR reale vs FDR nominale** (calibrazione dei p-value)
- **True Positive Rate (sensitivity/recall)**
- **Partial AUC** (per ROI a basso FDR)
- **Stabilità della gene list** su bootstrap resampling (Jaccard index, overlap coefficient)
- **Brier score** per la calibrazione
- Cerca se esiste un pacchetto R che già calcola queste metriche nel contesto DE benchmarking (es. `iCOBRA` — verifica se è ancora attivo e come si usa).

### 6. Paper correlati — gap analysis
Cerca i seguenti paper e paper simili recenti:
- **Soneson & Robinson, 2018** — benchmark di metodi DE (il paper di riferimento). DOI esatto, cosa hanno fatto, cosa NON hanno fatto che noi facciamo.
- Qualsiasi benchmark DE pubblicato 2022-2025. Cosa manca? Dove è il gap?
- Paper su robustezza di DESeq2 sotto violazioni specifiche.
- Paper che critica l'uso di p-value/FDR in DE analysis.
- Il paper controverso sul Wilcoxon test per RNA-seq (circa 2022-2024).

### 7. Aspetti computazionali
- Come parallelizzare le simulazioni in R (`BiocParallel`, `future`, `foreach`)?
- Stima realistica del tempo computazionale: se simulo 1000 dataset per scenario, con 10 scenari, e 8 metodi DE per dataset — quante ore/core servono?
- Come strutturare l'output delle simulazioni per analisi successiva (formato dati suggerito).

## Output richiesto

Sulla base di tutta la ricerca fatta, produci un **singolo prompt dettagliatissimo** destinato a Claude Code che contenga:

1. **Struttura del progetto**: directory structure esatta
2. **Lista completa dei pacchetti** R/Python necessari con versioni
3. **Script 1 — Download dati**: codice per scaricare tutti i dati pubblici necessari con error handling
4. **Script 2 — Stima parametri empirici**: estrazione delle distribuzioni empiriche dai dati GTEx
5. **Script 3 — Simulazione**: framework di simulazione con TUTTI gli scenari parametrizzati, con seed riproducibili
6. **Script 4 — Analisi DE**: applicazione di tutti i metodi a tutti i dataset simulati
7. **Script 5 — Metriche e valutazione**: calcolo di tutte le metriche
8. **Script 6 — Validazione su spike-in**: verifica su dati con ground truth
9. **Script 7 — Figure e tabelle**: codice per generare tutte le figure per il paper
10. **Un file di configurazione** centrale dove tutti i parametri dello studio sono definiti in un unico posto

Per ogni script, specifica:
- Funzione e ruolo esatto
- Input/output con path relativi
- Dipendenze da altri script
- Parametri configurabili
- Stima del tempo di esecuzione

Includi riferimenti concreti: URL di documentazione, DOI dei paper, accession number dei dataset, endpoint API esatti. **Non inventare nulla — ogni riferimento deve essere verificato nella tua ricerca.**

Il prompt deve essere sufficientemente dettagliato che Claude Code possa eseguirlo senza ambiguità e senza dover fare ricerche aggiuntive.
