# PROMPT — Deep Research: Bibliografia Annotata per Feature Selection Stability in Metabolomics

## Obiettivo
Produci un **file BibTeX completo** (`.bib`) con tutti i riferimenti necessari per il seguente articolo scientifico, più un **file di annotazioni** che per ogni entry spiega *cosa citare* e *dove nel paper*.

## Il nostro articolo
**"The stability illusion: A simulation-based assessment of feature selection methods for high-dimensional metabolomics data"**

Studio in silico: 9 dataset metabolomici pubblici → parametri empirici → 1590 dataset simulati (8 scenari) → 11 metodi di feature selection × 100 bootstrap × 30 repliche → metriche di stabilità (Nogueira), accuratezza (TPR/FDR), predizione (AUC), parsimonia. Target: Briefings in Bioinformatics.

### Risultati chiave
- Nessun metodo FS raggiunge stabilità alta (Nogueira 0.15–0.46)
- TPR <0.03 con effect size piccoli (fc=1.5), sale a 0.25–0.43 solo con fc≥2.0
- Knockoff: stabilità più alta ma converge solo con n > p (94% fallimenti)
- Wilcoxon FDR: miglior compromesso pratico
- Scenari semi-sintetici (spike-in su dati reali) danno risultati migliori degli MVN puri

### Metodi FS usati
wilcoxon_fdr, fold_change, volcano plot, LASSO, elastic net, Boruta, RF importance, stability selection, knockoff filter, spike-and-slab, SHAP (XGBoost)

### Dataset usati
7 CIMCB benchmark (Mendez et al. 2019): ST001047, MTBLS404, ST001000, MTBLS136, MTBLS92, ST000369, ST000496. 1 MetaboLights: MTBLS28. 1 MW: ST001706.

---

## Cosa cercare

Per OGNI riferimento trovato, fornisci:

1. **Entry BibTeX completa** (tipo @article/@inproceedings/@book, tutti i campi: author, title, journal, year, volume, pages, doi). Usa una chiave BibTeX leggibile (es. `nogueira2018stability`, `tibshirani1996lasso`).
2. **Citare per**: frase concisa che spiega PER COSA ESATTAMENTE citiamo questo paper (non un riassunto generico).
3. **Sezione del paper**: dove va citato (Introduction, Methods, Discussion, etc.)
4. **Affidabilità**: VERIFICATO (hai trovato il DOI funzionante) o DA VERIFICARE.

---

### CATEGORIA 1 — Metriche di stabilità (citare in Methods §Metrics + Introduction)

Cerca e verifica:

- **Nogueira, Sechidis & Brown (2018)** — JMLR. Citare per: definizione dell'indice di stabilità corretto per chance (Eq. 4), varianza jackknife, CI. È la nostra metrica primaria.
- **Kuncheva (2007)** — consistency index corretto per chance. Citare per: metrica secondaria, confronto con Nogueira.
- **Lustgarten et al. (2009)** — probabilmente il primo a proporre stabilità nel contesto FS biomedica. Citare per: contesto storico.
- **Jaccard (1901 o 1912)** — riferimento originale per indice di Jaccard. Citare per: metrica secondaria.
- **Dice (1945)** — coefficiente di Dice/Sørensen. Citare per: metrica secondaria.
- Cerca se **stabm** (pacchetto R CRAN) ha un paper associato. Citare per: implementazione metriche.

### CATEGORIA 2 — Metodi FS: paper originali (citare in Methods §Feature selection methods)

Per ciascun metodo usato nel nostro studio, trova il paper originale:

- **Wilcoxon rank-sum test** — Wilcoxon (1945) o Mann-Whitney (1947). Citare per: test statistico di base nel filtro univariato.
- **Benjamini & Hochberg (1995)** — FDR correction. Citare per: correzione per test multipli nel filtro Wilcoxon.
- **LASSO** — Tibshirani (1996), JRSS-B. Citare per: penalizzazione L1, selezione variabili embedded.
- **Elastic Net** — Zou & Hastie (2005), JRSS-B. Citare per: penalizzazione L1+L2, grouping effect su feature correlate.
- **glmnet** — Friedman, Hastie, Tibshirani (2010), JSS. Citare per: implementazione R usata.
- **Boruta** — Kursa & Rudnicki (2010), JMLR. Citare per: all-relevant feature selection via shadow features + RF.
- **Random Forest** — Breiman (2001), Machine Learning. Citare per: base learner di Boruta e RF importance.
- **ranger** — Wright & Ziegler (2017), JSS. Citare per: implementazione RF usata.
- **Stability selection** — Meinshausen & Bühlmann (2010), JRSS-B. Citare per: framework meta di subsampling + selezione, controllo PFER.
- **stabs** — Hofner, Boccuto, Göker (2015) o simile. Citare per: implementazione R di stability selection.
- **Knockoff filter** — Barber & Candès (2015), Annals of Statistics. Citare per: controllo FDR model-free.
- **Model-X Knockoffs** — Candès et al. (2018), JRSS-B. Citare per: estensione del knockoff filter, versione usata.
- **knockoff** (pacchetto R) — cerca se ha un paper associato o solo documentazione.
- **Spike-and-slab** — Mitchell & Beauchamp (1988) o George & McCulloch (1993). Citare per: prior bayesiano per variable selection.
- **spikeslab** (pacchetto R) — Ishwaran & Rao. Citare per: implementazione usata.
- **SHAP** — Lundberg & Lee (2017), NeurIPS. Citare per: Shapley values come importanza feature.
- **XGBoost** — Chen & Guestrin (2016), KDD. Citare per: base learner per SHAP.
- **Volcano plot** come criterio FS — cerca un riferimento formale. Potrebbe non esistere un singolo paper originale; in tal caso cerca un review che lo definisce nel contesto metabolomica.

### CATEGORIA 3 — Stabilità della FS nella letteratura (citare in Introduction + Discussion)

- **Saeys, Inza & Larrañaga (2007)** — review FS in bioinformatics, Bioinformatics. Citare per: tassonomia filter/wrapper/embedded, problema della instabilità menzionato.
- **He & Yu (2010)** — stable feature selection. Citare per: formalizzazione del tradeoff stabilità-accuratezza.
- **Boulesteix & Slawski (2009)** — stabilità FS in alta dimensionalità. Citare per: evidenza empirica che FS è instabile.
- **Haury et al. (2011)** — benchmark FS in genomica. Citare per: confronto tra metodi, evidenza instabilità.
- Cerca **paper 2020-2026** che studiano stabilità FS specificamente in metabolomica o proteomica — probabilmente pochi, questo è il nostro gap.
- Cerca **Guo et al.** o simili che hanno fatto benchmark FS in metabolomica recenti.

### CATEGORIA 4 — Metabolomica: dataset, riproducibilità, best practices (citare in Introduction + Methods)

- **Mendez, Reinke & Broadhurst (2019)** — Metabolomics 15:150. CIMCB benchmark datasets. Citare per: fonte di 7 dei 9 dataset usati.
- **Haug et al. (2020)** — MetaboLights database. Citare per: repository MTBLS28.
- **Sud et al. (2016)** — Metabolomics Workbench. Citare per: repository ST001706.
- **Sumner et al. (2007)** — MSI minimum reporting standards. Citare per: standard di riproducibilità.
- Cerca il paper che documenta **scarsa riproducibilità dei biomarker metabolomici** — "72% singleton findings" o statistica simile. Potrebbe essere **Trivedi et al.** o **Broadhurst & Kell (2006)**.
- **Broadhurst & Kell (2006)** — "Statistical strategies for avoiding false discoveries..." Metabolomics. Citare per: overfitting e falsi positivi in metabolomica.
- **Goodacre et al. (2007)** o **Dunn et al. (2011)** — best practices in metabolomics. Citare per: contesto preprocessing.

### CATEGORIA 5 — Simulazione (citare in Methods §Simulation)

- **mvtnorm** (pacchetto R) — Genz & Bretz. Citare per: generazione dati multivariati normali.
- **corpcor** — Schäfer & Strimmer (2005). Citare per: shrinkage della matrice di correlazione (Ledoit-Wolf).
- Cerca paper che descrivono **simulazione semi-sintetica** (spike-in su dati reali) in metabolomica o genomica. Chi lo ha proposto? Citare per: giustificazione approccio S8.
- Cerca paper sulle **proprietà distribuzionali dei dati metabolomici** (log-normalità). Citare per: giustificazione MVN → exp().
- Cerca **effect size tipici** in metabolomica caso-controllo — meta-analisi o review con range di fold change. Citare per: giustificazione scelta fc=1.2, 1.5, 2.0, 3.0.

### CATEGORIA 6 — Predizione e valutazione (citare in Methods §Evaluation metrics)

- **pROC** — Robin et al. (2011), BMC Bioinformatics. Citare per: calcolo AUC.
- **Matthews (1975)** — MCC. Citare per: metrica di accuratezza.
- **Chicco & Jurman (2020)** — vantaggi di MCC su F1. Citare per: giustificazione uso MCC.
- Cerca paper che discutono **AUC come metrica downstream** per valutare feature set. Citare per: giustificazione.

### CATEGORIA 7 — Software e riproducibilità (citare in Methods o Acknowledgements)

- **R Core Team** — citazione standard per R.
- **data.table** — Dowle & Srinivasan. Citare per: manipolazione dati.
- **ggplot2** — Wickham (2016). Citare per: visualizzazione.
- **here** — Müller (2020). Citare per: path management.
- Cerca se **Boruta**, **ranger**, **glmnet**, **stabs**, **knockoff**, **spikeslab** hanno paper JSS o JOSS associati.

### CATEGORIA 8 — Contesto generale, review influenti (citare in Introduction/Discussion)

- **Guyon & Elisseeff (2003)** — "An Introduction to Variable and Feature Selection", JMLR. Citare per: tassonomia classica.
- **Heinze, Wallisch & Dunkler (2018)** — variable selection. Citare per: problemi noti della variable selection in contesti predittivi.
- Cerca **review recenti (2023-2026) su feature selection in -omics**. Ci sono survey che coprono il landscape attuale?
- Cerca paper che discutono il **"researcher degrees of freedom"** nella feature selection metabolomica — il fatto che le scelte analitiche (preprocessing, metodo FS, soglie) influenzano drammaticamente i risultati.

---

## Formato output richiesto

### File 1: `references.bib`
File BibTeX standard con tutte le entry. Esempio formato:

```bibtex
@article{nogueira2018stability,
  author  = {Nogueira, Sarah and Sechidis, Konstantinos and Brown, Gavin},
  title   = {On the Stability of Feature Selection Algorithms},
  journal = {Journal of Machine Learning Research},
  year    = {2018},
  volume  = {18},
  number  = {174},
  pages   = {1--54},
  url     = {http://jmlr.org/papers/v18/17-514.html}
}
```

### File 2: `references_annotations.md`
Tabella markdown con annotazioni. Esempio formato:

```markdown
| Chiave BibTeX | Citare per | Sezione | Verificato |
|---------------|-----------|---------|------------|
| nogueira2018stability | Definizione indice di stabilità corretto per chance (Eq. 4); varianza jackknife e CI — nostra metrica primaria | Methods §2.4, Introduction §1.2 | ✅ |
| tibshirani1996lasso | Introduzione della penalizzazione L1 per selezione variabili; fondamento teorico del LASSO | Methods §2.2 | ✅ |
| mendez2019cimcb | Fonte di 7 dei 9 dataset metabolomici usati; benchmark validato per confronti ML in metabolomica | Methods §2.1 | ✅ |
```

---

## Regole critiche

1. **NON inventare DOI, URL, o numeri di pagina.** Se non riesci a verificare un campo, scrivi `note = {DOI DA VERIFICARE}` nella entry BibTeX e segnalo nella tabella annotazioni.
2. **Minimo 50 riferimenti**, massimo 80.
3. Ogni entry deve avere tutti i campi BibTeX standard compilati (author, title, journal/booktitle, year, volume, pages, doi).
4. Ordina le entry per categoria nel .bib (usa commenti `% === CATEGORIA N ===`).
5. Le chiavi BibTeX devono essere nel formato `cognome_primo_autoreANNOparola_chiave` (es. `nogueira2018stability`).
6. Nella tabella annotazioni, la colonna "Citare per" deve essere **specifica al nostro studio**, non un riassunto generico del paper.
