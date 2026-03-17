# Inventario Dataset Pubblici per Feature Selection Stability Study

> Ultimo aggiornamento: 2026-03-17
> Status: DRAFT — da verificare accesso effettivo prima di finalizzare

---

## 1. Database e API Disponibili

### 1.1 MetaboLights (EMBL-EBI)

- **Sito web:** https://www.ebi.ac.uk/metabolights/
- **API REST:** https://www.ebi.ac.uk/metabolights/ws/
- **Documentazione API (Swagger):** https://www.ebi.ac.uk/metabolights/ws/api/spec.html
- **GitHub:** https://github.com/EBI-Metabolights/MtblsWS-Py
- **Studi pubblici:** ~1358+ (marzo 2026)
- **Autenticazione:** Non richiesta per studi pubblici; API token disponibile per studi privati

**Pacchetto R: `metabolighteR`**
- Repository: https://github.com/aberHRML/metabolighteR
- CRAN: https://cran.r-project.org/web/packages/metabolighteR/
- Versione: 0.1.4 (development su GitHub)
- Funzionalità: query studi pubblici, download metadati, accesso file ISA
- Limitazione: search API restituisce max 100 items per response

**Pacchetto Python: `metabolights-utils`**
- PyPI: metabolights-utils
- Funzionalità: download metadati studi pubblici, ricerca, gestione file ISA

**Note operative:**
- Nessun token API richiesto per la maggior parte delle chiamate su studi pubblici
- Documentazione interattiva disponibile via Swagger
- Aggiornamento gennaio 2025: nuovo framework di validazione e sistema di accession
- Contatto supporto: metabolights-help@ebi.ac.uk

---

### 1.2 Metabolomics Workbench (NIH)

- **Sito web:** https://www.metabolomicsworkbench.org/
- **API REST base URL:** `https://www.metabolomicsworkbench.org/rest/[context]/[input]/[value]/[output]/[format]`
- **Documentazione API v1.2 (luglio 2025):** https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.2.pdf
- **Documentazione API v1.1 (novembre 2023):** https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.1.pdf
- **Pagina REST service:** https://www.metabolomicsworkbench.org/tools/mw_rest.php
- **Studi pubblici:** 4144 (su 4610 totali, marzo 2026)
- **Sponsor:** NIH Common Fund

**Endpoint principali per accesso dati:**

| Endpoint | Descrizione |
|----------|-------------|
| `/study/study_id/ST000001/summary` | Overview dello studio |
| `/study/study_id/ST000001/factors` | Campioni e variabili sperimentali |
| `/study/study_id/ST000001/analysis` | Informazioni analisi |
| `/study/study_id/ST000001/metabolites` | Metaboliti identificati |
| `/study/study_id/ST000001/data` | Misurazioni metaboliti (feature table) |
| `/study/analysis_id/AN000001/datatable/` | Risultati metaboliti identificati |
| `/study/study_id/ST/available` | Lista tutti gli studi pubblici |

**Contesti validi:** "study", "compound", "refmet", "gene", "protein", "moverz", "exactmass"
**Formati output:** JSON (default, raccomandato), text ("txt")

**Pacchetto R: `metabolomicsWorkbenchR`**
- Bioconductor: https://www.bioconductor.org/packages/release/bioc/html/metabolomicsWorkbenchR.html
- GitHub: https://github.com/computational-metabolomics/metabolomicsWorkbenchR
- Installazione: `BiocManager::install("metabolomicsWorkbenchR")`
- Versione: 1.0.0+
- Funzione principale: `do_query(context, input_item, input_value, output_item)`
- Output: converte in SummarizedExperiment, MultiAssayExperiment, DatasetExperiment
- Compatibilità: struct e structToolbox packages

---

### 1.3 GNPS / MassIVE (UC San Diego)

- **GNPS:** https://gnps.ucsd.edu/
- **MassIVE:** https://massive.ucsd.edu/ProteoSAFe/
- **Documentazione API:** https://ccms-ucsd.github.io/GNPSDocumentation/api/
- **Dataset pubblici:** 1800+ (490,000+ file MS, 1.2+ miliardi spettri tandem MS)

**Accesso dati:**
- Download via FTP link per ogni dataset (subdirectory "ccms_peak" per dati processed)
- Library spectra: `https://external.gnps2.org/gnpslibraryjson`
- Formati: JSON, MGF, MSP, CSV
- Supporto USI (Universal Spectrum Identifier) per spettri individuali

**Limitazioni per il nostro studio:**
- Principalmente dati MS/MS, meno feature tables caso-controllo processed
- Nessun pacchetto R dedicato — accesso via HTTP diretto
- Utile come fonte supplementare, non primaria

---

### 1.4 Zenodo (Dataset Benchmark)

- **URL:** https://zenodo.org/
- **Risorse identificate:**
  - G-Aligner benchmark datasets: https://zenodo.org/records/7995790
  - MS-Net resources: https://zenodo.org/records/17669288
  - Viime visualization/integration: https://zenodo.org/records/4075591
  - Public MS/MS libraries: https://zenodo.org/records/11363475

---

## 2. Criteri di Selezione Dataset

Per il nostro studio, ogni dataset deve soddisfare **TUTTI** i seguenti criteri:

| Criterio | Requisito |
|----------|-----------|
| Design | Caso-controllo o two-group comparison |
| Feature table | Matrice metaboliti × campioni (processed, non raw spectra) |
| Campioni per gruppo | ≥ 30 |
| Metadata cliniche | Almeno età, sesso; preferibilmente anche BMI |
| Piattaforma | Mix: almeno 2 NMR e 2-3 LC-MS nel totale |
| Accessibilità | Scaricabile programmaticamente via API |
| Formato | Tabellare (CSV, TSV, mwTab) |

---

## 3. Dataset con Ground Truth Nota (Benchmark)

Questi dataset sono cruciali per validare le metriche di TPR/FDR perché hanno un insieme noto di features "vere".

### 3.1 MTBLS733 (MetaboLights)

- **Accession:** MTBLS733
- **Piattaforma:** Thermo Q Exactive HF (LC-MS)
- **Design:** Due miscele standard (SA, SB) da Piper nigrum, 5 replicati ciascuna
- **Features totali:** 1100 composti con rapporti di concentrazione definiti
- **Ground truth:** 130 composti con variazione significativa di concentrazione; 836 features uniche identificate
- **Uso:** Benchmark per feature detection, quantificazione, marker selection
- **Riferimento:** https://www.omicsdi.org/dataset/metabolights_dataset/MTBLS733
- **Paper associato:** BMC Journal of Cheminformatics, 2023 (DOI da verificare)

### 3.2 Apple Spike-in Dataset

- **Design:** 10 controllo + 10 spiked samples
- **Features totali:** 1632
- **Ground truth:** 22 features da molecole spiked note
- **Spiked compounds:** Includono acido aspartico, acido malico, altri
- **Uso:** Benchmark per software evaluation (MS-Dial, MZmine 2, XCMS, MarkerView, Compound Discoverer)
- **Accession/URL:** Da confermare — cercare su MetaboLights o Metabolomics Workbench

---

## 4. Dataset Caso-Controllo da Identificare

### TODO: Identificazione Specifica

I dataset specifici devono essere identificati interrogando le API. La strategia è:

1. **Query Metabolomics Workbench:**
   ```
   GET /study/study_id/ST/available
   ```
   Filtrare per: study_type = "case-control" o "disease", n_samples >= 60

2. **Query MetaboLights:**
   Usare `metabolighteR` per listare studi e filtrare per criteri

3. **Patologie target** (con maggiore probabilità di avere multipli dataset indipendenti):
   - Diabete tipo 2 (molto probabile multipli dataset NMR e LC-MS)
   - Cancro (vari tipi — colorettale, mammella, polmone)
   - Malattie cardiovascolari
   - Obesity/sindrome metabolica
   - NAFLD/NASH

4. **Obiettivo:** 5-8 da MetaboLights + 3-5 da Metabolomics Workbench

### Nota sull'accesso effettivo

MetaboLights ha avuto problemi di accesso in passato. Per ogni dataset identificato, bisogna verificare:
- [ ] Feature table effettivamente scaricabile
- [ ] Formato compatibile (non solo raw spectra)
- [ ] Metadata complete (almeno gruppo caso/controllo + covariate demografiche)

---

## 5. Cross-Database Validation — Prerequisiti

Per lo Script 07 (validazione cross-database), servono **coppie di dataset indipendenti della stessa patologia**.

**Evidenze dalla letteratura:**
- Studio cross-biobank su 700,217 partecipanti (Nature Communications 2024, DOI: 10.1038/s41467-024-54357-0): NMR metabolomico replicato across tre biobank nazionali — dimostra fattibilità
- La replicabilità dei biomarker metabolomici è bassa: 72% dei metaboliti significativi riportati in un solo studio (PMC11999569)

**Patologie con multipli dataset probabili:**
1. **Diabete tipo 2** — la più studiata in metabolomica, alta probabilità di trovare ≥2 dataset indipendenti
2. **Cancro colorettale** — diversi studi metabolomici pubblicati
3. **Malattie cardiovascolari** — studi di coorte con metabolomica

---

## 6. Riferimenti Chiave

- Nature Communications 2025: Pan-repository universal identifiers per metabolomica (DOI: 10.1038/s41467-025-60067-y)
- MetaboLights: Haug et al. (2020) Nucleic Acids Research
- Metabolomics Workbench: Sud et al. (2016) Nucleic Acids Research
