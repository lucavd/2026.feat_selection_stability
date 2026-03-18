# Inventario Dataset Pubblici

> Ultimo aggiornamento: 2026-03-18
> Status: FINALIZZATO — 9 dataset processati con successo

---

## 1. Fonti Dati

### 1.1 CIMCB Benchmark (fonte primaria — 7 dataset)

Pre-processed Excel files da **Mendez et al. (2019)**, "A comparative evaluation of the generalised predictive ability of eight machine learning algorithms across ten clinical metabolomics data sets for binary classification" (*Metabolomics* 15:150, DOI: 10.1007/s11306-019-1612-4).

- **Repository:** https://github.com/CIMCB/MetabComparisonBinaryML
- **Formato:** Excel con sheet "Data" (SampleID, Class, M1..Mn) + sheet "Peak" (nomi metaboliti)
- **Vantaggi:** Dati già puliti, classi binarie pre-validate, usati in benchmark pubblicati

### 1.2 MetaboLights (1 dataset)

- **API REST:** https://www.ebi.ac.uk/metabolights/ws/
- **Download:** FTP + API per MAF files (m_*.tsv) e sample sheets (s_*.txt)

### 1.3 Metabolomics Workbench (1 dataset)

- **API REST:** https://www.metabolomicsworkbench.org/rest/
- **Endpoint dati:** `/study/study_id/{ID}/data` (JSON nested)
- **Endpoint fattori:** `/study/study_id/{ID}/factors` (JSON nested, pipe-delimited)

---

## 2. Dataset Processati

| # | ID | Piattaforma | Fonte | n | p | Classe 0 | Classe 1 | Condizione |
|---|-----|------------|-------|--:|---:|----------|----------|------------|
| 1 | ST001047 | NMR | CIMCB | 83 | 149 | HE (40) | GC (43) | Gastric cancer vs healthy — urine |
| 2 | ST001706 | NMR | MW API | 256 | 50 | Control (174) | RCC (82) | Renal cell carcinoma vs control — urine |
| 3 | MTBLS28 | LC-MS | MetaboLights | 1005 | 1359 | Control (536) | Case (469) | NSCLC vs control — urine |
| 4 | MTBLS404 | LC-MS | CIMCB | 184 | 120 | Female (83) | Male (101) | Sacurine sex classification — urine |
| 5 | ST001000 | LC-MS | CIMCB | 107 | 1533 | CD (49) | UC (58) | IBD subtype (UC vs CD) — stool |
| 6 | MTBLS136 | LC-MS | CIMCB | 668 | 787 | E-only (331) | E+P (337) | Hormone use — serum |
| 7 | MTBLS92 | LC-MS | CIMCB | 253 | 138 | Pre-chemo (142) | Post-chemo (111) | Breast cancer treatment — plasma |
| 8 | ST000369 | GC-MS | CIMCB | 80 | 181 | Control (31) | Cancer (49) | Lung adenocarcinoma — serum |
| 9 | ST000496 | GC-MS | CIMCB | 100 | 69 | Pre (50) | Post (50) | Periodontal debridement — saliva |

### Copertura piattaforme
- **NMR:** 2 dataset (ST001047, ST001706)
- **LC-MS:** 5 dataset (MTBLS28, MTBLS404, ST001000, MTBLS136, MTBLS92)
- **GC-MS:** 2 dataset (ST000369, ST000496)

### Range dimensionali
- **Campioni (n):** 80 — 1005
- **Features (p):** 50 — 1533
- **Rapporto p/n:** 0.6 — 14.3

---

## 3. Dataset Esclusi

| ID | Motivo esclusione |
|----|------------------|
| MTBLS374 | MAF NMR con solo 3 features (dati bin non nel MAF); reprocessing con rDolphin necessario |
| MTBLS97 | API restituisce 403 Forbidden |
| MTBLS537 | API restituisce 403 Forbidden |
| MTBLS1 | MAF NMR formato non standard (no colonne campione numeriche) |
| MTBLS352 | MAF con solo 4 righe, 0 colonne numeriche |
| MTBLS733 | Benchmark spike-in con solo 2 colonne numeriche nel MAF |
| ST000388 | Solo 7 metaboliti via REST API |
| ST001047 (MW) | Factors solo "Sample_Type:Sample" — nessuna classe disease |
| ST001386 | File dati non trovato via REST API |

---

## 4. Riferimenti Benchmark

1. **Mendez, Reinke & Broadhurst (2019).** A comparative evaluation of the generalised predictive ability of eight machine learning algorithms across ten clinical metabolomics data sets for binary classification. *Metabolomics* 15:150. DOI: 10.1007/s11306-019-1612-4
2. **Labory, Njomgue-Fotso & Bottini (2024).** Feature selection in metabolomics. *Comput. Struct. Biotechnol. J.* 23:1274–1287
3. **Bifarin (2023).** Interpretable machine learning with SHAP for metabolomics. *PLoS ONE*
4. **Bifarin et al. (2021).** RCC ML prediction. *J. Proteome Research* 20(7):3629–3641

---

## 5. Note sul Preprocessing

- **CIMCB datasets:** già pre-processati; filtro per classi valide (escludere QC, NA), poi QC standard (missing rate, zero-variance)
- **MetaboLights (MTBLS28):** parsing MAF → matching sample sheet via "Sample Name" → binarizzazione classi
- **MW (ST001706):** JSON nested → matrice features × samples; factors pipe-delimited → estrazione classe "SampleType"
- **Normalizzazione:** log2 + autoscaling (default `log_auto` da config) applicata a tutti
- **Imputazione:** median imputation per valori mancanti residui
