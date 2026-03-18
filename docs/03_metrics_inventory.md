# Inventario Metriche di Valutazione

> Ultimo aggiornamento: 2026-03-18
> Status: DRAFT

---

## Overview

Le metriche sono organizzate in 4 categorie:
1. **Stabilità** — quanto i risultati cambiano al resampling
2. **Accuratezza della selezione** — TPR, FDR (quando ground truth è nota)
3. **Performance predittiva** — AUC, balanced accuracy delle features selezionate
4. **Parsimonia** — numero di features selezionate

---

## 1. METRICHE DI STABILITÀ

### Pacchetto chiave: `stabm` (CRAN)

- **CRAN:** https://cran.r-project.org/web/packages/stabm/
- **Versione:** aggiornato 2025-07-23
- **Implementa:** 20 misure di stabilità per feature selection
- **Supporta:** Sottoinsiemi di cardinalità diversa tra i resampling runs
- **Feature similarity adjustments:** disponibili per alcune misure

### 1.1 Nogueira Stability Measure (stato dell'arte)

- **Funzione R:** `stabm::stabilityNogueira()`
- **Formula:**

  Dato un insieme di B selezioni, ciascuna un sottoinsieme di {1, ..., p}, con h_j = frequenza di selezione della feature j:

  ```
  Φ̂ = 1 - (1/p · Σ_j [B/(B-1) · ĥ_j/B · (1 - ĥ_j/B)]) / (q̄/(B·p) · (1 - q̄/(B·p)))
  ```

  dove q̄ = media del numero di features selezionate.

- **Range:** [0, 1], dove 1 = stabilità perfetta
- **Vantaggi chiave:**
  - Generalizza Kuncheva index
  - Gestisce sottoinsiemi di dimensione diversa
  - **Primo metodo a permettere intervalli di confidenza e test di ipotesi sulla stabilità**
  - Varianza stimabile → comparazione rigorosa tra metodi
- **Come calcolare CI:**
  ```r
  library(stabm)
  # features_list: lista di B vettori di indici selezionati
  stab_value <- stabilityNogueira(features_list)
  # Per CI: bootstrap sulla lista di selezioni o formula analitica dal paper
  ```
- **Paper:** Nogueira, Sechidis & Brown (2018), JMLR 18(174):1-54
- **URL paper:** https://jmlr.org/papers/v18/17-514.html

### 1.2 Kuncheva Index

- **Funzione R:** `stabm::stabilityKuncheva()` o `OmicsMarkeR::kuncheva()`
- **Formula:** Per due sottoinsiemi A, B di dimensione k da p features:
  ```
  KI(A, B) = (|A ∩ B| · p - k²) / (k · (p - k))
  ```
  Per B selezioni: media su tutte le coppie.
- **Range:** [-1, 1] (originale) o [0, 1] (scalato in OmicsMarkeR)
- **Limitazione:** Richiede sottoinsiemi della **stessa dimensione** (non gestisce |A| ≠ |B|)
- **Paper:** Kuncheva (2007), IASTED International Conference on AI and Applications

### 1.3 Jaccard Index

- **Funzione R:** `stabm::stabilityJaccard()`
- **Formula:** Per due sottoinsiemi A, B:
  ```
  J(A, B) = |A ∩ B| / |A ∪ B|
  ```
  Per B selezioni: media su tutte le coppie.
- **Range:** [0, 1]
- **Limitazione:** Non corretto per chance — set grandi hanno Jaccard alto per caso

### 1.4 Altre Misure in `stabm`

| Misura | Funzione | Note |
|--------|----------|------|
| Dice | `stabilityDice()` | Simile a Jaccard, peso diverso |
| Hamming | `stabilityHamming()` | Basata su disagreement bit-a-bit |
| Ochiai | `stabilityOchiai()` | Coefficiente di similarità |
| Phi | `stabilityPhi()` | Correlazione Phi tra vettori di inclusione |
| Davis | `stabilityDavis()` | - |
| Lustgarten | `stabilityLustgarten()` | - |
| Yu | `stabilityYu()` | - |
| Zucknick | `stabilityZucknick()` | Include feature similarity |
| Intersection count | `stabilityIntersectionCount()` | Semplice conteggio |
| Intersection greedy | `stabilityIntersectionGreedy()` | - |
| Intersection mean | `stabilityIntersectionMean()` | - |
| Intersection MBM | `stabilityIntersectionMBM()` | - |

### 1.5 Average Tanimoto Index

- **Formula:** Equivalente al Jaccard index (Tanimoto e Jaccard sono identici per insiemi binari)
- **Implementazione:** Usare `stabilityJaccard()` da `stabm`

### 1.6 Spearman Correlation of Feature Rankings

- **Non è in `stabm`** — implementazione custom necessaria
- **Implementazione:**
  ```r
  # rankings_matrix: matrice B × p, ogni riga è il ranking delle p features in un resampling run
  spearman_stability <- function(rankings_matrix) {
    B <- nrow(rankings_matrix)
    cors <- combn(B, 2, function(idx) {
      cor(rankings_matrix[idx[1], ], rankings_matrix[idx[2], ], method = "spearman")
    })
    mean(cors)
  }
  ```
- **Range:** [-1, 1]
- **Vantaggi:** Misura la stabilità del ranking completo, non solo del set selezionato

### 1.7 Proportion of Stably Selected Features

- **Implementazione custom:**
  ```r
  # freq: vettore di frequenze di selezione per ogni feature (su B runs)
  # Proporzione di features selezionate in >X% dei runs:
  prop_stable <- function(freq, B, threshold = 0.5) {
    sum(freq / B > threshold) / length(freq)
  }
  # Calcolare per threshold = 0.50, 0.80, 0.95
  ```
- **Output:** Tre valori (prop50, prop80, prop95) — descrivono il "core" stabile

---

## 2. METRICHE DI ACCURATEZZA DELLA SELEZIONE (con ground truth)

Applicabili solo ai dati simulati dove conosciamo le features "vere".

### 2.1 True Positive Rate (TPR / Sensitivity / Recall)

```r
TPR <- sum(selected %in% true_features) / length(true_features)
```

### 2.2 False Discovery Rate (FDR)

```r
FDR <- sum(!(selected %in% true_features)) / max(length(selected), 1)
```

### 2.3 F1 Score

```r
precision <- 1 - FDR
recall <- TPR
F1 <- 2 * precision * recall / (precision + recall)
```

### 2.4 Matthew's Correlation Coefficient (MCC)

- Più robusto di F1 per classi sbilanciate (pochi veri positivi su molte features)
```r
# tp, fp, tn, fn calcolati dal confronto selected vs true_features
MCC <- (tp * tn - fp * fn) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
```

---

## 3. METRICHE DI PERFORMANCE PREDITTIVA

Valutano se le features selezionate producono un buon classificatore.

### 3.1 AUC (Area Under ROC Curve)

```r
library(pROC)
# Dopo aver selezionato features e trainato un classificatore:
roc_obj <- roc(y_test, pred_probs)
auc_value <- auc(roc_obj)
```

### 3.2 Balanced Accuracy

```r
bal_acc <- (sensitivity + specificity) / 2
```

### 3.3 Protocollo di valutazione

**IMPORTANTE:** La performance predittiva deve essere valutata su un test set **completamente indipendente** dalla feature selection:

```
Per ogni resampling run b in 1:B:
  1. Split stratificato per classe: training (70%) + test (30%)
  2. Feature selection su training ONLY
  3. Train classificatore (es. logistic regression) su training con features selezionate
  4. Predict su test set
  5. Calcola AUC e balanced accuracy
```

**Note implementative (2026-03-18):**
- Lo split train/test è stratificato per classe (garantisce entrambe le classi in entrambi i set)
- Se il numero di features selezionate supera n/3, si applica pruning basato su `feature_priority` (importance media o frequenza di selezione)
- Guard per dataset degenerati (classe unica, troppo pochi campioni per classe)

---

## 4. METRICHE DI PARSIMONIA

### 4.1 Numero di Features Selezionate

```r
n_selected <- length(selected)
```

### 4.2 Proporzione di Features Selezionate

```r
prop_selected <- length(selected) / p
```

**Interpretazione clinica:** Un metodo che seleziona 200 features su 1000 non è utile per un pannello diagnostico. L'obiettivo clinico tipico è 5-20 biomarker.

---

## 5. FRAMEWORK DI VALUTAZIONE INTEGRATO

### Metriche primarie (da riportare sempre):

| Categoria | Metrica | Perché |
|-----------|---------|--------|
| Stabilità | **Nogueira** | Gold standard, con CI |
| Stabilità | **Jaccard** | Intuitivo, confrontabile |
| Accuratezza | **TPR** | Quante vere features trovate |
| Accuratezza | **FDR** | Quanti falsi positivi |
| Predizione | **AUC** | Performance classificazione |
| Parsimonia | **n_selected** | Dimensione del pannello |

### Metriche secondarie (supplementary):

- Kuncheva, Dice, Phi (stabilità alternativa)
- F1, MCC (accuratezza alternativa)
- Balanced accuracy (predizione)
- Spearman correlation of rankings
- Proportion stably selected (50%, 80%, 95%)

### Visualizzazioni chiave:

1. **Heatmap di stabilità:** metodi × scenari, colore = Nogueira index
2. **TPR vs FDR curves:** per ogni scenario, come varia al variare della soglia/iperparametro
3. **Radar/spider charts:** confronto multi-dimensionale dei metodi (stabilità, TPR, 1-FDR, AUC, 1/n_selected)
4. **UpSet plots:** overlap tra features selezionate da metodi diversi
5. **Stability vs prediction trade-off:** scatterplot Nogueira vs AUC per ogni metodo × scenario

---

## 6. Riferimenti

1. Nogueira, Sechidis & Brown (2018). On the Stability of Feature Selection Algorithms. JMLR 18(174):1-54
2. Kuncheva (2007). A stability index for feature selection. IASTED AIAP
3. Jaccard (1912). The distribution of the flora in the alpine zone. New Phytologist 11(2):37-50
4. Meinshausen & Bühlmann (2010). Stability selection. JRSS-B 72(4):417-473
5. He & Yu (2010). Stable Feature Selection for Biomarker Discovery. Comp Bio Chem 34(4):215-225
