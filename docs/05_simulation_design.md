# Design Sperimentale della Simulazione

> Ultimo aggiornamento: 2026-03-18
> Status: DRAFT

---

## 1. Filosofia della Simulazione

La simulazione deve essere **empirically-grounded**: parametri estratti da dati reali, non scelti arbitrariamente. Ogni scelta deve essere giustificabile con evidenza dalla letteratura o da analisi esplorativa dei dataset pubblici.

**Approccio a due livelli:**
1. **Livello 1 (Empirico):** Estrarre parametri da dataset reali (distribuzione marginale, struttura correlazione, missing pattern)
2. **Livello 2 (Simulazione):** Generare dati sintetici con parametri controllati, variando sistematicamente uno scenario alla volta

---

## 2. Estrazione Parametri da Dati Reali (Script 03)

### 2.1 Parametri da estrarre per ogni dataset reale

| Parametro | Come stimare | Pacchetto R |
|-----------|-------------|-------------|
| Distribuzione marginale | Fit log-normale, gamma, normale per ogni metabolita | `fitdistrplus` |
| Struttura di correlazione | Matrice di correlazione shrinkage (Ledoit-Wolf) | `corpcor::cor.shrink()` |
| Precision matrix | Ridge-penalized estimation | `rags2ridges::ridgeP()` |
| Proporzione missing | Per-feature e per-sample | Base R |
| Pattern missing | Test MCAR/MAR/MNAR | Custom + `naniar` |
| Effect size empirici | Fold-change e Cohen's d per ogni metabolita caso vs controllo | Base R |
| Numero features | p per dataset | Base R |
| Rapporto p/n | p / n_samples | Base R |

### 2.2 Stima covarianza per p >> n

**Problema:** La matrice di covarianza campionaria è singolare quando p > n. Soluzioni:

1. **Shrinkage estimator (Ledoit-Wolf):**
   ```r
   library(corpcor)
   Sigma_shrink <- cor.shrink(X)  # Shrinkage verso matrice diagonale
   # Converte in covarianza:
   Sigma <- diag(sd_vec) %*% Sigma_shrink %*% diag(sd_vec)
   ```

2. **Ridge estimation della precision matrix:**
   ```r
   library(rags2ridges)
   P_ridge <- ridgeP(S = cov(X), lambda = optPenalty.LOOCV(X)$optLambda)
   Sigma_ridge <- solve(P_ridge)  # Inversione per ottenere covarianza
   ```

3. **Factor model:** Stima con pochi fattori latenti + rumore diagonale
   ```r
   # PCA-based factor model:
   pca <- prcomp(X, center = TRUE, scale. = TRUE)
   k <- 10  # Numero di fattori
   L <- pca$rotation[, 1:k] %*% diag(pca$sdev[1:k])
   Sigma_factor <- L %*% t(L) + diag(diag(var(X) - L %*% t(L)))
   ```

**Raccomandazione:** Usare shrinkage (Ledoit-Wolf) come default; ridge estimation come sensibilità check.

---

## 3. Framework di Simulazione (Script 04)

### 3.1 Generazione Dati Base

**Step 1: Generare matrice Z multivariata normale**
```r
library(mvtnorm)
Z <- rmvnorm(n, mean = rep(0, p), sigma = Sigma)
```

**Step 2: Trasformare in distribuzione log-normale (tipica per metaboliti)**
```r
# Z è N(0, Sigma) → exp(Z) è log-normale con struttura di correlazione approssimata
X <- exp(Z * sd_vec + mean_vec)
# Dove mean_vec e sd_vec sono parametri empirici dal dataset reale
```

**Step 3: Impiantare segnale (shift nei "veri" metaboliti)**
```r
# Definire S = indici dei metaboliti "veri" (con signal)
# Per i campioni del gruppo "caso":
delta <- log(FC_vector)  # FC_vector contiene i fold-change desiderati
Z_case[, S] <- Z_case[, S] + delta  # Shift in scala log
X_case <- exp(Z_case * sd_vec + mean_vec)
```

**Nota critica sull'impianto del segnale:**
Il segnale deve essere impiantato PRIMA della trasformazione esponenziale, in scala log, per preservare la struttura di correlazione. Impiantare dopo la trasformazione altererebbe le correlazioni.

**Step 4: Aggiungere missing values (MNAR)**
```r
# LOD-based MNAR: valori sotto il k-esimo percentile hanno probabilità di essere missing
for (j in 1:p) {
  lod <- quantile(X[, j], probs = 0.05, na.rm = TRUE)
  prob_missing <- ifelse(X[, j] < lod, 0.8, 0.02)  # Alta prob sotto LOD
  missing_idx <- which(runif(n) < prob_missing)
  X[missing_idx, j] <- NA
}
```

### 3.2 Parametri del Segnale Impiantato

| Parametro | Range | Giustificazione |
|-----------|-------|-----------------|
| Numero vere features (p_true) | 10, 20, 50 | Pannelli biomarker realistici |
| Effect size (FC) | 1.2, 1.5, 2.0, 3.0 | FC < 1.5 dominante in letteratura |
| Distribuzione FC tra le vere features | Uniforme o decrescente | Pochi metaboliti con FC grande, molti con FC piccolo |
| Direzione | Mix up/down | Realistico: metaboliti sia up che down-regolati |

---

## 4. Scenari di Simulazione

### Scenario Base (default)
- n = 100 (50 caso + 50 controllo)
- p = 500
- p_true = 20 (4% del totale)
- FC = 1.5 per le vere features
- Correlazione: struttura empirica da dataset reale
- Missing: 5% MNAR
- Preprocessing: log-transform + auto-scaling

### Scenari Variati (uno alla volta dal base)

| # | Scenario | Variazione | Valori | Evidenza |
|---|----------|-----------|--------|----------|
| S1 | **p/n ratio** | Varia p e n | p/n = 5, 10, 20, 50 | Range realistico studi pubblicati |
| S2 | **Effect size** | Varia FC | 1.2, 1.5, 2.0, 3.0 | FC < 1.5 dominante (Metabolomics 2019) |
| S3 | **Multicollinearità** | Scala correlazioni + tipo struttura | r_max = 0.3, 0.5, 0.7, 0.9; `correlation_source`: empirical (default), ar1, block | NMR: r fino a 0.9; LC-MS: r tipico 0.3-0.6. AR(1) e block per sensitivity |
| S4 | **Non-linearità** | Aggiungi interazioni | 0, 5, 10 interazioni | Tree-based dovrebbero catturarle, lineari no |
| S5 | **Confounders** | Aggiungi covariate | 0, 1, 3 confounders | Età, BMI, farmaci (letteratura) |
| S6 | **Missing data** | Varia % missing | 0%, 5%, 15%, 30% | 5-30% tipico in metabolomica |
| S7 | **Preprocessing** | Varia normalizzazione | log, pareto, auto-scaling, PQN | Impatto documentato (Metabolites 2022) |

**Totale scenari:** 7 × ~4 livelli ciascuno = ~28 condizioni sperimentali

### 4.1 Dettaglio Scenario S1: p/n ratio

| Condizione | n | p | p/n | p_true |
|-----------|---|---|-----|--------|
| S1a | 100 | 500 | 5 | 20 |
| S1b | 50 | 500 | 10 | 20 |
| S1c | 50 | 1000 | 20 | 20 |
| S1d | 20 | 1000 | 50 | 20 |

### 4.2 Dettaglio Scenario S4: Non-linearità

```r
# Aggiungere interazioni al modello generativo:
# logit(P(Y=1)) = Xβ + X_i * X_j * γ_ij
# Dove (i,j) sono coppie di features interagenti
n_interactions <- c(0, 5, 10)
interaction_strength <- 0.5  # Forza dell'interazione
```

### 4.3 Dettaglio Scenario S5: Confounders

```r
# Generare covariate confondenti:
age <- rnorm(n, mean = 55, sd = 10)
bmi <- rnorm(n, mean = 27, sd = 4)
# Confounders influenzano sia Y che X:
# P(Y=1) dipende anche da age/bmi
# X[, subset] dipende anche da age/bmi
confounder_effect_on_X <- 0.3  # Coefficiente
confounder_effect_on_Y <- 0.5
```

---

## 5. Bootstrap Resampling per Stabilità

### Protocollo

```
Per ogni scenario s in S:
  Per ogni replicazione r in 1:R (R = 50 dataset simulati per scenario):
    Genera dataset D_r con parametri dello scenario s
    Per ogni bootstrap b in 1:B (B = 200 resamples):
      Genera D_r_b = bootstrap sample di D_r (o subsampling 63.2%)
      Per ogni metodo m in M:
        Applica m a D_r_b → set selezionato S_{r,b,m}
    Calcola metriche di stabilità per il metodo m sulla replicazione r:
      - Nogueira index su {S_{r,1,m}, ..., S_{r,B,m}}
      - Jaccard medio
      - TPR, FDR (confronto con ground truth)
      - AUC su test set
  Aggrega metriche su R replicazioni → media ± CI
```

### Parametri del Resampling

| Parametro | Valore | Giustificazione |
|-----------|--------|-----------------|
| B (bootstrap) | 200 | Compromesso tra stabilità stimata e tempo; Nogueira (2018) usa 100-200 |
| R (replicazioni) | 50 | Per ottenere CI robusti sulle metriche aggregate |
| Tipo resampling | Subsampling senza replacement (63.2%) | Raccomandato da Meinshausen & Bühlmann; riduce bias rispetto a bootstrap con replacement |
| Stratificazione | Sì, per gruppo caso/controllo | Mantiene proporzione gruppi |

### Alternativa: Subsampling vs Bootstrap

- **Subsampling (0.632):** Più teoricamente giustificato per stabilità; usato da stability selection
- **Bootstrap (con replacement):** Più comune nella pratica; duplicati possono creare bias
- **Raccomandazione:** Usare **subsampling** come default, bootstrap come sensitivity check

---

## 6. Imputazione Missing Values

### Metodi da confrontare (nello Scenario S6)

| Metodo | Implementazione R | Tipo |
|--------|------------------|------|
| Half-minimum | `min(x, na.rm=T) / 2` | Semplice, comune |
| KNN imputation | `impute::impute.knn()` | MAR-oriented |
| Random Forest | `missForest::missForest()` | Non-parametrico |
| Probabilistic minimum | `imputeLCMD::impute.MinProb()` | MNAR-specific |

### Protocollo nello scenario S6
Per ogni livello di missing (5%, 15%, 30%):
1. Imputare con metodo default (half-minimum)
2. Applicare tutti i 12 metodi di feature selection
3. **Sensitivity:** ripetere con KNN e MinProb per verificare se il metodo di imputazione cambia le conclusioni

---

## 7. Preprocessing Pipeline (Scenario S7)

### Opzioni di Normalizzazione

| Metodo | Formula | Quando usare |
|--------|---------|-------------|
| Log transformation | x → log2(x + 1) | Default per metabolomica |
| PQN (Probabilistic Quotient) | Normalizza al metabolita di riferimento mediano | NMR, robusto a outlier |
| TIC (Total Ion Current) | x_ij / Σ_j x_ij | LC-MS, normalizza per intensità totale |
| Median normalization | x_ij / median_j(x_ij) | Semplice, robusto |

### Opzioni di Scaling

| Metodo | Formula | Effetto |
|--------|---------|--------|
| Auto-scaling | (x - mean) / sd | Ogni feature ha varianza 1 |
| Pareto scaling | (x - mean) / sqrt(sd) | Riduce importanza features ad alta varianza senza equalizzare |
| Range scaling | (x - min) / (max - min) | Normalizza a [0, 1] |

### Protocollo nello scenario S7
Combinazioni da testare:
- log + auto-scaling (default)
- log + pareto
- PQN + auto-scaling
- nessun preprocessing (raw)

---

## 8. Seed e Riproducibilità

```r
# Master seed:
set.seed(42)

# Generare seed per ogni combinazione scenario × replicazione:
all_seeds <- sample.int(.Machine$integer.max, size = n_scenarios * n_replications)

# Ogni simulazione usa il suo seed dedicato → riproducibilità completa
```

---

## 9. Struttura Output della Simulazione

```
data/
├── simulated/
│   ├── scenario_S1a/
│   │   ├── rep_001/
│   │   │   ├── X.rds          # Matrice features × samples
│   │   │   ├── y.rds          # Vettore outcome
│   │   │   ├── true_features.rds  # Indici features vere
│   │   │   └── params.rds     # Parametri usati
│   │   ├── rep_002/
│   │   └── ...
│   ├── scenario_S1b/
│   └── ...
```

---

## 10. Dimensionamento dello Studio

| Componente | Quantità |
|-----------|----------|
| Scenari | 7 categorie × ~4 livelli = 28 |
| Replicazioni per scenario | 50 |
| Bootstrap per replicazione | 200 |
| Metodi | 12 |
| **Totale fit di feature selection** | **28 × 50 × 200 × 12 = 3,360,000** |

→ Necessaria parallelizzazione massiva (vedi [06_computational_plan.md](06_computational_plan.md))
