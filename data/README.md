# data

Reference values transcribed from McGilchrist & Aisbett (1991), *Biometrics* 47(2), 461–466.

| File | Contents |
| :--- | :--- |
| `paper_table2.csv` | Regression coefficient estimates and standard errors for age, sex, GN, AN and PKD, plus the frailty variance σ² = 0.3821. |

The kidney infection data itself is not stored here. It ships with the R `survival` package as the `kidney` dataset and is loaded directly in the analysis:
The `frail` column has been checked against the published Table 1 and matches for all 38 patients, so it is used directly as the reference for the frailty comparison rather than a separate CSV.

```r
library(survival)
data(kidney)
```

Columns: `id` (patient), `time` (recurrence time), `status` (1 = infection, 0 = censored), `age`, `sex` (1 = male, 2 = female), `disease` (Other, GN, AN, PKD), `frail` (the paper's Table 1 frailty estimate for that patient, repeated on both of their rows).
