# bayesian-cox-replication

A fully Bayesian replication of McGilchrist & Aisbett (1991), *Regression with Frailty in Survival Analysis*, implemented in Stan. The paper fits a Cox proportional hazards model with a shared log-normal frailty per patient to recurrence times of infection in kidney dialysis patients. We reproduce the regression coefficients, frailty variance and patient-level frailties, and aim to contribute the resulting posterior to [posteriorDB](https://github.com/stan-dev/posteriordb).

Project for *Bayesian Statistics and Data Analysis* (Uppsala University). See the [course project instructions](https://github.com/MansMeg/BSDA/blob/main/project/project_work.pdf).

## Links

- **Report (Overleaf, WIP):** https://www.overleaf.com/project/6ab2864eaeb612620349f6de
- **Paper:** McGilchrist & Aisbett (1991), *Biometrics* 47(2), 461–466. [JSTOR](https://www.jstor.org/stable/2532138) · [doi:10.2307/2532138](https://doi.org/10.2307/2532138)
- **Course project repo:** https://github.com/MansMeg/BSDA/
- **Contributing to posteriorDB:** https://stan-dev.r-universe.dev/posteriordb/doc/contributing.html

## Data

The `kidney` dataset from the R `survival` package: 38 patients, 76 recurrence intervals, with covariates age, sex and disease type (GN, AN, PKD, with Other as the baseline). Age is standardised before fitting.

## Model

Cox partial likelihood with Breslow handling of ties and a shared frailty per patient:

- `h_i(t) = h_0(t) · exp(xᵢᵀβ + u_{pid(i)})`
- `β ~ Normal(0, τ²)` with fixed τ = 2.5
- `u ~ Normal(0, σ²)`, constrained to sum to zero
- `σ ~ Half-Normal(0, 1)`
- Frailties reported as `Z = exp(u)`, matching the paper's Table 1

## Reproducing

Requires R with `rstan` and `survival`.

```r
install.packages(c("rstan", "survival", "rmarkdown"))
rmarkdown::render("analysis/model_run.Rmd")
```

Fitting takes about a minute with four chains.

## Repository layout

```
stan/          Stan model
analysis/      R Markdown notebook that fits the model and compares to the paper
report/        Written report
```

## Authors

Thauri Dattadeen
Isac Sandberg Kanon
Julius Searancke

## References

McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis. *Biometrics*, 47(2), 461–466.
