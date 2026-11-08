# bayesian-cox-replication

A fully Bayesian replication of McGilchrist & Aisbett (1991), *Regression with Frailty in Survival Analysis*, implemented in Stan. The paper fits a Cox proportional hazards model with a shared log-normal frailty per patient to recurrence times of infection in kidney dialysis patients. We reproduce the regression coefficients, frailty variance and patient-level frailties, and aim to contribute the resulting posterior to [posteriorDB](https://github.com/stan-dev/posteriordb).

Project for *Bayesian Statistics and Data Analysis* (Uppsala University). See the [course project instructions](https://github.com/MansMeg/BSDA/blob/main/project/project_work.pdf).

## Links

- **Report:** [`report/Final_project_BSDA.pdf`](report/Final_project_BSDA.pdf) (Overleaf source: https://www.overleaf.com/project/6ab2864eaeb612620349f6de)
- **Paper:** McGilchrist & Aisbett (1991), *Biometrics* 47(2), 461–466. [JSTOR](https://www.jstor.org/stable/2532138) · [doi:10.2307/2532138](https://doi.org/10.2307/2532138)
- **Course project repo:** https://github.com/MansMeg/BSDA/
- **Contributing to posteriorDB:** https://stan-dev.r-universe.dev/posteriordb/doc/contributing.html

## Data

The `kidney` dataset from the R `survival` package: 38 patients, 76 recurrence intervals, with covariates age, sex and disease type (GN, AN, PKD, with Other as the baseline). Covariates are used exactly as in the paper: age in years, sex coded 1 = male / 2 = female, and 0/1 disease indicators. Reference values from the paper are in [`data/`](data/).

## Model

Cox partial likelihood (Breslow handling of ties, observations sorted by time) with a shared log-frailty per patient, matching the setup in the paper:

- `h_i(t) = h_0(t) · exp(xᵢᵀβ + u_{pid(i)})`
- `β`: flat (improper uniform) prior, as in the paper's penalised partial likelihood
- `u ~ Normal(0, σ²)`, constrained to sum to zero, with σ fixed at 1 (passed in as data)
- `σ²` is reported as the plug-in estimate `uᵀu / (M − 1)` in `generated quantities`, mirroring the paper's REML-style update
- Frailties reported as `Z = exp(u)`, matching the paper's Table 1

The model is in [`stan/frailty_cox_replication.stan`](stan/frailty_cox_replication.stan). The notebook fits it twice: with NUTS (posterior means/medians) and with the Newton optimizer (posterior mode / MAP, with Hessian-based standard errors), and compares both to the paper's Tables 1 and 2.

## Reproducing

Requires R with `rstan` and `survival`.

```r
install.packages(c("rstan", "survival", "rmarkdown"))
rmarkdown::render("analysis/paper_replication_model_run.Rmd")
```

Sampling runs 4 chains of 50,000 iterations (thinned by 4), which takes a few minutes; the optimizer run is near-instant.

## Repository layout

```
data/          Reference values transcribed from the paper (Table 2) + data notes
stan/          frailty_cox_replication.stan: the Stan model
analysis/      paper_replication_model_run.Rmd: fits the model (NUTS + optimizer) and
               compares to the paper; frailties_comparison.png is the resulting figure
posteriordb/   PosteriorDB contribution (model, data, posterior info, contribution.R)
report/        Final_project_BSDA.pdf: the written report
```

## Authors

- Thauri Dattadeen
- Julius Searancke
- Isac Sandberg Kanon

## References

McGilchrist, C. A. and Aisbett, C. W. (1991). Regression with frailty in survival analysis. *Biometrics*, 47(2), 461–466.
