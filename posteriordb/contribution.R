# This file did not commit in properly, this document is WIP


library(posteriordb)
# Remember to set WD to repo
setwd("~/Documents/GitHub//bayesian-cox-replication/")

## Libraries

library(rstan)
library(survival)
library(tidyverse)
options(mc.cores = parallel::detectCores())

## Reference for PosteriorDBs bibtex file:
"@article{mcgilchrist1991regression,
  author  = {McGilchrist, C. A. and Aisbett, C. W.},
  title   = {Regression with Frailty in Survival Analysis},
  journal = {Biometrics},
  year    = {1991},
  volume  = {47},
  number  = {2},
  pages   = {461--466},
  doi     = {10.2307/2532138}
}"

#### PosteriorDB Data Object (Same For Both) ####

## Data 

## Importing it
kidney
k <- kidney

## Covariates defined to match the paper: age in years, sex coded 1 = male, 2 = female,
## 0/1 indicators for disease types GN, AN, PKD ("Other" is the baseline).
X <- cbind(age = k$age,
           sex = k$sex,
           GN  = as.integer(k$disease == "GN"),
           AN  = as.integer(k$disease == "AN"),
           PKD = as.integer(k$disease == "PKD"))

## List format to align with the stan model
stan_data <- list(N = length(k$id),
                  M = length(unique(k$id)),
                  P = ncol(X),
                  t = k$time,
                  C = k$status,
                  pid = k$id,
                  X = X)

## Data info

data_name = "kidney_data"

data_info <- as.pdb_data_info(list(
  name        = data_name,
  title       = "Kidney Survival Data",
  description = "Data on kidney catheter infections and their reaccurance times. From the survival package in r",
  keywords    = c("survival", "kidney", "infection", "weibull"),
  urls        = "https://search.r-project.org/CRAN/refmans/survival/html/kidney.html",
  references  = NULL,
  added_by    = "Isac Kanon, Julius Searancke, Thauri Dattadeen",
  added_date  = Sys.Date()
))


pdb_data <- as.pdb_data(stan_data, info = data_info)
info(pdb_data)


#### Weibull Replication ####

### PosteriorDB model object

## Getting model

stan_path <- "stan/Weibull_v1.stan"

stan_code <- paste(
  readLines(stan_path_wei),
  collapse = "\n"
)

cat(stan_code)

## Model info

model_name <- "weibull_frailty_model"

model_info <- as.pdb_model_info(list(
  name = model_name_wei,
  title = "Weibull estimation for frailty term hazard regression",
  description = "Weibull estimation of the frailty hazard regression formulated by C. A. McGilchrist and C. W. Aisbett (1991) for use with the Kidney dataset from the survival package in r",
  keywords = c("hazard regression", "frailty", "weibull", "regression"),
  urls = c("https://doi.org/10.2307/2532138",
           "https://search.r-project.org/CRAN/refmans/survival/html/kidney.html"),
  references = "mcgilchrist1991regression", # From bibtex above
  framework = "stan",
  added_by = "Isac Kanon, Julius Searancke, Thauri Dattadeen",
  added_date = Sys.Date()
))


pdb_model <- as.pdb_model_code(stan_code, info = model_info, framework = "stan")
info(pdb_model)


### PosteriorDB posterior object

pdb_posterior <- as.pdb_posterior(list(
  pdb_model_code = pdb_model,
  pdb_data = pdb_data,
  keywords = c("hazard regression", "frailty", "weibull", "regression", "survival", "kidney", "infection"),
  urls = c("https://doi.org/10.2307/2532138",
           "https://search.r-project.org/CRAN/refmans/survival/html/kidney.html"),
  references = "mcgilchrist1991regression",
  dimensions = list(
    shape = 1,
    log_scale_0 = 1,
    betav= as.integer(stan_data$P),
    u = as.integer(stan_data$M),
    sigma = 1,
    tau = 1
  ),
  added_by = "Isac Kanon, Julius Searancke, Thauri Dattadeen",
  added_date = Sys.Date()
))


### PosteriorDB write

Sys.setenv(
  PDB_PATH = "~/posteriordb/weibull"
)

pdbl <- pdb_local()

write_pdb(pdb_data, pdbl, overwrite = TRUE)

write_pdb(pdb_model, pdbl, overwrite = TRUE)

write_pdb(pdb_posterior, pdbl, overwrite = TRUE)
