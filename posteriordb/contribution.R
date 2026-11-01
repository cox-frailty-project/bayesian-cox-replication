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

#### PosteriorDB Data Object####

## Data 

k <- survival::kidney
k <- k[order(k$time), ]

# Covariates exactly as in the paper: age in years, sex coded 1 = male, 2 = female,
# 0/1 indicators for disease types GN, AN, PKD (Other is the baseline).
X <- cbind(age = k$age,
           sex = k$sex,
           GN  = as.integer(k$disease == "GN"),
           AN  = as.integer(k$disease == "AN"),
           PKD = as.integer(k$disease == "PKD"))

stan_data <- list(N   = nrow(k),
                  M   = length(unique(k$id)),
                  P   = ncol(X),
                  t   = k$time,
                  D   = k$status,
                  pid = k$id,
                  X   = X,
                  sigma = 1)

## Data info

data_name = "kidney_data"

data_info <- as.pdb_data_info(list(
  name        = data_name,
  title       = "Kidney Survival Data",
  description = "Data on kidney catheter infections and their reaccurance times. From the survival package in r",
  keywords    = c("survival", "kidney", "infection"),
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

stan_path <- "stan/frailty_cox_replication.stan"

stan_code <- paste(
  readLines(stan_path),
  collapse = "\n"
)

cat(stan_code)

## Model info

model_name <- "frailty_cox_replication"

model_info <- as.pdb_model_info(list(
  name = model_name,
  title = "Replication of Cox Frailty Partial Likelihood including Frailties",
  description = "A stan replication of the model by C. A. McGilchrist and C. W. Aisbett (1991) for use with the Kidney dataset from the survival package in r",
  keywords = c("hazard regression", "frailty", "survival analysis", "regression"),
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
  keywords = c("hazard regression", "frailty", "survival analysis", "survival", "kidney", "infection"),
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
  PDB_PATH = "/Users/my_files/Documents/GitHub/bayesian-cox-replication/posteriordb/paper_replication"
)

pdbl <- pdb_local()

write_pdb(pdb_data, pdbl, overwrite = TRUE)

write_pdb(pdb_model, pdbl, overwrite = TRUE)

write_pdb(pdb_posterior, pdbl, overwrite = TRUE)
