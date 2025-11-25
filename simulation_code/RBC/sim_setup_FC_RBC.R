# Load required libraries
library(parallel)
library(progress)
library(pbmcapply)
library(logspline)
library(e1071)
library(ranger)
library(knitr)
library(kableExtra)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(purrr)
library(glmnet)
library(patchwork)
library(magrittr)

# plot settings
theme_set(theme_bw())
theme_update(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# source functions
source("/media/disk2/multivariateBWAS/code_for_review/convenience_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/influence_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/corr_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/build_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/results_funcs.R")
source("/media/disk2/multivariateBWAS/code_for_review/nice_figures.R")
source("/media/disk2/multivariateBWAS/code_for_review/sim_funcs_3param.R")
source("/media/disk2/multivariateBWAS/code_for_review/sim_funcs_4param.R")

# estimators
all_types = c("pf", "pm", "sf", "sm", "ct", "pl", "sl", "pq", "sq", "cq", "pz", "sz", "pfr", "pmr", "plr", "slr", "pqr", "sqr", "pzr", "szr")

# harmonized functional data
dat = as.data.frame(read.csv("/media/disk2/RBC_version0.1/FC_network17_harmonized_covariates.csv"))


# FC variables
X = dat[,grep("_to_", colnames(dat), value = TRUE)]

# View missingness pattern
# mice::md.pattern(X)

# Phenotypes we want to predict: age, internalizing, externalizing
dat %<>% select(study, study_site, age, internalizing_mcelroy_harmonized_all_samples,
                externalizing_mcelroy_harmonized_all_samples, p_factor_mcelroy_harmonized_all_samples, attention_mcelroy_harmonized_all_samples)
# View missingness pattern
# mice::md.pattern(dat)

# add imaging data as a matrix to the dataframe
dat$img = as.matrix(X)
rm(X)

# data.frame for simulations
colnames(dat) = c("study", "study_site", "age", "internal", "external", "p_factor", "attention", "img")
