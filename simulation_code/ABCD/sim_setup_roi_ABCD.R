# This file loads the ABCD structural data and sets up for the simulations

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

# Read in ROI ABCD data
# harmonized structural data
dat = as.data.frame(read.csv("/media/disk2/multivariateBWAS/ABCD_sims/ABCD_struct_harmonized.csv"))


# Use only total GMV for the 68 ROIs for X
X = dat[,grep("Vol", colnames(dat))]
# View missingness pattern
# mice::md.pattern(X)

# Phenotypes we want to predict: age, fluid intelligence, depression score
dat %<>% select(age_years, nihtbx_fluidcomp_fc, cbcl_scr_dsm5_depress_t)
# View missingness pattern
# mice::md.pattern(dat)

# add imaging data as a matrix to the dataframe
dat$img = as.matrix(X)
rm(X)

# data.frame for simulations
colnames(dat) = c("age", "fluid", "depression", "img")
