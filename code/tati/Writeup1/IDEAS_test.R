library(ideas)
library(foreach)
library(doRNG)
RNGkind("L'Ecuyer-CMRG")
library(doParallel)
registerDoParallel(cores=6)

simu_data_rds = "sim_data_ncase_10_nctrl_10_ncell_120_fold_mean_1.2_var_1.5.rds"
sim_data = readRDS("~/Documents/subject-de/git/subject-de/code/tati/Writeup1/sim_data_ncase_10_nctrl_10_ncell_120_fold_mean_1.2_var_1.5.rds")

names(sim_data)
sim_data$count_matrix[1:5,1:5]
head(sim_data$meta_cell)
head(sim_data$meta_ind)

count_matrix = sim_data$count_matrix[1:100,]
meta_cell    = sim_data$meta_cell
meta_ind     = sim_data$meta_ind

var2test      = "phenotype"
var2adjust    = "RIN"
var2test_type = "binary"
var_per_cell  = "cell_rd"

dist1 = ideas_dist(count_matrix, meta_cell, meta_ind, 
                   var_per_cell, var2test, var2test_type, 
                   d_metric = "Was", fit_method = "nb")

pval_ideas = permanova(dist1, meta_ind, var2test, var2adjust, 
                       var2test_type, n_perm=999, r.seed=903)