# Bayes_Pred
This repository contain all codes to implement the methods introduced in the paper "Asymptotics of predictive distributions driven by sample means and variances" (Garelli, Leisen, Pratelli, Rigo) and replicate the experiments shown in there.

The file Convergence.cpp contains the definition of all functions needed to replicate the experiments of section 4.1 of the paper. The most relevant functions are gauss_path_1 and cop_path_1, which compute the predictive density and c.d.f. at each step of the predictive resampling algorithm for both the Gaussian predictive distribution and the copula-based predictive distribution. The file Convergence.R contains the R script that reproduces the experiments of section 4.1. 
