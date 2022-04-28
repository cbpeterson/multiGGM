# multiGGM
## Bayesian inference of multiple Gaussian graphical models

Author: Christine B. Peterson
Contact: cbpeterson@gmail.com

The given Matlab files for Bayesian inference of multiple graphical models 
are associated with the following publication:

Peterson, C., Stingo, F. and Vannucci, M. (2015). Bayesian inference of multiple
Gaussian graphical models. *Journal of the American Statistical Association*.
110(509): 159—174.

These scripts rely on the Matlab code for G-Wishart sampling by Hao Wang associated with the following publication:
Wang, H. and Li, S. (2012). Efficient Gaussian graphical model determination
under G-Wishart prior distributions. *Electronic Journal of Statistics*.
6: 168—198.

Please cite both publications if you use this code. Thanks!


## OVERVIEW OF FILES


## Example_multiple_graphs.m
Basic example of running MCMC sampler and generating results summaries
on a simple setting with 3 groups with identical dependence structure


## MCMC_multiple_graphs.m
Code for running MCMC sampler


## calc_mrf_C.m
Helper function for calculating normalizing constant for MRF prior


## generate_sim1_input.m
Script to generate matrices similar to those used as input to Simulation 1


## fix_matrix.m
Helper function to ensure that the random matrices generated as simulation
input are in fact positive definite

Please also see [here](https://github.com/cbpeterson/scalable_multiGGM) for a more scalable approach.
