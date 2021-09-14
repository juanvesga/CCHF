# Work version  #  Juan Vesga 5 Aug 2021
# Script to plot moel fits from posterior samples

rm(list = ls()) 

library(Hmisc)
library(reshape2)
library(zoo)
library(plyr)
library(ggplot2)
library(gmodels)
library(data.table)
library(devtools)
library(fitR)
library(coda)
library(deSolve)
library(mvtnorm)
library(ISOweek)
library(lubridate)
library(dplyr)
library(zoo)
library(stringr)
library(here)
library (Rcpp)
library(matrixStats)
library(splines)
library(data.table)
library(broom)

country<- "AFG"


########################################################################################
# Part 1. Import set up: Load data and model parameters and functions
########################################################################################
source(here("src","setup_data.R"))
source(here("src","setup_model.R"))
source(here("src","make_model.R"))
source(here("src","get_sim_results.R"))
source(here("src","spline_functions.R"))
#Rcpp vectorized functions
source(here("src","get_objective_vectorized.R"))
source(here("src","goveqs_basis_rcpp.R"))
source(here("src","plot_fits_function.R"))
sourceCpp(here("src","compute_model_arma.cpp"))
#R vectorized options
# source(here("src","get_objective_vectorized.R"))
# source(here("src","goveqs_basis.R"))


run_existing<-1

##################################

# Load posterior samples
posteriors <- read.table(here("output",country,"posteriors.txt"), sep = "\t" )
theta<-posteriors%>%
  select(A,F_risk, O_factor,imm_p, knot1, knot2, beta1)


if (run_existing==0){
## Run model 

sim<-get_sim_results(theta)

saveRDS(sim, file=here("output",country,"simulations.rds"))
}else{
# ... or load previously saved sim results
sim<-readRDS(here("output",country,"simulations.rds"))
}

# Plot 


  # pdf(here("output",country,"model_fits.pdf")) 

plot_fits_function(sim,observations)

   # dev.off() 






