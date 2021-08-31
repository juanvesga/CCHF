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

########################################################################################
# Part 1. Import set up: Load data and model parameters and functions
########################################################################################
source(here("src","setup_data.R"))
source(here("src","setup_model.R"))
source(here("src","make_model.R"))
source(here("src","get_sim_results.R"))
#Rcpp vectorized functions
source(here("src","get_objective_vectorized.R"))
source(here("src","goveqs_basis_rcpp.R"))
source(here("src","plot_fits_function.R"))
sourceCpp(here("src","compute_model_arma.cpp"))
#R vectorized options
# source(here("src","get_objective_vectorized.R"))
# source(here("src","goveqs_basis.R"))

country<- "AFG"

##################################

# Load posterior samples
posteriors <- read.table(here("output",country,"posteriors_2008_2014.txt"), sep = "\t" )
theta_in<-posteriors%>%
  select(A,F_risk, O_factor)


## Run model 
sim<-get_sim_results(theta_in)
saveRDS(sim, file=here("output",country,"simulations.rds"))

# ... or load previously saved sim results
sim<-readRDS(here("output",country,"simulations.rds"))

pdf(here("output",country,"model_fits.pdf")) 

plot_fits_function(sim,observations)

dev.off() 


