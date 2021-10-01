# Work version  #  Juan Vesga 5 Aug 2021
# Script to plot moel fits from posterior samples

rm(list = ls()) 
library("pacman")                  # Load pacman package
p_load(Hmisc,reshape2,zoo,plyr,ggplot2,gmodels,data.table,devtools,fitR,
       coda,deSolve,mvtnorm,ISOweek,lubridate,dplyr,zoo,stringr,here,Rcpp,psych,splines,
       data.table,broom,matrixStats)



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



##################################

# Load posterior samples
posteriors <- read.table(here("output",country,"posteriors_combined.txt"), sep = "\t" )
theta<-posteriors%>%
  select(A,F_risk, O_factor,imm_p,RRreport,knot1, knot2, beta1)


## Run model 

 sim<-get_sim_results(theta)
# 
# saveRDS(sim, file=here("output",country,"simulations.rds"))

# ... or load previously saved sim results
  sim<-readRDS(here("output",country,"simulations_combined.rds"))
# 
 # pdf(here("output",country,"model_fits_newdata_revised.pdf")) 

plot_fits_function(sim,observations)

 # dev.off() 


 





