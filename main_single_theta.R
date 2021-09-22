# Work version  #  Juan Vesga 5 Aug 2021
# Script to run a single set of theta values against observed data 

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
library(profvis)
library(matrixStats)
library(splines)
library(data.table)
library(broom)
library(psych)

country<- "AFG"
test_mode=1 # Set to 1 to run test parameter values

########################################################################################
# Part 1. Import set up: Load data and model parameters and functions
########################################################################################
# All required functions
source(here("src","setup_data.R"))
source(here("src","setup_model.R"))
source(here("src","get_sim_results.R"))
source(here("src","plot_fits_function.R"))
source(here("src","spline_functions.R"))



#Rcpp vectorized functions
source(here("src","make_model.R"))
source(here("src","get_objective_vectorized.R"))
source(here("src","goveqs_basis_rcpp.R"))
sourceCpp(here("src","compute_model_arma.cpp"))

#R vectorized options
# source(here("src","make_model.R"))
# source(here("src","get_objective_vectorized.R"))
# source(here("src","goveqs_basis.R"))

#R yearly age ODEs )not vectorized)
# source(here("src","get_objective.R"))
# source(here("src","seir_model_ageing.R"))


###########################################################################################################
# PART 1. Define theta 
###########################################################################################################
if (test_mode==1){
  theta <- data.frame(
    A=0.5546031, # driving temperature dependent force of infection
    F_risk=0.1, # risk for farmers
    O_factor=0.1,
    imm_p=0.649,
    animal_in=20000,
    knot1=140/30,
    knot2=300/30,
    beta1=0.4)#
    
}else{

# Load posterior samples
posteriors <- read.table(here("output",country,"posteriors.txt"), sep = "\t" )

theta_in<-posteriors%>%
  select(A,F_risk, O_factor, imm_p,RRreport)

theta <- data.frame(
  A= median(theta_in$A)	, # driving temperature dependent force of infection
  F_risk=median(theta_in$F_risk)	, # risk for farmers
  O_factor=median(theta_in$O_factor),
  imm_p=median(theta_in$imm_p),
  RRreport=median(theta_in$RRreport),
  knot1=median(knot1$RRreport),
  knot2=median(knot2$RRreport),
  beta1=median(beta1$RRreport)
  )#
}


# profvis({  

sim<-get_sim_results(theta)

 
# })


  # pdf(here("output",country,"model_proj_influx_inf.pdf")) 

plot_fits_function(sim,observations)


  # dev.off() 









