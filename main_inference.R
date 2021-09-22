# Work version  #  Juan Vesga 5 Aug 2021
# CCHF MODEL  Births - Afghanistan 2008 - 2014
# SIR age stratified livestock population + with human demographics
# Livestock pop structure from Mauritania
# imm t0 from Bulgaria (high endemic situation)
# Tick transmission amongst livestock driven by temperature 
# Use of temperature: fitting a multiplication factor (A) for temperature to create Rt value. 
# Spillover to humans SEIR: farming and other occupations

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6805735/
##### set up

rm(list = ls()) 

# install.packages("pacman")         # Install pacman package
library("pacman")                  # Load pacman package
p_load(Hmisc,reshape2,zoo,plyr,ggplot2,gmodels,data.table,devtools,fitR,
coda,deSolve,mvtnorm,ISOweek,lubridate,dplyr,zoo,stringr,here,Rcpp,psych,splines,
data.table,broom, psych)
########################################################################################
# Part 0. select country
########################################################################################
country<-"AFG"


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
sourceCpp(here("src","compute_model_arma.cpp"))
#R vectorized options
# source(here("src","get_objective_vectorized.R"))
# source(here("src","goveqs_basis.R"))


# Set MCMC important parameters 

chain<-"chain1.csv"
n_iterations<-5000


###########################################################################################################
# PART 2. INFERENCE: Model: prior, posterior and inference
###########################################################################################################
theta <- data.frame(
  A=0.08, # driving temperature dependent force of infection
  F_risk=4.92 , # risk for farmers
  O_factor=0.99,
  imm_p=0.69,
  RRreport=0.4,
  knot1=140/30,
  knot2=300/30,
  beta1=0.4)#
  
# 2.1. MY PRIOR
###########################################################################################################
my_prior <- function(theta) {
  ## uniform prior on A
  log.prior.A <- dunif(theta[["A"]], min = 0, max = 1, log = TRUE) 
  ## uniform prior on risk to farmers infected
  log.prior.F_risk <- dunif(theta[["F_risk"]],  min = 0, max = 5, log = TRUE) #0.01 / 5
  ## uniform prior on risk to others: scalar for risk to other occupations
  log.prior.O_factor <- dunif(theta[["O_factor"]], min = 0, max= 1, log = TRUE)
  log.prior.imm_p <- dunif(theta[["imm_p"]], min = 0, max= 1, log = TRUE)
  log.prior.RRreport <- dunif(theta[["RRreport"]], min = 0, max= 1, log = TRUE)
  log.prior.knot1 <- dunif(theta[["knot1"]], min = 1, max= 9, log = TRUE)
  log.prior.knot2 <- dunif(theta[["knot2"]], min = 9.1, max= 14, log = TRUE)
  log.prior.beta1 <- dunif(theta[["beta1"]], min = 0, max= 1, log = TRUE)
  
  return(log.prior.A  + log.prior.F_risk + log.prior.O_factor +
           log.prior.imm_p  + log.prior.RRreport + log.prior.knot1+
           log.prior.knot2+log.prior.beta1)
}
my_prior(theta)

# 2.2. MY POSTERIOR
###########################################################################################################
my_posterior <- function(theta) { 
  
 
  theta_in<-bind_rows(theta)  
  # 1. Run model and get output
  #########################################################################################################  
   
  sim<-get_sim_results(theta_in)
  
  # 3. Estimate LogLikelihood and output posterior density 
  #########################################################################################################  
  loglik_human_reported_month<-  sum(dpois(x = round(observations$cases_human_mo),
                                       lambda=as.numeric(sim$h_inc_month[observations$index_mo_cases]),
                                       log=TRUE),
                               na.rm=TRUE)
  
  loglik_human_reported_year<-  sum(dpois(x = round(observations$cases_human_yr),
                                           lambda=as.numeric(sim$h_inc_year [observations$index_yr_cases]),
                                           log=TRUE),
                                     na.rm=TRUE)
  

  
    
  loglik_liv_prev_age <- sum(dbinom(x= observations$prev_liv_age_pos_IgG, size=observations$prev_liv_age_denom,
                                    prob=sim$l_prev_age, log=TRUE), na.rm=TRUE)
  
  loglik_liv_prev_all <- dbinom(x=observations$prev_liv_all_pos_IgG, size=observations$prev_liv_all_denom,
                                prob=sim$l_prev_all , log=TRUE)
  
  loglik_prev_farmer <- dbinom(x=observations$prev_farmer_Posa, size=observations$prev_farmer_denoma,
                               prob=sim$h_prev_farmer, log=TRUE)
  
  loglik_prev_other <- dbinom(x=observations$prev_other_Posa, size=observations$prev_other_denoma,
                              prob=sim$h_prev_other, log=TRUE)
  
  
  if (country=="AFG"){
  
  log.likelihood <- 
    loglik_human_reported_month +
    loglik_human_reported_year
    loglik_liv_prev_age + 
    loglik_liv_prev_all + 
    loglik_prev_farmer +  
    loglik_prev_other
}else{
  
  log.likelihood <-  
    loglik_liv_prev_age + 
    loglik_liv_prev_all + 
    loglik_prev_farmer +  
    loglik_prev_other
  }
  log.prior <- my_prior(theta)
  
  log.posterior <- log.prior + log.likelihood 
  
  return(log.density=log.posterior)
  
}


 #  profvis({   
 # # 
   my_posterior(theta)
 #  })
###########################################################################################################
#  PART 3. Inference MCMC-MH
###########################################################################################################
init.theta <-c(  
    A=0.54 , # driving temperature dependent force of infection
    F_risk=0.05 , # risk for farmers
    O_factor=0.33 ,
    imm_p=0.97 ,
    RRreport=0.86 ,
    knot1=4.49 ,
    knot2=12.23,
    beta1= 0.57)#
     
proposal.sd <- init.theta/12

n.iterations <- n_iterations
print.info.every <- 50

limits=list(lower=c(A= 0, 
                    F_risk=0,
                    O_factor=0,
                    imm_p=0,
                    RRreport=0,
                    knot1=1, 
                    knot2=9.1,
                    beta1=0),
            
            upper=c(A= Inf, 
                    F_risk=Inf, 
                    O_factor=Inf, 
                    imm_p=1, 
                    RRreport=1,
                    knot1=9, 
                    knot2=14,
                    beta1=Inf)) #inf!

adapt.size.start <- 500
adapt.size.cooling <- 0.999
adapt.shape.start <- 1000
acceptance.rate.weight <- NULL 
acceptance.window <- NULL
#max.scaling.sd <- 5
covmat <- NULL

# 3.1. Chains
##############################################################################################################
start_time <- Sys.time()


res_mcmc <- mcmcMH(target=my_posterior, 
                   init.theta=init.theta, 
                   proposal.sd=proposal.sd, 
                   n.iterations=n.iterations, 
                   limits=limits, 
                   adapt.size.start=adapt.size.start, 
                   adapt.size.cooling=adapt.size.cooling, 
                   adapt.shape.start=adapt.shape.start, 
                   print.info.every=print.info.every,
                   covmat=covmat,
                   #max.scaling.sd= max.scaling.sd,
                   verbose=FALSE) 
end_time <- Sys.time()
timelag<-end_time-start_time
##############################################################################################################
# 3.2. Export results ########################################################################################
##############################################################################################################
vc <- res_mcmc$trace
trace <- apply(vc, 2, unlist)
trace <- data.frame(trace)
acc_rate<- res_mcmc$acceptance.rate
cov_mat <- res_mcmc$covmat.empirical

write.csv(trace, here("output",country,paste("trace_",chain,sep = "")))
write.csv(cov_mat, here("output",country,paste("cov_mat_",chain,sep = "")))
write.csv(acc_rate,here("output",country,paste("acc_rate_",chain, sep = "")))

#########################################################################################################
### END OF CODE
#########################################################################################################



