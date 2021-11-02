#Merge chains  CCHF May 2021

rm(list = ls()) 
library(Hmisc)
library(reshape2)
library(zoo)
library(plyr)
library(ggplot2)
library(gmodels)
library(data.table)
library(devtools)
#install_github("sbfnk/fitR")
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

library(coda)

country<- "AFG"


burn<- 0.5 # at what point in trace to chop burn-in runs
thin_n<- 10 # sample every n runs (reduce autocorrelation) 
# analyse MCMC samples


# Import chains
trace1_1 <- read.csv(here("output",country,"trace_chain4_combinedMCMC.csv"))
trace1_2 <- read.csv(here("output",country,"trace_chain3.csv"))

######################################################################
# 1. chain 1
head(trace1_1)
mcmc.trace1_1 <- mcmc(trace1_1)
summary(mcmc.trace1_1)        
acceptanceRate <- 1 - rejectionRate(mcmc.trace1_1)
acceptanceRate
effectiveSize(mcmc.trace1_1)
mcmc.trace.burned1_1 <- burnAndThin(mcmc.trace1_1, 
                                    burn = round(nrow(trace1_1)*burn))
mcmc.trace.burned.thinned1_1<- 
  burnAndThin(mcmc.trace.burned1_1, thin = thin_n) 

summary(mcmc.trace.burned.thinned1_1)

#################
# 2. chain 2
head(trace1_2)
mcmc.trace1_2 <- mcmc(trace1_2)
summary(mcmc.trace1_2)        
acceptanceRate <- 1 - rejectionRate(mcmc.trace1_2)
acceptanceRate
effectiveSize(mcmc.trace1_2)
mcmc.trace.burned1_2 <- burnAndThin(mcmc.trace1_2, burn = round(nrow(trace1_2)*burn))
mcmc.trace.burned.thinned1_2<- burnAndThin(mcmc.trace.burned1_2, thin = thin_n) 
summary(mcmc.trace.burned.thinned1_2)

#################
# 1.3. Combine Chains 
mcmc1_1 <- data.frame(mcmc.trace.burned.thinned1_1)
mcmc1_2 <- data.frame(mcmc.trace.burned.thinned1_2)
df1_1 <- data.frame(mcmc1_1)
df1_2 <- data.frame(mcmc1_2)
df_both <- rbind(df1_1, df1_2)
trace_both <- mcmc(df_both)
autocorr.plot(df_both)
plot(trace_both)

write.table(df1_1, here("output",country,"posteriors_combMCMC.txt"), sep = "\t" )
#####################################################################################
windows()
plot(mcmc.trace1_1)
plot(mcmc.trace.burned1_1)
autocorr.plot(mcmc.trace.burned1_1)
autocorr.plot(mcmc.trace.burned.thinned1_1)
plot(mcmc.trace.burned.thinned1_1)

plot(mcmc.trace1_2)
plot(mcmc.trace.burned1_2)
autocorr.plot(mcmc.trace.burned1_2)
autocorr.plot(mcmc.trace.burned.thinned1_2)
plot(mcmc.trace.burned.thinned1_2)
