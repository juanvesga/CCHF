#Merge chains  CCHF May 2021

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

country<- "SA"

# analyse MCMC samples


# Import chains
trace1_1 <- read.csv(here("output",country,"trace_chain1.csv"))
trace1_2 <- read.csv(here("output",country,"trace_chain2.csv"))

######################################################################
# 1. chain 1
head(trace1_1)
mcmc.trace1_1 <- mcmc(trace1_1)
summary(mcmc.trace1_1)        

acceptanceRate <- 1 - rejectionRate(mcmc.trace1_1)
acceptanceRate
effectiveSize(mcmc.trace1_1)
plot(mcmc.trace1_1)
mcmc.trace.burned1_1 <- burnAndThin(mcmc.trace1_1, burn = round(nrow(trace1_1)*0.5))
plot(mcmc.trace.burned1_1)
autocorr.plot(mcmc.trace.burned1_1)
mcmc.trace.burned.thinned1_1<- burnAndThin(mcmc.trace.burned1_1, thin = 50) 
autocorr.plot(mcmc.trace.burned.thinned1_1)
plot(mcmc.trace.burned.thinned1_1)
summary(mcmc.trace.burned.thinned1_1)

#################
# 2. chain 2
head(trace1_2)
mcmc.trace1_2 <- mcmc(trace1_2)
summary(mcmc.trace1_2)        

acceptanceRate <- 1 - rejectionRate(mcmc.trace1_2)
acceptanceRate
effectiveSize(mcmc.trace1_2)
plot(mcmc.trace1_2)
mcmc.trace.burned1_2 <- burnAndThin(mcmc.trace1_2, burn = round(nrow(trace1_2)*0.25))
plot(mcmc.trace.burned1_2)
autocorr.plot(mcmc.trace.burned1_2)
mcmc.trace.burned.thinned1_2<- burnAndThin(mcmc.trace.burned1_2, thin = 20) 
autocorr.plot(mcmc.trace.burned.thinned1_2)
plot(mcmc.trace.burned.thinned1_2)
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

write.table(df1_1, here("output",country,"posteriors.txt"), sep = "\t" )
#####################################################################################


