###########################################################################################################
# Part 1:  BUILD MODEL STRUCTURE
###########################################################################################################
source(here("src","get_addresses.R"))


# Model states
gps<-list( age=c('a1','a2','a3','a4','a5'),
           occ=c('farmer','other'))

states0 <- c('L_Ri')
states1 <- c('L_S','L_I','L_R')
states2 <- c('S','E','I','R')


# use addresses function to map states and indexes
i<-list()
s<-list()
groups<-list(states0)
ref<- get_addresses(groups, i, s, 0)
groups<-list(states1,gps$age)
ref<- get_addresses(groups, ref$i, ref$s, ref$i$nstates)
groups<-list(states2,gps$occ)
ref<-  get_addresses(groups, ref$i, ref$s, ref$i$nstates)

# Include the auxiliaries
auxnames <- c('inc','mort')
auxinds  <- c(  3  , 1    )
ref$i$aux<-list()
lim_all<-ref$i$nstates
for (ii in 1:length(auxnames)){
  inds <- lim_all + (1:auxinds[ii])
  ref$i$aux[[auxnames[ii]]] <- inds
  lim_all <- inds[length(inds)]
}
ref$i$nx <- lim_all


#States
ref$s$infectious_l <- c(ref$s$L_I)
ref$s$livestock   <- c(ref$s$L_Ri, ref$s$L_S,ref$s$L_I,ref$s$L_R)
ref$s$humans      <- c(ref$s$farmer, ref$s$other)
ref$s$prevalent_l <- c(ref$s$L_R)
ref$s$prevalent_f <- c(ref$i$R$farmer)
ref$s$prevalent_o <- c(ref$i$R$other)

ref$s$nstates<-(1:ref$i$nstates)

#__________________________________________________________________________
#  Prepare selectors and aggregartors of model output
#__________________________________________________________________________

#forbidden transitions
tmp <- matrix(1,ref$i$nstates,ref$i$nstates)
tmp[ref$s$livestock,ref$s$humans]<-0
tmp[ref$s$humans,ref$s$livestock]<-0
tmp[ref$s$farmer,ref$s$other]<-0
tmp[ref$s$other,ref$s$farmer]<-0

check<-tmp - diag(diag(tmp))

#  Incidence: 1.livestock 2.farmer 3. other
tmp <- matrix(0,3,ref$i$nstates)
tmp[1,ref$s$L_I] <- 1
tmp[2,ref$i$I$farmer] <- 1
tmp[3,ref$i$I$other] <- 1
agg <-list(inc= tmp)

tmp <-matrix(0,ref$i$nstates,ref$i$nstates)
tmp[ref$s$L_I,] <- 1
tmp[ref$s$I,] <- 1
tmp<-tmp*check
sel<-list(inc= tmp - diag(diag(tmp)))

# Selectors mortality
tmp <- matrix(0,1,ref$i$nstates)
tmp[1,ref$s$livestock] <- 1
agg$mort <-tmp

ref$agg<- agg
ref$sel<- sel

###########################################################################################################
# Part 2:  MODEL PARAMETERS
###########################################################################################################

params<-list(
  ###########################################################################################################
  # Time steps in days starts 1st April 2005 ends 31st March 2014
  ###########################################################################################################
  nt =  nrow(temp_month),
  # 2.1. Livestock parameters
  ###########################################################################################################
  # Number of livestock age groups
  ya = 5,
  # Population size N animals based on FAO report from region in 2008
  N_liv = 15193,
  
  # livestock population age structure  and death rates
  # proportion of each age-group in the population
    pop_st = c(0.4918, 0.2499, 0.1270, 0.0646, 0.0667),
   # pop_st = c(0.4234375, 0.2890625, 0.1375, 0.11875, 0.03125), # Mauritania
   # pop_st = pop <- c(0.2743, 0.2654, 0.186, 0.0973, 0.177),
   # pop_st=c(0.25714286, 0.09028571, 0.23657143, 0.25142857, 0.16457143),
   mu_1=0.002643826*30, 
  
  #daily death rate
  # deathd = c(0.002643826,0.002643826,0.002643826,0.002643826,0.002643826),
   deathd=c(0.07618946, 0.07436794, 0.07466935, 0.07447458, 0.07470417) , 
  # deathd=c(1/(5*12), 1/(4*12), 1/(3*12), 1/(2*12), 1/(1*12)) , 
  # Proportion of immune livestock at t0. Source?
  #imm_t0 <- c(0.05, 0.08,0.1,0.12,0.15)*4
   imm_t0_bulgaria = c(0.29,0.48, 0.8, 0.87, 0.87), # bulgaria paper  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4434116/ Table 2
   # imm_t0 = c(0.3,0.3, 0.3, 0.3, 0.3), # bulgaria paper  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4434116/ Table 2
  imm_t0=c(0.29,0.48, 0.8, 0.87, 0.87)*0.5,
  #imm_t0<- prev_liv$imm_t0
  # Transovarial transmission for seeding infection every year on April 1st
  TOT = 1,
  
  # Duration of infectiousness in livestock (in days)
  D_inf_L = 7/30,
  
  D_lact_liv=6, #(months)
  
  D_imm_liv=(5*12),
  
  # 2.2. Parameters for humans
  ###########################################################################################################
  # total population of farmers checked from world bank and USAID report from 2008 - 30% pop farmers
  N_hum = 25382,
  prop_o = 0.7, # proportion of other occupations
  
  #7615
  
  # Proportion of farmers immune at time 0. from data
  prop_F_R= 0.095,#0.133333, # 0.0952
  # Proportion of other occupations immune at time 0. from data
  prop_O_R= 0.041,#0.046875, ###0.0411
  
  # Reporting fraction. fraction of human infections reported by the surveillance system. Source?
  RR = 0.10,
  RRreport=0.75, # Increse of reporting over the yeras (2013 to 2019)
  
  # Latent period in humans in days (to calculate the rate at which human move from latent to infectious)
  D_lat_H = 4/30,
  # Infectious period in humans in days (to calculate the rate at which human move from infectious to immune)
  D_inf_H = 9/30,
  # Duration immunity in humans
  D_imm_H = 10*12,
  # case fatality rate
  CFR = 0.33,
  # Clinical fraction
  clin_frac = 1-0.8, # Turkey data (Bodur et al.) 
  # mortality and birth rates in human for constant pop size
  L = 61.5, # life span in years  between 2008 and 2014 https://data.worldbank.org/indicator/SP.DYN.LE00.IN?locations=AF 
  b_d = (1/(61.5*12)), # daily birth and death rate for constant pop 
  
  sp_hz=10, # spline hazard multiplier. This value will be multiplied by spline 
  
  
  ## baseline vaccination paramts
  # Baseline conditions
  start_vax= 70, 
  stop_vax = 73,
  vacc_eff = 0.9,
  time_imm = 0.5,
  pop_coverage_f = 0,
  pop_coverage_o = 0,
  p_vacc_livestock= 0,
  liv_vacc_yearly=FALSE
  
  
)

# Rate at which livestock and humans becomes infectious and immune
params$time_immune_livestock <- 1/params$D_inf_L # rate at which infectious livestock become immune
params$time_to_infous <- 1/params$D_lat_H # rate at which latent human become infectious
params$time_immune_human <- 1/params$D_inf_H # rate at which infectious human become immune
params$time_susceptible_human <- 1/params$D_imm_H # rate at which humans become susceptible again
params$time_passimm_loss_livestock<- 1/params$D_lact_liv # rate of moving from passive immnity at birth to susceptible
params$time_susceptible_livestock<- 1/params$D_imm_liv # rate of moving from passive immnity at birth to susceptible
params$time_protection<- 1/params$time_imm
# Added parameters defined with others 
params$indices_infL<- c(6,7,8,9,10)-1
params$indices_livestock<-seq(1,15)-1


# total population of other occupations
params$N_O <- round(params$N_hum*params$prop_o, digit=0)
#17767
# total farming population
params$N_F <- params$N_hum - params$N_O
params$Birth_F <- params$N_F*params$b_d
params$Birth_O <- params$N_O*params$b_d
params$Birth <-sum(params$pop_st*params$deathd) 
params$Ageing <- (1-params$deathd)/12
params$recover <- params$time_immune_livestock# (1-params$deathd-params$Ageing)*params$time_immune_livestock

# pass references
params$ref<-ref

