# SET UP CCHF - starting April 2008
# adapted by J vesga , using here path package 
# 5 ages groups
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4434116/ Table 2
# Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4117442/ 
# "To avoid colostral immunity, animals were selected as follows: cattle were between 10 months and one-year of age, 
# and small ruminants were between three to eight months of age. "
# source("/home/metras/CCHF/setup_CCHF.R")

## 1. Human data incidence
###########################################################################################################
if(country=="AFG"){
  temp_start_date <-as.Date("2008-04-01")
  end_date<- as.Date("2008-10-31")
  temp_start_mo <-strftime(temp_start_date, format = "%Y-%m")
  end_mo<- strftime(end_date, format = "%Y-%m")
  fwork_start<-as.Date("2009-08-30")
  fwork_end <-as.Date("2009-08-15")
  fwork_mid <-as.Date("2009-08-31")
  
  fwork_starty<-as.Date("2008-09-01")
  fwork_endy <-as.Date("2009-07-31")
  fwork_midy <-as.Date("2009-02-28")
  format_date<-"%d/%m/%Y"
  
  
  
  # Get weekly data to montly 
  human_inc1 <- read.csv(here("data",country,"cchf_human_weekly2008.csv"), header=TRUE)#, sep=,)
  human_inc1$date<- as.Date(human_inc1$date, format="%d/%m/%Y")
  human_inc1$month<-strftime(human_inc1$date, format = "%Y-%m")
  human_inc1<-human_inc1%>%
    group_by(month)%>%
    summarise(cases=sum(obs))
  human_inc1<-human_inc1[human_inc1$month>=temp_start_mo,]
  human_inc1$time<-seq(1,NROW(human_inc1$cases),1)
  
  #Get monthly data
  human_inc2 <- read.csv(here("data",country,"cchf_human_monthly2017_2018.csv"), header=TRUE)
  human_inc2$date<- as.Date(human_inc2$date, format="%d/%m/%Y")
  human_inc2$month<-strftime(human_inc2$date, format = "%Y-%m")
  human_inc2$time<-1+(interval(temp_start_date,human_inc2$date) %/% months(1)) 
  
  human_inc2 <- human_inc2%>%
    rename(cases=obs)%>%
    select(month,cases,time)
  # Final monthly cases data
  human_month<-rbind(human_inc1,human_inc2)
  human_month<-human_month%>%
    mutate(cases=ifelse(cases==0,NA,as.numeric(paste(cases))))
  
  # Get Yearly data 
  human_inc3 <- read.csv(here("data",country,"cchf_human_yearly2009_2016.csv"), header=TRUE)
  human_inc3$date<- as.Date(human_inc3$date, format="%d/%m/%Y")
  human_inc3$year<-strftime(human_inc3$date, format = "%Y")
  human_year<-human_inc3%>%
    rename(cases=obs)
  
  human_year<-human_year%>%
    mutate(cases=ifelse(cases==0,NA,as.numeric(paste(cases))))
  human_year<-human_year[human_year$year>=2015,]
  
  years_fit<-human_year$year
  
}else if(country=="SA"){ 
  temp_start_date <-as.Date("2016-01-01")
  end_date<- as.Date("2017-12-31")
  fwork_start<-as.Date("2017-05-01")
  fwork_end <-as.Date("2017-11-30")
  fwork_mid <-as.Date("2017-07-15")
  
  fwork_starty<-as.Date("2017-05-01")
  fwork_endy <-as.Date("2018-05-31")
  fwork_midy <-as.Date("2017-11-15")
  format_date<-"%Y-%m-%d"
  
}

year_vec<-seq(2008,2019,by=1)

# define the dates of the field work August 2009
field_work_start<-(interval(temp_start_date,fwork_start) %/% months(1))
field_work_end <- (interval(temp_start_date,fwork_end) %/% months(1)) 
field_work_mid <- (interval(temp_start_date,fwork_mid) %/% months(1)) 

# year_work_start<- human_inc[human_inc$dates==fwork_starty,"time"] 
# year_work_end<- human_inc[human_inc$dates==fwork_endy,"time"] 
# year_work_mid<- human_inc[human_inc$dates==fwork_midy,"time"] 

# 2. Import human prevalence by occupation farming and other occupations (from CCHF 2009 dataset)
#########################################################################################################################
prev_other <- read.csv(here("data",country,'others.csv'), header=TRUE)
prev_other$low_ci <- round(binconf(prev_other$Posa, prev_other$denoma)[2], digit=4)
prev_other$up_ci <- round(binconf(prev_other$Posa, prev_other$denoma)[3], digit=4)

prev_farmer <- read.csv(here("data",country,'farmers.csv'), header=TRUE)
prev_farmer$low_ci <- round(binconf(prev_farmer$Posa, prev_farmer$denoma)[2], digit=4)
prev_farmer$up_ci <- round(binconf(prev_farmer$Posa, prev_farmer$denoma)[3], digit=4)

# 3. Livestock age-stratified prevalence from 2009 (5 age-groups)
#########################################################################################################################
prev_liv<-read.csv(here("data",country,'Live_Sero.csv'), header=TRUE)
prev_liv$low_ci <- round(binconf(prev_liv$pos_igg, prev_liv$denom)[,2], digit=4)
prev_liv$up_ci <- round(binconf(prev_liv$pos_igg, prev_liv$denom)[,3], digit=4)

prev_liv_all <- read.csv(here("data",country,'Live_Sero_wo_age.csv'), header=TRUE)
prev_liv_all$low_ci <- round(binconf(prev_liv_all$pos_igg, prev_liv_all$denom)[,2], digit=4)
prev_liv_all$up_ci <- round(binconf(prev_liv_all$pos_igg, prev_liv_all$denom)[,3], digit=4)

# 4. Import Soil Temperature April 2005 - March 2014 (informs the probability of tick transmission)
#########################################################################################################################
# csv file for Ingil district, Afghanistan - Median- Source ERA5 corpernicus soil temperature level 1 (0-7cm)

temp_day<- read.csv(here("data",country,'soil_temp_extended.csv'), header=TRUE, sep=',')
temp_day$date <- as.Date(temp_day$vector_days,  format=format_date)
temp_day <- temp_day[temp_day$date >= temp_start_date,]
temp_day$month <- strftime(temp_day$date, format = "%Y-%m")

temp_month<-temp_day%>%
  group_by(month)%>%
  summarise(soil_t=mean(Temperature))
temp_month$time <- c(1:nrow(temp_month))

soil_t<- temp_month$soil_t

# # transmission seasons starts (when outbreaks start April 1st) and ends on 31/10
# Transm_Start <- temp_day[temp_day$date==temp_start_date, "time"]
# Transm_End <- temp_day[temp_day$date==end_date, "time"] # 2004-10-31
# 
# # Date of spring start 2008 (1st April 2008)
# spring_2008_start_date <- temp_day[temp_day$date==temp_start_date, "time"] # 1097
# Spring08 <- spring_2008_start_date



# Index for monthly human cases
tmp<-which(!is.na(human_month$cases))




observations<-list(
  
  start_date=temp_start_date,
  start_month=temp_start_mo,
  index_mo_cases=human_month$time[tmp],
  cases_human_mo = human_month$cases[tmp],
  
  index_yr_cases=which(year_vec%in%years_fit),
  cases_human_yr = human_year$cases,
  
  prev_liv_age_pos_IgG=prev_liv$pos_igg,
  prev_liv_age_denom=prev_liv$denom,
  prev_liv_age=prev_liv$prop,
  prev_liv_age_low_ci=prev_liv$low_ci,
  prev_liv_age_up_ci=prev_liv$up_ci,
  
  prev_liv_all_pos_IgG=prev_liv_all$pos_igg,
  prev_liv_all_denom=prev_liv_all$denom,
  prev_liv_all=prev_liv_all$prop,
  prev_liv_all_low_ci=prev_liv_all$low_ci,
  prev_liv_all_up_ci=prev_liv_all$up_ci,
  
  prev_farmer_Posa=prev_farmer$Posa,
  prev_farmer_denoma=prev_farmer$denoma,
  prev_farmer=prev_farmer$prev,
  prev_farmer_low_ci=prev_farmer$low_ci,
  prev_farmer_up_ci=prev_farmer$up_ci,
  
  prev_other_Posa=prev_other$Posa,
  prev_other_denoma=prev_other$denoma,
  prev_other=prev_other$prev,
  prev_other_low_ci=prev_other$low_ci,
  prev_other_up_ci=prev_other$up_ci
  
  
)

tmp<-data.frame(
  time=c(observations$index_mo_cases,observations$index_yr_cases*12),
  obs=c(observations$cases_human_mo,observations$cases_human_yr)
)
data1_human <- tmp[order(tmp$time),]

data2_human<-data.frame(
  time=c(field_work_start:field_work_end),
  num_f=observations$prev_farmer_Posa,
  deno_f=observations$prev_farmer_denoma,
  num_o=observations$prev_other_Posa,
  deno_o=observations$prev_other_denoma)




#Functions

# Find FOI factor according to Soil temperature
temp_cap<-30
temp_foi_func<-function(celsius,temp_foi_factor,min_limit=min(soil_t)){
  if (celsius<=temp_cap){
    x<-(temp_foi_factor*(celsius - min_limit))

  }else if(celsius>temp_cap){
    x<-temp_foi_factor*(temp_cap-(celsius-temp_cap)-min_limit)
  }
  return(x)
}
# 
# temp_foi_func<-function(celsius,temp_foi_factor,min_limit=min(soil_t)){
#   x<-(temp_foi_factor*(min(celsius,30) - min_limit))
#   return(x)
# }
