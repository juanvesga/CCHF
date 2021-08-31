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
  fwork_start<-as.Date("2009-08-30")
  fwork_end <-as.Date("2009-08-15")
  fwork_mid <-as.Date("2009-08-31")
  
  fwork_starty<-as.Date("2008-09-01")
  fwork_endy <-as.Date("2009-07-31")
  fwork_midy <-as.Date("2009-02-28")
  format_date<-"%d/%m/%Y"
  
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



human_inc2 <- read.csv(here("data",country,"cchf_human05_14.csv"))#, sep=,)
head(human_inc2, 15)

human_inc2$time <- c(1:nrow(human_inc2))
human_inc2$dates <- as.Date(human_inc2$date, format="%d/%m/%Y")
human_inc2$week <- strftime(human_inc2$dates, format = "%Y-W%V")
head(human_inc2)
human_inc <- human_inc2[human_inc2$dates>=temp_start_date,]
human_inc$time <- c(1:nrow(human_inc))

# define the dates of the field work August 2009
field_work_start<- human_inc[human_inc$dates==fwork_start,"time"] 
field_work_end <- human_inc[human_inc$dates==fwork_end,"time"] 
field_work_mid <- human_inc[human_inc$dates==fwork_mid,"time"] 

year_work_start<- human_inc[human_inc$dates==fwork_starty,"time"] 
year_work_end<- human_inc[human_inc$dates==fwork_endy,"time"] 
year_work_mid<- human_inc[human_inc$dates==fwork_midy,"time"] 

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


temp_day2<- read.csv(here("data",country,'soil_temp.csv'), header=TRUE, sep=',')
temp_day2$date <- as.Date(temp_day2$vector_days,  format=format_date)
temp_day <- temp_day2[temp_day2$date >= temp_start_date,]
temp_day$time <- c(1:nrow(temp_day))

head(temp_day)

# transmission seasons starts (when outbreaks start April 1st) and ends on 31/10
Transm_Start <- temp_day[temp_day$date==temp_start_date, "time"]
Transm_End <- temp_day[temp_day$date==end_date, "time"] # 2004-10-31

# Date of spring start 2008 (1st April 2008)
spring_2008_start_date <- temp_day[temp_day$date==temp_start_date, "time"] # 1097
Spring08 <- spring_2008_start_date

soil_t<- temp_day$Temperature


# Data object for calibartion and plotting
tmp<-data.frame(date = seq( temp_start_date, by=1, len=length(human_inc$obs)),
                   time = c(1:(length(human_inc$obs))),
                   cases=human_inc$obs,
                   week=rep(seq(length(human_inc$obs)), each=7, length.out=length(human_inc$obs)))

tmp2<-tmp%>%
  group_by(week)%>%
  summarise(reported=sum(cases),.groups = 'drop')%>%
  mutate(reported=ifelse(reported==0,NA,as.numeric(paste(reported))))

ind_h<-which(!is.na(tmp2$reported))



observations<-list(
  h_inc_week.x = ind_h,
  h_inc_week = tmp2$reported[ind_h],
  
  
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




#Functions

# Find FOI factor according to Soil temperature
temp_foi_func<-function(x,temp_foi_factor,min_limit=min(soil_t)){
  x<-(temp_foi_factor*(min(x,30) - min_limit))
  return(x)
}
