plot_fits_function<-function(sim,observations, to_save=0){
  
  
  
  # if (nrow(sim$h_inc_month)>1)
  #   # #################################################################
  # # PLot mutiple runs (uncertanty and median)
  # # #################################################################
  # 
  # {
  #   
  ############### Livestock age prevalence vs survey Data
  x_d<-c(1, 3, 5, 7, 9)
  
  df_d<-data.frame(x=x_d,
                   prev=observations$prev_liv_age,
                   low=observations$prev_liv_age_low_ci,
                   up=observations$prev_liv_age_up_ci)
  
  df_s<-data.frame(x=x_d,prev=as.numeric(sim$l_prev_age))
  
  l_prev_qtls <-as.data.frame(rowQuantiles(t(sim$l_prev_age),
                                           probs=c(0.025,0.5,0.975)))
  
  l_prev_qtls_peak <-as.data.frame(rowQuantiles(t(sim$l_prev_age_peak),
                                                probs=c(0.025,0.5,0.975)))
  l_prev_qtls_low <-as.data.frame(rowQuantiles(t(sim$l_prev_age_low),
                                               probs=c(0.025,0.5,0.975)))
  
  l_prev_qtls$x<-x_d
  l_prev_qtls_peak$x<-x_d
  l_prev_qtls_low$x<-x_d
  
  p1 <- ggplot(data=df_d,aes(x=x))+ 
    geom_point(aes(y=prev))+
    geom_errorbar(aes(ymin=low, ymax=up), 
                  width=.2,position=position_dodge(.9)) + theme_bw()+
    geom_ribbon(data=l_prev_qtls,aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", 
                alpha=0.2)+
    geom_line(data=l_prev_qtls, aes(x = x, y=`50%`), col="blue")+
    
    geom_ribbon(data=l_prev_qtls_peak,aes(ymin=`2.5%`,ymax=`97.5%`), fill="green", 
                alpha=0.2)+
    geom_line(data=l_prev_qtls_peak, aes(x = x, y=`50%`), col="green")+
    
    geom_ribbon(data=l_prev_qtls_low,aes(ymin=`2.5%`,ymax=`97.5%`), fill="red", 
                alpha=0.2)+
    geom_line(data=l_prev_qtls_low, aes(x = x, y=`50%`), col="red")+
    
    ylab("Prevalence") +
    scale_x_discrete(name ="Livestock age (yrs)", 
                     limits=c("0-1","","1-2","","2-3","","3-4","","4+",""))+
    ylim(c(0,1))
  
  
  ################### Prevalence by age
  
  df0<-as.data.frame(rowQuantiles(t(sim$l_prev_0_long),probs=c(0.025,0.5,0.975)))
  df0$x<-seq(1,params$nt)
  df1<-as.data.frame(rowQuantiles(t(sim$l_prev_1_long),probs=c(0.025,0.5,0.975)))
  df1$x<-seq(1,params$nt)
  df2<-as.data.frame(rowQuantiles(t(sim$l_prev_2_long),probs=c(0.025,0.5,0.975)))
  df2$x<-seq(1,params$nt)
  df3<-as.data.frame(rowQuantiles(t(sim$l_prev_3_long),probs=c(0.025,0.5,0.975)))
  df3$x<-seq(1,params$nt)
  df4<-as.data.frame(rowQuantiles(t(sim$l_prev_4_long),probs=c(0.025,0.5,0.975)))
  df4$x<-seq(1,params$nt)
  
  
  
  p2<-ggplot(df0, aes(x=x))+
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", alpha=0.2)+
    geom_line(aes(y=`50%`, color="Age0_1"))+
    geom_ribbon(data=df1,aes(ymin=`2.5%`,ymax=`97.5%`), fill="red", alpha=0.2)+
    geom_line(data=df1,aes(y=`50%`,color="Age1_2"))+
    geom_ribbon(data=df2,aes(ymin=`2.5%`,ymax=`97.5%`), fill="orange", alpha=0.2)+
    geom_line(data=df2,aes(y=`50%`,color="Age2_3"))+
    geom_ribbon(data=df3,aes(ymin=`2.5%`,ymax=`97.5%`), fill="purple", alpha=0.2)+
    geom_line(data=df3,aes(y=`50%`, color="Age3_4"))+
    geom_ribbon(data=df4,aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
    geom_line(data=df4,aes(y=`50%` ,color="Age4+"))+
    scale_colour_manual("", 
                        breaks = c("Age0_1", "Age1_2", "Age2_3","Age3_4","Age4+"),
                        values = c("blue", "red", "orange","purple","darkgreen")) +
    ylab("Seroprevalence") + xlab("Month")+ylim(0,1)
  
  
  
  # R effective in livestock 
  
  r_temp_livestock<-matrix(NA,params$nt,nrow(theta))
  
  for (jj in 1:nrow(theta) ){
    
    for (tt in 1:length(soil_t)){
       r_temp_livestock[tt,jj]<-envrdriver_foi_func(environment_factor[tt],theta$A[jj])*sim$sus_l_frac[jj,tt]
      
    }
  }
  
  
  Rtl <-as.data.frame(rowQuantiles((r_temp_livestock),probs=c(0.025,0.5,0.975)))
  Rtl$x<-seq(1,params$nt)
  Rtl$group<-"livestock"
  
  p3 <- ggplot(data=Rtl,aes(x=x))+ 
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`, color=group), alpha=0.3)+
    geom_line(aes(y=`50%`, color=group))+
    ylab("R_t") + xlab("day")+xlim(1,130) +ylim(0,10)
  
  
 
  ################## Human reported cases 
  mo<-seq( 1, by=1, len=ncol(sim$h_prev_farmer_long))
  df_s <-as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
                                    probs=c(0.025,0.5,0.975)))
  df_s$x<-mo
  df_d<-data.frame(x=observations$index_mo_cases,cases=observations$cases_human_mo)
  
  df_sim<-data.frame(x=mo, t(sim$h_incidence_reported_month))
  dat_sim<-reshape2::melt(df_sim, id="x")
  
  p4<- ggplot(data=df_s, aes(x=x))+
    geom_line(data=dat_sim, aes(x=x, y=value, group=variable), col="grey", 
              alpha=0.2, lwd=0.4) +
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
    geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
    geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
    # geom_vline(xintercept=which.max(df_s$`50%`), color="red") +
    xlab("Month")+ylab("Monthly cases Reported in humans")
  
  
  ####################### Human reported cases yearly 
  yr<-seq( 1, by=1, len=length(unique(as.Date(temp_month$month,format="%Y"))))+2007
  df_s <-as.data.frame(rowQuantiles(t((sim$h_incidence_reported_year)),
                                    probs=c(0.025,0.5,0.975)))
  df_s$x<-yr
  df_d<-data.frame(x=observations$index_yr_cases+2007,cases=observations$cases_human_yr)
  df_sim<-data.frame(x=yr, t(sim$h_incidence_reported_year))
  dat_sim<-reshape2::melt(df_sim, id="x")
  
  p5<- ggplot(data=df_s, aes(x=x))+
    geom_line(data=dat_sim, aes(x=x, y=value, group=variable), col="grey", 
              alpha=0.2, lwd=0.4) +
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
    geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
    geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
    xlab("Year")+ylab("Cases Reported in humans")
  
  ##=========== Prevalence in farmers and others
  df_sf <- as.data.frame(rowQuantiles(t(sim$h_prev_farmer_long),probs=c(0.025,0.5,0.975)))
  df_sf$x<- seq(1:params$nt)
  df_so <- as.data.frame(rowQuantiles(t(sim$h_prev_other_long),probs=c(0.025,0.5,0.975)))
  df_so$x<- seq(1:params$nt)
  
  
  p6<- ggplot(data=df_sf,aes(x=x))+
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
    geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
    geom_ribbon(data=df_so,aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", alpha=0.2)+
    geom_line(data=df_so, aes(y=`50%`), col="blue", lwd=0.4) +
    geom_point(aes(x=field_work_mid,y=prev_farmer$prev))+
    geom_errorbar(aes(x=field_work_mid,ymin=prev_farmer$low_ci, ymax=prev_farmer$up_ci), width=.2,
                  position=position_dodge(.9)) + 
    geom_point(aes(x=field_work_mid,y=prev_other$prev))+
    geom_errorbar(aes(x=field_work_mid,ymin=prev_other$low_ci, ymax=prev_other$up_ci), width=.2,
                  position=position_dodge(.9))+ theme_bw()+
    ylab("IgG seroprevalence (humans)") + xlab("Month")
  
  
  
  
  ####################### Human fatalities yearly 
  yr<-seq( 1, by=1, len=length(unique(as.Date(temp_month$month,format="%Y"))))+2007
  df_s <-as.data.frame(rowQuantiles(t((sim$h_fatality_year)),
                                    probs=c(0.025,0.5,0.975)))
  df_s$x<-yr
  df_d<-data.frame(x=observations$index_yr_deaths+2007,cases=observations$deaths_human_yr)
  df_sim<-data.frame(x=yr, t(sim$h_fatality_year))
  dat_sim<-reshape2::melt(df_sim, id="x")
  
  p7<- ggplot(data=df_s, aes(x=x))+
    geom_line(data=dat_sim, aes(x=x, y=value, group=variable), col="grey", 
              alpha=0.2, lwd=0.4) +
    geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
    geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
    geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
    xlab("Year")+ylab("CCHF fatal cases in humans")+
    xlim(2008, 2018)
  
  
  ################## Human spectrum of incidence 
  

  df_d<-data.frame(x=observations$index_mo_cases,cases=observations$cases_human_mo)
  mo<-seq( 1, by=1, len=ncol(sim$h_prev_farmer_long))
  #reported
  df_rep <-as.data.frame(rowQuantiles(t((sim$h_incidence_reported_month)),
                                    probs=c(0.025,0.5,0.975)))
  df_rep$x<-mo
  # Symptomatic
  df_clin <-as.data.frame(rowQuantiles(t((sim$h_incidence_clinical_month)),
                                      probs=c(0.025,0.5,0.975)))
  df_clin$x<-mo
  # All incidence
  df_all <-as.data.frame(rowQuantiles(t((sim$h_incidence_month)),
                                       probs=c(0.025,0.5,0.975)))
  df_all$x<-mo

  longx<-seq( 1,ncol(sim$h_prev_farmer_long), by=1/30)
  
  df_all2<-approxm(df_all,length(longx),"linear")
  df_rep2<-approxm(df_rep,length(longx),"linear")
  df_clin2<-approxm(df_clin,length(longx),"linear")
  
  
  
  
  p8<- ggplot(data=df_all2, aes(x=x))+
    geom_ribbon(aes(ymin=`2.5%`*0,ymax=`97.5%`, fill="All"), alpha=1)+
    geom_ribbon(data=df_clin2,aes(ymin=`2.5%`*0,ymax=`97.5%`, fill="Clinical"), alpha=1)+
    geom_ribbon(data=df_rep2,aes(ymin=`2.5%`*0,ymax=`97.5%`, fill="Reported"), alpha=1)+
    geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
    xlim(0,130)+
    scale_fill_manual("Incidence",values=c("All"="dodgerblue4", "Clinical"="grey52", "Reported"="orange3"))+
    xlab("Month")+ylab("Incidence")
  
  
  
  
  if(to_save==0){
    windows()
  }
  
  gridExtra::grid.arrange(p1,p2,p3)
  
  if(to_save==0){
    windows()
  }
  
  gridExtra::grid.arrange(p4,p5,p6,p7)
  
  if(to_save==0){
    windows()
  }
  
  gridExtra::grid.arrange(p8)
  
  
}


