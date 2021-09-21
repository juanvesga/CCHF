plot_fits_function<-function(sim,observations){
  
  # windows()
  
  if (nrow(sim$h_inc_month)>1)
    # #################################################################
  # PLot mutiple runs (uncertanty and median)
  # #################################################################
  
  {
    
    # Livestock age prevalence
    x_d<-c(1, 3, 5, 7, 9)
    
    df_d<-data.frame(x=x_d,
                     prev=observations$prev_liv_age,
                     low=observations$prev_liv_age_low_ci,
                     up=observations$prev_liv_age_up_ci)
    
    df_s<-data.frame(x=x_d,prev=as.numeric(sim$l_prev_age))
    
    l_prev_qtls <-as.data.frame(rowQuantiles(t(sim$l_prev_age),
                                             probs=c(0.025,0.5,0.975)))
    l_prev_qtls$x<-x_d
    
    p1 <- ggplot(data=df_d,aes(x=x))+ 
      geom_point(aes(y=prev))+
      geom_errorbar(aes(ymin=low, ymax=up), 
                    width=.2,position=position_dodge(.9)) + theme_bw()+
      geom_ribbon(data=l_prev_qtls,aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", 
                  alpha=0.2)+
      geom_line(data=l_prev_qtls, aes(x = x, y=`50%`), col="blue")+
      ylab("Prevalence") +
      scale_x_discrete(name ="Livestock age (yrs)", 
                       limits=c("0-1","","1-2","","2-3","","3-4","","4+",""))+
      ylim(c(0,1))
    
    
    
    # Human cases 
    mo<-seq( 2, by=1, len=ncol(sim$h_prev_farmer_long)-1)
    df_s <-as.data.frame(rowQuantiles(t((sim$h_inc_mo)),
                                      probs=c(0.025,0.5,0.975)))
    df_s$x<-mo
    df_d<-data.frame(x=observations$index_mo_cases,cases=observations$cases_human_mo)
    
    df_sim<-data.frame(x=mo, t(sim$h_inc_mo))
    dat_sim<-reshape2::melt(df_sim, id="x")
    
    p2<- ggplot(data=df_s, aes(x=x))+
      geom_line(data=dat_sim, aes(x=x, y=value, group=variable), col="grey", 
                alpha=0.2, lwd=0.4) +
      geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
      geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
      geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
      # geom_vline(xintercept=which.max(df_s$`50%`), color="red") +
      xlab("Month")+ylab("Monthly cases Reported in humans")
    
    # Human cases yearly 
    yr<-seq( 1, by=1, len=length(unique(as.Date(temp_month$month,format="%Y"))))+2007
    df_s <-as.data.frame(rowQuantiles(t((sim$h_inc_year)),
                                      probs=c(0.025,0.5,0.975)))
    df_s$x<-yr
    df_d<-data.frame(x=observations$index_yr_cases+2007,cases=observations$cases_human_yr)
    df_sim<-data.frame(x=yr, t(sim$h_inc_year))
    dat_sim<-reshape2::melt(df_sim, id="x")
    
    p3<- ggplot(data=df_s, aes(x=x))+
      geom_line(data=dat_sim, aes(x=x, y=value, group=variable), col="grey", 
                alpha=0.2, lwd=0.4) +
      geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
      geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
      geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
      xlab("Year")+ylab("Cases Reported in humans")+
      xlim(2008, 2018)
    
    ## Prevalence
    
    df_sf <- as.data.frame(rowQuantiles(t(sim$h_prev_farmer_long),probs=c(0.025,0.5,0.975)))
    df_sf$x<- seq(1:params$nt)
    df_so <- as.data.frame(rowQuantiles(t(sim$h_prev_other_long),probs=c(0.025,0.5,0.975)))
    df_so$x<- seq(1:params$nt)
    
    
    p4<- ggplot(data=df_sf,aes(x=x))+
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
    
    # # R_t dependent on temp
    # 
    # 
    # r_temp<-matrix(NA,params$nt,nrow(theta))
    # for (jj in 1:nrow(theta) ){
    #   for (tt in 1:length(soil_t)){
    #     r_temp[tt,jj]<-temp_foi_func(soil_t[tt],theta$A[jj])
    #   }
    # }
    # 
    # 
    # Rt <-as.data.frame(rowQuantiles((r_temp),probs=c(0.025,0.5,0.975)))
    # Rt$x<-seq(1,params$nt)
    # Rt$soil<-soil_t
    # 
    # 
    # 
    # p5<-    ggplot(data=Rt,aes(x=x))+ 
    #   geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", alpha=0.2)+
    #   geom_line(aes(y=`50%`), col="blue")+
    #   ylab("R_t") + xlab("Month")
    
    # ggplot(data=Rt,aes(x=soil))+ 
    #   geom_point(aes(y=`50%`), col="red")+
    #   ylab("R_t") + xlab("Celsius")
    # R_t dependent on temp
    
    
    r_temp<-matrix(NA,params$nt,nrow(theta))
    for (jj in 1:nrow(theta) ){
      
      x<-seq(0,1,length.out = params$nt)
      knots <- as.numeric(c(theta$knot1[jj]/params$nt,theta$knot2[jj]/params$nt ,1000/params$nt ))
      betas = c(0,theta$beta1[jj], 0.1, 0.1, 0.1, 0.1)
      sdata <- genSpline(x, knots, 2, betas)
      
      
      for (tt in 1:length(soil_t)){
        r_temp[tt,jj]<-temp_foi_func(soil_t[tt],theta$A[jj])*
          sdata$dt$y.spline[tt]*params$sp_hz*params$D_inf_L*sim$sus_l_frac[jj,tt]
      }
    }
    
    
    Rt <-as.data.frame(rowQuantiles((r_temp),probs=c(0.025,0.5,0.975)))
    Rt$x<-seq(1,params$nt)
    Rt$soil<-soil_t
    
    
    
    p5<-    ggplot(data=Rt,aes(x=x))+ 
      geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", alpha=0.2)+
      geom_line(aes(y=`50%`), col="blue")+
      ylab("R_t") + xlab("day")
    
    
    gridExtra::grid.arrange(p1,p2,p3,p4,p5)
    
  }
  # #################################################################
  else # Run a single run vs data
    # #################################################################
  {
    # Livestock age prevalence
    
    x_d<-c(1, 3, 5, 7, 9)
    x_s<- x_d+0.5
    
    df_d<-data.frame(x=x_d,
                     prev=observations$prev_liv_age,
                     low=observations$prev_liv_age_low_ci,
                     up=observations$prev_liv_age_up_ci)
    df_s<-data.frame(x=x_d,prev=as.numeric(sim$l_prev_age))
    
    p1 <- ggplot(data=df_d,aes(x=x))+ 
      geom_point(aes(y=prev))+
      geom_errorbar(aes(ymin=low, ymax=up), width=.2,position=position_dodge(.9)) + 
      geom_line(data=df_s, aes(x = x, y=prev), col="blue")+
      ylab("Prevalence") +
      scale_x_discrete(name ="Livestock age (yrs)", 
                       limits=c("0-1","","1-2","","2-3","","3-4","","4+",""))+
      ylim(c(0,1))
    
    
    
    
    # Human cases 
    mo<-seq( 2, by=1, len=length(sim$h_prev_farmer_long)-1)
    df_s<-data.frame(x=mo,cases=as.numeric(sim$h_inc_mo))
    df_d<-data.frame(x=observations$index_mo_cases,cases=observations$cases_human_mo)
    
    
    p2<- ggplot(data=df_s, aes(x=x))+
      geom_line(aes(y=cases), col="darkgreen", lwd=0.4) +
      geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
      geom_vline(xintercept=which.max(df_s$cases), color="red") +
      ylab("Reported cases in humans")+xlab("Month")
    
    
    # Human cases yearly 
    yr<-seq( 1, by=1, len=length(unique(as.Date(temp_month$month,format="%Y"))))+2007
    df_s<-data.frame(x=yr,cases=as.numeric(sim$h_inc_year))
    df_d<-data.frame(x=observations$index_yr_cases+2007,cases=observations$cases_human_yr)
    
    
    p3<- ggplot(data=df_s, aes(x=x))+
      geom_line(aes(y=cases), col="darkgreen", lwd=0.4) +
      geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
      ylab("Reported cases in humans")+xlab("Year")+
      xlim(2008, 2018)
    
    
    ## Prevalence farmers and others
    
    df_sf<- data.frame(x=seq(1:params$nt), prev=as.numeric(sim$h_prev_farmer_long))
    df_so<- data.frame(x=seq(1:params$nt), prev=as.numeric(sim$h_prev_other_long))
    
    p4<- ggplot(data=df_sf,aes(x=x))+
      geom_line(aes(y=prev, colour = "Farmers"), lwd=0.4) +
      geom_line(data=df_so, aes(y=prev,colour = "Others"), lwd=0.4) +
      scale_colour_manual("", 
                          breaks = c("Farmers", "Others"),
                          values = c("darkgreen","blue"))+
      geom_point(aes(x=field_work_mid,y=prev_farmer$prev))+
      geom_errorbar(aes(x=field_work_mid,ymin=prev_farmer$low_ci, ymax=prev_farmer$up_ci), width=.2,
                    position=position_dodge(.9), col="darkgreen") + 
      geom_point(aes(x=field_work_mid,y=prev_other$prev))+
      geom_errorbar(aes(x=field_work_mid,ymin=prev_other$low_ci, ymax=prev_other$up_ci), width=.2,
                    position=position_dodge(.9),col="blue")+ theme_bw()+
      ylab("IgG seroprevalence (humans)") + xlab("Month")
    
    
    # R_t dependent on temp
    
    
    # r_temp<-matrix(NA,params$nt,nrow(theta))
    # for (jj in 1:nrow(theta) ){
    #   for (tt in 1:length(soil_t)){
    #     r_temp[tt,jj]<-temp_foi_func(soil_t[tt],theta$A[jj])
    #   }
    # }
    # 
    # 
    # Rt <-as.data.frame(rowQuantiles((r_temp),probs=c(0.025,0.5,0.975)))
    # Rt$x<-seq(1,params$nt)
    # Rt$soil<-soil_t
    # 
    # 
    # 
    # p5<-    ggplot(data=Rt,aes(x=x))+ 
    #   geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", alpha=0.2)+
    #   geom_line(aes(y=`50%`), col="blue")+
    #   ylab("R_t") + xlab("Month")
    # 
    # ggplot(data=Rt,aes(x=soil))+ 
    #   geom_point(aes(y=`50%`), col="red")+
    #   ylab("R_t") + xlab("Celsius")
    
    # R_t dependent on temp
    
    
    r_temp<-matrix(NA,params$nt,nrow(theta))
    for (jj in 1:nrow(theta) ){
      
      x<-seq(0,1,length.out = params$nt)
      knots <- as.numeric(c(theta$knot1[jj]/params$nt,theta$knot2[jj]/params$nt ,1000/params$nt ))
      betas = c(0,theta$beta1[jj], 0.1, 0.1, 0.1, 0.1)
      sdata <- genSpline(x, knots, 2, betas)
      
      for (tt in 1:length(soil_t)){
        r_temp[tt,jj]<-temp_foi_func(soil_t[tt],theta$A[jj])*
          sdata$dt$y.spline[tt]*params$sp_hz*params$D_inf_L*sim$sus_l_frac
      }
    }
    
    
    Rt <-as.data.frame(rowQuantiles((r_temp),probs=c(0.025,0.5,0.975)))
    Rt$x<-seq(1,params$nt)
    Rt$soil<-soil_t
    
    
    
    p5<-    ggplot(data=Rt,aes(x=x))+ 
      geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="blue", alpha=0.2)+
      geom_line(aes(y=`50%`), col="blue")+
      ylab("R_t") + xlab("day")
    
    gridExtra::grid.arrange(p1,p2,p3,p4,p5)
    
    
    
    
  }
  
}


