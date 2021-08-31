plot_fits_function<-function(sim,observations){
  
  if (nrow(sim$h_inc_week)>1)
    # #################################################################
  # PLot mutiple runs (uncertanty and median)
  # #################################################################
  
  {
    
    
    # Livestock age prevalence
    
    x_d<-c(1, 3, 5, 7, 9)
    x_s<- x_d+0.5
    
    df_d<-data.frame(x=x_d,
                     prev=observations$prev_liv_age,
                     low=observations$prev_liv_age_low_ci,
                     up=observations$prev_liv_age_up_ci)
    
    df_s<-data.frame(x=x_s,prev=as.numeric(sim$l_prev_age))
    
    
    l_prev_qtls <-as.data.frame(rowQuantiles(t(sim$l_prev_age),probs=c(0.025,0.5,0.975)))
    l_prev_qtls$x<-x_s
    
    p1 <- ggplot(data=df_d,aes(x=x))+ 
      geom_point(aes(y=prev))+
      geom_errorbar(aes(ymin=low, ymax=up), width=.2,position=position_dodge(.9)) + theme_bw()+
      geom_segment(data=df_s, aes(x = x, y=prev , xend = x+0.5, yend = prev), col="grey",alpha=0.4)+
      geom_point(data=l_prev_qtls, aes(x = x+0.25, y=`50%`), col="red")+
      geom_errorbar(data=l_prev_qtls,aes(x=x+0.25,ymin=`2.5%`, ymax=`97.5%`), width=.2,
                    position=position_dodge(.9), col="red")+
      ylab("Prevalence") +
      scale_x_discrete(name ="Livestock age (yrs)", limits=c("1","","2","","3","","4","","5",""))+
      ylim(c(0,1))
    
    
    
    # Human cases 

    
    df_s <-as.data.frame(rowQuantiles(t(sim$h_inc_week),probs=c(0.025,0.5,0.975)))
    df_s$x<-seq(1:(params$nt/7))
    df_d<-data.frame(x=observations$h_inc_week.x,cases=observations$h_inc_week)
    
    
    p2<- ggplot(data=df_s, aes(x=x))+
      geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`), fill="darkgreen", alpha=0.2)+
      geom_line(aes(y=`50%`), col="darkgreen", lwd=0.4) +
      geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
      geom_vline(xintercept=which.max(df_s$`50%`), color="red") +
      geom_text(x=150, y=8, label="peak (week 74)", size=3, col="red") +
      geom_text(x=160, y=7, label="25/08/2009-31/08/2009", size=3, col="red") 
    
    
    ## Prevalence
    
    df_sf <- as.data.frame(rowQuantiles(t(sim$h_prev_farmer_long),probs=c(0.025,0.5,0.975)))
    df_sf$x<- seq(1:params$nt)
    df_so <- as.data.frame(rowQuantiles(t(sim$h_prev_other_long),probs=c(0.025,0.5,0.975)))
    df_so$x<- seq(1:params$nt)
    
    
    p3<- ggplot(data=df_sf,aes(x=x))+
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
      ylab("IgG seroprevalence") + xlab("day")
    
    
    gridExtra::grid.arrange(p1,p2,p3)
    
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
    df_s<-data.frame(x=x_s,prev=as.numeric(sim$l_prev_age))
    
    p1 <- ggplot(data=df_d,aes(x=x))+ 
      geom_point(aes(y=prev))+
      geom_errorbar(aes(ymin=low, ymax=up), width=.2,position=position_dodge(.9)) + 
      geom_segment(data=df_s, aes(x = x, y=prev , xend = x+0.5, yend = prev), col="red")+
      ylab("Prevalence") +
      scale_x_discrete(name ="Livestock age (yrs)", limits=c("1","","2","","3","","4","","5",""))+
      ylim(c(0,1))
    
    
    
    
    # Human cases 
    df_s<-data.frame(x=seq(1:(params$nt/7)),cases=as.numeric(sim$h_inc_week))
    df_d<-data.frame(x=observations$h_inc_week.x,cases=observations$h_inc_week)
    
    
    p2<- ggplot(data=df_s, aes(x=x))+
      geom_line(aes(y=cases), col="darkgreen", lwd=0.4) +
      geom_point(data =df_d, aes(x=x, y=cases) ) + theme_bw()+
      geom_vline(xintercept=which.max(df_s$cases), color="red") +
      geom_text(x=150, y=8, label="peak (week 74)", size=3, col="red") +
      geom_text(x=160, y=7, label="25/08/2009-31/08/2009", size=3, col="red") +
      ylab("Reported cases in humans")+xlab("week")
    
    
    ## Prevalence farmers and others
  
    df_sf<- data.frame(x=seq(1:params$nt), prev=as.numeric(sim$h_prev_farmer_long))
    df_so<- data.frame(x=seq(1:params$nt), prev=as.numeric(sim$h_prev_other_long))
    
    p3<- ggplot(data=df_sf,aes(x=x))+
      geom_line(aes(y=prev), col="darkgreen", lwd=0.4) +
      geom_line(data=df_so, aes(y=prev), col="blue", lwd=0.4) +
      geom_point(aes(x=field_work_mid,y=prev_farmer$prev))+
      geom_errorbar(aes(x=field_work_mid,ymin=prev_farmer$low_ci, ymax=prev_farmer$up_ci), width=.2,
                    position=position_dodge(.9), col="darkgreen") + 
      geom_point(aes(x=field_work_mid,y=prev_other$prev))+
      geom_errorbar(aes(x=field_work_mid,ymin=prev_other$low_ci, ymax=prev_other$up_ci), width=.2,
                    position=position_dodge(.9),col="blue")+ theme_bw()+
      ylab("IgG seroprevalence") + xlab("day")
    
    gridExtra::grid.arrange(p1,p2,p3)
    
    
    
  }
  
}


