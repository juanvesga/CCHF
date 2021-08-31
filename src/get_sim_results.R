get_sim_results<- function(theta){
  
  nruns<- nrow(theta)
  day_length<-params$nt
  week_length<-params$nt/7
  
  # Allocate memory
  h_inc_week <- matrix(NA, nrow=nruns, ncol=week_length)
  l_prev_age <- matrix(NA, nrow=nruns, ncol=5)
  l_prev_all <- matrix(NA, nrow=nruns, ncol=1)
  h_prev_farmer <- matrix(NA, nrow=nruns, ncol=1)
  h_prev_other  <- matrix(NA, nrow=nruns, ncol=1)
  l_prev_all_long <- matrix(NA, nrow=nruns, ncol=day_length)
  h_prev_farmer_long <- matrix(NA, nrow=nruns, ncol=day_length)
  h_prev_other_long  <- matrix(NA, nrow=nruns, ncol=day_length)
  
  
  
  for (jj in 1:nruns){
    
    
    # Pass calibrated parameters 
    params$theta <- c(A=theta$A[jj], 
                      F_risk=theta$F_risk[jj], 
                      O_factor=theta$O_factor[jj])
    # Risk of transmission in other occupations
    params$risk_O <-params$theta[["O_factor"]]*params$theta[["F_risk"]]
    
    # 1. Call model function 
    #########################################################################################################  
    
    out<-get_objective(params,"euler")
    
    
    # 2. Process and format model output 
    #########################################################################################################  
    
    
    #Reported cases in humans
    tmp<-data.frame(date = seq( temp_start_date, by=1, len=params$nt-1),
                    time = c(1:(params$nt-1)),group = rep("sim",params$nt-1),
                    cases_tmp=round(diff(out$daily_F_incidence+out$daily_O_incidence)* params$RR,digits = 4),
                    week=rep(seq(params$nt), each=7, length.out=params$nt-1))
    tmp$cases<- ifelse(tmp$cases_tmp==0, 1e-6, as.numeric(paste(tmp$cases_tmp)))
    tmp$cases_tmp<-NULL
    
    human_inc_weekly<-tmp%>%
      group_by(group,week)%>%
      summarise(reported=sum(cases),.groups = 'drop')%>%
      mutate(reported=ifelse(reported==0,NA,as.numeric(paste(reported))))
    
    h_inc_week[jj,]<-human_inc_weekly$reported 
    
    
    # Livestock prevalence by age 
    L_1<-out$L_S1+out$L_I1+out$L_R1+out$L_Ri
    L_2<-out$L_S2+out$L_I2+out$L_R2
    L_3<-out$L_S3+out$L_I3+out$L_R3
    L_4<-out$L_S4+out$L_I4+out$L_R4
    L_5<-out$L_S5+out$L_I5+out$L_R5
    
    simu_age_R<-c(
      mean(out$L_R1[(field_work_start:field_work_end)]+out$L_Ri[(field_work_start:field_work_end)]),
      mean(out$L_R2[(field_work_start:field_work_end)]),
      mean(out$L_R3[(field_work_start:field_work_end)]),
      mean(out$L_R4[(field_work_start:field_work_end)]),
      mean(out$L_R5[(field_work_start:field_work_end)]))
    
    simu_age_N<-c(
      mean(L_1[(field_work_start:field_work_end)]),
      mean(L_2[(field_work_start:field_work_end)]),
      mean(L_3[(field_work_start:field_work_end)]),
      mean(L_4[(field_work_start:field_work_end)]),
      mean(L_5[(field_work_start:field_work_end)]))
    
    l_prev_age[jj,] <- simu_age_R / simu_age_N
    
    simu_all_R<-c(out$L_Ri+out$L_R1+out$L_R2+out$L_R3+out$L_R4+out$L_R5)
    simu_all_N<-c(L_1+L_2+L_3+L_4+L_5)
    l_prev_all_long[jj,]<-simu_all_R/simu_all_N
    simu_all_R<-mean(simu_all_R[(field_work_start:field_work_end)])
    simu_all_N<-mean(simu_all_N[(field_work_start:field_work_end)])
    l_prev_all[jj]<-simu_all_R/simu_all_N
    
    # Human prevalence (farmers and others)
    
    simu_F_R<-  mean(out$F_R[(field_work_start:field_work_end)])
    simu_F_N<-  c(out$F_S + out$F_E + out$F_I + out$F_R)
    simu_F_N<-  mean(simu_F_N[(field_work_start:field_work_end)])
    h_prev_farmer[jj]<- simu_F_R/simu_F_N
    h_prev_farmer_long[jj,]<- out$F_R/c(out$F_S + out$F_E + out$F_I + out$F_R)
    
    
    simu_O_R<-  mean(out$O_R[(field_work_start:field_work_end)])
    simu_O_N<-  c(out$O_S + out$O_E + out$O_I + out$O_R)
    simu_O_N<-  mean(simu_O_N[(field_work_start:field_work_end)])
    h_prev_other[jj]<- simu_O_R/simu_O_N
    h_prev_other_long[jj,]<- out$O_R/c(out$O_S + out$O_E + out$O_I + out$O_R)
    
    
    
  }
  
  
  
  sim<-list(
    
    h_inc_week=h_inc_week,
    l_prev_age = l_prev_age,
    l_prev_all = l_prev_all,
    h_prev_farmer = h_prev_farmer,
    h_prev_other  = h_prev_other,
    l_prev_all_long = l_prev_all_long,
    h_prev_farmer_long = h_prev_farmer_long,
    h_prev_other_long  = h_prev_other_long
    
  )
  
  
  return(sim)
  
  
  
  
}