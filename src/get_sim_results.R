get_sim_results<- function(theta){
  
  nruns<- nrow(theta)
  mo_length<-params$nt
  yr_length<-length(unique(as.Date(temp_month$month,format="%Y")))
  
  # Allocate memory
  h_inc_month <- matrix(NA, nrow=nruns, ncol=mo_length-1)
  h_inc_year <- matrix(NA, nrow=nruns, ncol=yr_length)
  l_prev_age <- matrix(NA, nrow=nruns, ncol=5)
  l_prev_all <- matrix(NA, nrow=nruns, ncol=1)
  sus_l_frac <- matrix(NA, nrow=nruns, ncol=mo_length)
  h_prev_farmer <- matrix(NA, nrow=nruns, ncol=1)
  h_prev_other  <- matrix(NA, nrow=nruns, ncol=1)
  l_prev_all_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  h_prev_farmer_long <- matrix(NA, nrow=nruns, ncol=mo_length)
  h_prev_other_long  <- matrix(NA, nrow=nruns, ncol=mo_length)
  
  
  for (jj in 1:nruns){
    
    
    # Pass calibrated parameters 
    params$theta <- c(A=theta$A[jj], 
                      F_risk=theta$F_risk[jj], 
                      O_factor=theta$O_factor[jj],
                      imm_p=theta$imm_p[jj],
                      RRreport=theta$RRreport[jj],
                      knot1=theta$knot1[jj],
                      knot2=theta$knot2[jj],
                      beta1=theta$beta1[jj])
    # Risk of transmission in other occupations
    params$risk_O <-params$theta[["O_factor"]]*params$theta[["F_risk"]]
    params$imm_t0<-params$imm_t0_bulgaria*params$theta[["imm_p"]]
    # 1. Call model function 
    #########################################################################################################  
    
    out<-get_objective(params,"lsoda")
    
    
    # 2. Process and format model output 
    #########################################################################################################  
    
    
    # #Reported cases in humans
    
    # Reporting vvector 
    # RRyear<- c(rep(params$RR,2013-2008),seq(params$RR, params$theta[["RRreport"]],
    #   length.out=yr_length-(2013-2008)))
    # 
    
    
    
    RRyear<-logistic(seq(1,yr_length),d=yr_length*0.72, 
                  a=2,
                  c=params$RR/params$theta[["RRreport"]], z=1)*params$theta[["RRreport"]]
    
    # plot(seq(1,yr_length),RRyear)
    
    # RRmonth<- c(rep(params$RR,(2013-2008)*12),seq(params$RR, params$theta[["RRreport"]],
    #                                     length.out=(mo_length-(2013-2008)*12)-1))
  
    RRmonth<-logistic(seq(1,mo_length-1),d=(mo_length-1)*0.72, 
                     a=0.2*(mo_length-1)/100,
                     c=params$RR/params$theta[["RRreport"]], z=1)*params$theta[["RRreport"]]
    
   # plot(seq(1,mo_length-1),RRmonth)
    
    
    
      # RRmonth<- seq(params$RR, params$theta[["RRreport"]], length.out=mo_length-1)
    # RRyear<- seq(params$RR, params$theta[["RRreport"]], length.out=yr_length)
    dates<-temp_month$month[2:params$nt]
    tmp<-data.frame(date = dates ,
                     cases_tmp=diff(out$daily_F_incidence+out$daily_O_incidence))
    tmp$cases<- ifelse(tmp$cases_tmp==0, 1e-6, as.numeric(paste(tmp$cases_tmp)))
    tmp$cases_tmp<-NULL
    # 
    human_inc_yr<-tmp%>%
       mutate(year=as.Date(date,format="%Y"))%>%
       group_by(year)%>%
       summarise(reported=sum(cases))%>%
       mutate(cases=ifelse(reported==0,NA,as.numeric(paste(reported))))
    
    h_inc_month[jj,]<- tmp$cases* RRmonth
    h_inc_year[jj,]<-  human_inc_yr$cases * RRyear
    
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
    simu_all_N_long<-c(L_1+L_2+L_3+L_4+L_5)
    l_prev_all_long[jj,]<-simu_all_R/simu_all_N_long
    simu_all_R<-mean(simu_all_R[(field_work_start:field_work_end)])
    simu_all_N<-mean(simu_all_N_long[(field_work_start:field_work_end)])
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
    
    
    # susceptible fraction of livestock
    sus_l_frac[jj,]<-c(out$L_S1+out$L_S2+out$L_S3+out$L_S4+
                        out$L_S5)/simu_all_N_long
    
    
  }
  
  
  
  sim<-list(
    
    h_inc_month=h_inc_month,
    h_inc_year=h_inc_year,
    l_prev_age = l_prev_age,
    l_prev_all = l_prev_all,
    h_prev_farmer = h_prev_farmer,
    h_prev_other  = h_prev_other,
    l_prev_all_long = l_prev_all_long,
    h_prev_farmer_long = h_prev_farmer_long,
    h_prev_other_long  = h_prev_other_long,
    sus_l_frac = sus_l_frac
    # 
  )
  
  
  return(sim)
  
  
  
  
}
