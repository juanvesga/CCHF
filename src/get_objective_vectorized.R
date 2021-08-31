get_objective<-function(params,methd_intgr="euler"){
  
  # Initial conditions
  ###########################################################################################################
  init_prev<-sum(params$imm_t0*params$pop_st*params$N_liv)/params$N_liv
  prev_vec<-rep(0,5) 
  prev_vec[1]<-init_prev
  
  states_ini <- c(
    L_Ri= prev_vec[1]*(1-params$imm_t0[1])*params$pop_st[1]*params$N_liv , 
    L_S = (1-prev_vec)*(1-params$imm_t0)*params$pop_st*params$N_liv - params$pop_st,
    L_I = params$pop_st,
    L_R = params$imm_t0*params$pop_st*params$N_liv,
    
    # Humans
    # Farmers
    F_S = params$N_F - params$N_F*params$prop_F_R,
    O_S =params$N_O - params$N_O*params$prop_O_R,
    F_E = 0,
    O_E = 0,
    F_I = 0,
    O_I = 0,
    F_R = params$N_F*params$prop_F_R,
    O_R = params$N_O*params$prop_O_R,
    #Model outpput
    daily_L_incidence = 0,
    daily_F_incidence = 0,
    daily_O_incidence = 0,
    new_deaths_per_day= 0
    
  )
  
  # Create model transition matrices

  M0 <- make_model(params)
  params$M<-M0
  
  # Seeding infectious livestock each year (transovarial transmission)
  ###########################################################################################################
  #Define vector for seeding transovarial event (every 365 days) 
  modulo_vec<- which((seq(1:params$nt)%%(365))==0)
  seed_vec  <- modulo_vec+1
  
  seed_event<-data.frame(var =rep( c("L_I1","L_I2","L_I3","L_I4","L_I5"),length(seed_vec)),
                         time = rep(seed_vec,each=5) ,
                         value = rep(params$pop_st,length(seed_vec)),
                         method = rep(c("add", "add", "add", "add","add"),length(seed_vec)))
  
  
  # Calculating R value for livestock temperature dependent
  ###########################################################################################################
  soil_data<-data.frame(times=seq(1,params$nt), soil_temp=soil_t)
  
  foi_temp_factor_df<-soil_data%>%
    group_by(times)%>%
    mutate(temp_R_fac= temp_foi_func(soil_temp,params$theta[["A"]]))%>%
    select(times,temp_R_fac)
  
  #Define interp function to find temp factor for t
  foi_temp_factor_func <- approxfun(foi_temp_factor_df, rule = 2)
  params$foi_temp_factor_func<-foi_temp_factor_func
  # Call model function and run with current Theta
  #########################################################################################################  
  
  
  
  fx<- goveqs_basis
  times_new<- seq(1,params$nt)
  
  out <- as.data.frame(ode(y = states_ini, times = times_new, 
                           func = fx, parms = params,method=methd_intgr,
                            events = list(data = seed_event)))
  
  return(out)
}