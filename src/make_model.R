make_model<-function(p){
  
  s<- params$ref$s
  i<- params$ref$i

  
  
  m <- matrix(0,i$nstates,i$nstates)
  
  
  # Lost of passive immunity (~ to average lactation period)
  source <- s$L_Ri; destin <- intersect(s$L_S,s$a1); rate <- p$time_passimm_loss_livestock
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  

  # Livestock immunity after infection
  source <- s$L_I; destin <- s$L_R; rate <- p$recover
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  #Ageing
  source <- s$a1; destin <- s$a2; rate <- p$Ageing[1]
  m[ cbind(destin, source) ] <- m[cbind(destin, source)]+rate
  
  source <- s$a2; destin <- s$a3; rate <- p$Ageing[2]
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  source <- s$a3; destin <- s$a4; rate <- p$Ageing[3]
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  source <- s$a4; destin <- s$a5; rate <- p$Ageing[4]
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  #Human transitions in farmers and others
  source <- s$E; destin <- s$I; rate <- p$time_to_infous
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  source <- s$I; destin <- s$R; rate <- p$time_immune_human
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  source <- s$R; destin <- s$S; rate <- p$time_susceptible_human
  m[ cbind(destin, source) ] <- m[ cbind(destin, source) ]+rate
  
  M<-list()
  M$lin <- m - diag(colSums(m))    
  
  ## Matrices for non-linear transitions (CCHFV transmission) in livestock
  M$nlin<-list()
  
  m <- matrix(0,i$nstates,i$nstates)
  
  source <- s$L_S; destin <- s$L_I; rate <- 1
  m[ cbind(destin, source) ] <- m[ cbind(destin, source)] + rate
  
  M$nlin$l <- m-diag(colSums(m))
  
  mf <- matrix(0,i$nstates,i$nstates)
  mo <- matrix(0,i$nstates,i$nstates)
  
  source <- i$S$farmer; destin <- i$E$farmer; rate <- 1
  mf[ cbind(destin, source) ] <- mf[ cbind(destin, source)] + rate
  M$nlin$f <- mf-diag(colSums(mf))
  
  source <- i$S$other; destin <- i$E$other; rate <- 1
  mo[ cbind(destin, source) ] <- mo[ cbind(destin, source)] + rate
  M$nlin$o  <- mo-diag(colSums(mo))
  
  
  # --- Getting force-of-infection for all groups
  tmp <-matrix(0,3,i$nstates)    # 1.Livestock, 2.Farmer, 3.Other
  tmp[1,s$L_I]  <- 1 - exp(-(1/(p$D_inf_L*p$N_liv)))
  tmp[2,s$L_I]  <- 1 - exp((-p$theta[["F_risk"]]/1e6))
  tmp[3,s$L_I]  <- 1 - exp((-p$risk_O/1e6))
  
  M$lambda <- tmp
  
  # Mortality
  m <- matrix(0,1,i$nstates)
  m[s$a1]<-p$deathd[1]
  m[s$a2]<-p$deathd[2]
  m[s$a3]<-p$deathd[3]
  m[s$a4]<-p$deathd[4]
  m[s$a5]<-p$deathd[5]
  m[s$humans]<-p$b_d
  M$mortvec <- t(m)
  
  return(M)
  
}
