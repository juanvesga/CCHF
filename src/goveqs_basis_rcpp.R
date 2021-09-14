goveqs_basis <- function (t, state, parameters) {
  
  
  with(as.list(c(state,parameters)),             
       {
         i<- ref$i
         s<- ref$s
         sel<- ref$sel
         agg<- ref$agg
         
         invec <- as.numeric(state[1:i$nstates])
         n_aux=i$nx-i$nstates
         
         dx<-compute_model(invec,
                        M$lambda,
                        M$lin,
                        M$nlin$l,
                        M$nlin$f,
                        M$nlin$o,
                        foi_temp_factor_func,
                        foi_event_func,
                        M$mortvec,
                        Birth,
                        Birth_F,
                        Birth_O,
                        n_aux,
                        i$L_S$a1-1,
                        i$S$farmer-1,
                        i$S$other-1,
                        i$aux$inc-1,
                        i$aux$mort-1,
                        s$livestock-1,
                        s$prevalent_l-1,
                        i$L_Ri-1,
                        agg$inc,
                        agg$mort,
                        sel$inc,
                        t)
         list(dx)
         
         
       }
  )
}