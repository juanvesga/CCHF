#include <RcppArmadillo.h>


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec compute_model(
    const arma::vec& invec, 
    const arma::mat& lambda, 
    const arma::mat& linear, 
    const arma::mat& nlin_livestock, 
    const arma::mat& nlin_farmer,
    const arma::mat& nlin_other,
    Rcpp::Function func,
    Rcpp::Function func_event,
    const arma::vec& mortvec,
    int n_aux,
    arma::uvec sus_a1_liv_indx,
    arma::uvec sus_farmer_indx,
    arma::uvec sus_other_indx, 
    arma::uvec aux_inc_indx,
    arma::uvec aux_mort_indx,
    arma::uvec livest_indx,
    arma::uvec livest_prev_indx,
    arma::uvec livest_pass_imm_indx,
    arma::uvec farmer_indx,
    arma::uvec other_indx,
    const arma::mat& agg_inc,
    const arma::mat& agg_mort,
    const arma::mat& sel_inc,
    double t) {
  
  // time dependent Functions to num vector
  double N_livestock =sum(invec(livest_indx));
  Rcpp::NumericVector temperature_r_factor = func(t);
  Rcpp::NumericVector spline_event_factor = func_event(t);
  
  // -- Get force of infection
  arma::vec lam = lambda*invec/N_livestock;
  
  // -- Get all model components together
  arma::mat allmat (linear +  
    lam(0) * nlin_livestock * temperature_r_factor[0] +
    lam(1) * spline_event_factor[0]*nlin_farmer +
    lam(2) * spline_event_factor[0]*nlin_other);
  
  
  arma::vec dx = allmat*invec;
  
  // Implement deaths
  arma::vec morts = mortvec%invec;
  dx = dx - morts;
  
  //  Implement births
  
  double liv_prev= sum(invec(livest_prev_indx))/N_livestock;
  
  dx(sus_a1_liv_indx) = dx(sus_a1_liv_indx)+sum(morts(livest_indx))*(1-liv_prev) ;
  dx(livest_pass_imm_indx) = dx(livest_pass_imm_indx)+sum(morts(livest_indx))*liv_prev ;
  
  dx(sus_farmer_indx) = dx(sus_farmer_indx)+sum(morts(farmer_indx)) ;
  dx(sus_other_indx)  = dx(sus_other_indx)+sum(morts(other_indx)) ;
  
  
  // Get the auxiliaries
  arma::vec ou(n_aux, arma::fill::zeros);
  arma::vec out;  
  out= arma::join_cols(dx,ou);
  out(aux_inc_indx)         = agg_inc*(sel_inc%allmat)*invec;
  out(aux_mort_indx)        = agg_mort*morts;
  
  return out;
  
}