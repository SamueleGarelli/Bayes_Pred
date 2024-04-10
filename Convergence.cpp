#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Compute the Cholesky decomposition of 1x1 and 2x2 matrices. In particular, it return an upper-triangular matrix that, if multiplied by its transpose,
// yields the inputed matrix.
// [[Rcpp::export]]
arma::mat cholesky(arma::mat sigma){
  
  if (sigma.n_cols == 1){
    
    return(arma::sqrt(sigma));
    
  }
  
  else{
    
    arma::mat L(2,2);
    
    L(0,0) = sqrt(sigma(0,0) - pow(sigma(0,1),2)/sigma(1,1));
    L(0,1) = sigma(0,1)/sqrt(sigma(1,1));
    L(1,0) = 0;
    L(1,1) = sqrt(sigma(1,1));
    
    return L;
    
  }
  
}

// Return the value of the univariate and bivariate Gaussian density.
// [[Rcpp::export]]
double dens_mvnorm(arma::colvec x, arma::colvec mean, arma::mat cov) {
  
  int k = x.n_rows;
  double det_cov = arma::det(cov);
  
  if (k == 1) {  // For univariate normal
    double sigma = sqrt(cov(0, 0));
    return 1 / (sqrt(2 * M_PI) * sigma) * exp(-0.5 * pow((x(0) - mean(0)) / sigma, 2));
  } else {  // For multivariate normal
    arma::mat diff = x - mean;
    arma::mat inv_cov = arma::inv(cov);
    double exponent = -0.5 * arma::as_scalar(diff.t() * inv_cov * diff);
    double constant = pow(2 * M_PI, -k / 2) * pow(det_cov,-0.5);
    return constant * exp(exponent);
  }
}

// Sample one observation from a multivariate Gaussian distribution.
// [[Rcpp::export]]
arma::rowvec rmvnorm(arma::rowvec mu, arma::mat sigma) {
  
  // Dimension of the multivariate normal distribution
  int p = mu.size();
  
  // Correct the covariance matrix: make it positive definite. This guarantees the existence of a unique Cholesky decomposition
  arma::vec eigenval;
  arma::mat eigenvec;
  arma::eig_sym(eigenval, eigenvec, sigma);
  
  // The arma::chol function fails even with eigenvalues that are positive but very close to zero (usually between 0 and 10^-6).
  // Hence, check the condition "< 10^-6" instead of "< 0"
  if (arma::any(eigenval <= 1e-6)){
    
    arma::uvec indices = arma::find(eigenval <= 1e-6); 
    int s = indices.size();
    for (int i = 0; i < s; ++i) {
      eigenval(indices(i)) = 1e-6;
    }
    
    arma::mat sigma_prime = eigenvec * arma::diagmat(eigenval) * eigenvec.t();
    
    // Create a matrix to store the samples
    arma::rowvec samples(p);
    
    arma::mat L = cholesky(sigma_prime);
    
    // Generate standard normal random samples
    arma::rowvec Z = arma::randn(1,p);
    
    // Transform the standard normal samples to multivariate normal samples
    samples = mu + Z * L;
    
    return samples;
    
  }
  
  else{
    
    // Create a matrix to store the samples
    arma::rowvec samples(p);
    
    // Cholesky decomposition of the covariance matrix
    arma::mat L = cholesky(sigma);
    
    // Generate standard normal random samples
    arma::rowvec Z = arma::randn(1,p);
    
    // Transform the standard normal samples to multivariate normal samples
    samples = mu + Z * L;
    
    return samples;
    
  }
}

// Compute the predictive density and cdf for the Gaussian predictive on a grid of points at each iteration of the predictive resampling algorithm 
// for a given N and B = 1. The
// output is a list of two matrices, the first one containing the density and the the second one containing the cdf.
// [[Rcpp::export]]
List gauss_path_1(int N, arma::vec grid){
  
  // N = size of the population built with predictive resampling and n° of densities/cdf's we are taking
  // grid = grid of points over which we compute the predictive densities/cdf's

  int m = grid.size();
  
  arma::mat p(N,m);
  arma::mat P(N,m);
  
  // compute the predictive density/cdf at the first step and then sample
  p.row(0) = arma::normpdf(grid).t();
  P.row(0) = arma::normcdf(grid).t();
  double x = arma::randn();
  
  // compute the mean, the predictive density/cdf at the second step and then sample (sd is still 1)
  double mu = x;
  p.row(1) = arma::normpdf(grid,mu,1).t();
  P.row(1) = arma::normcdf(grid,mu,1).t();
  x = R::rnorm(mu,1);
  
  double s = sqrt((pow(mu,2)+pow(x,2))/2 - pow((mu+x)/2,2)); // define and compute the sd outside the loop because it is needed inside of it (at this stage mu is just the first
                                                             // observation and x is the second)
  
  for (int i = 2; i < N; ++i){
    
    mu = (mu*(i-1) + x)/i;
    p.row(i) = arma::normpdf(grid,mu,s).t();
    P.row(i) = arma::normcdf(grid,mu,s).t();
    x = R::rnorm(mu,s);
    
    // compute in this iteration, although used in the next one, because we need the "old" mu and not the updated version
    s = sqrt((((pow(s,2) + pow(mu,2))*i + pow(x,2))/(i+1)) -  pow((mu*i + x)/(i+1),2)); 
                                                                          
  }
  
  Rcpp::List result;
  result["p"] = p;
  result["P"] = P;

return result;}

// Same as the above function but with the copula-based predictive distribution. Weights of this predictive distribution are taken as in "Martingale 
// posterior distributions" (Fong, Holmes and Walker), i.e. alpha_n = (2 - 1/n)/(n+1) for n = 1,2,... . Please note that here weights are called
// alpha_n while in our paper we name them r_n.
// [[Rcpp::export]]
List cop_path_1(int N, arma::vec grid, double rho){
  
  // N = size of the population built with predictive resampling and n° of densities/cdf's we are taking
  // grid = grid of points over which we compute the predictive densities/cdf's
  // rho = copula correlation parameter
  
  int m = grid.size();
  
  arma::mat p(N,m);
  arma::mat P(N,m);
  
  // compute p_0/P_0 and then sample from it
  p.row(0) = arma::normpdf(grid).t();
  P.row(0) = arma::normcdf(grid).t();
  double x = arma::randn();
  
  for (int i = 1; i < N; ++i){
    
    // update
    double q_v = R::qnorm5(P(i-1,arma::index_min(abs(grid-x))),0,1,true,false);
    
    for (int c = 0; c < m; ++c){
      
      double q_c = R::qnorm5(P(i-1,c),0,1,true,false);
      
      double numerator = (1/(2*M_PI*sqrt(1-rho*rho)))*exp(-0.5/(1-pow(rho,2))*(pow(q_c,2)+pow(q_v,2)-2*q_c*q_v*rho));
      
      double denominator = R::dnorm4(q_c,0,1,false)*R::dnorm4(q_v,0,1,false); 
      
      p(i,c) = (1-(2-1.0/i)/(i+1))*p(i-1,c)+(2-1.0/i)/(i+1)*(numerator/denominator)*p(i-1,c);
      
      P(i,c) = (1-(2-1.0/i)/(i+1))*P(i-1,c)+((2-1.0/i)/(i+1))*R::pnorm5((q_c - rho*q_v)/sqrt(1-pow(rho,2)),0,1,true,false);
      
    }
    
    // sampling
    double v = arma::randu();
    x = grid(arma::index_min(abs(P.row(i)-v)));
    
  }

  Rcpp::List result;
  result["p"] = p;
  result["P"] = P;
  
return result;}

// The following function performs the same task as the above function but different weights, i.e. we use alpha_n = (2 - 1/n)/((n+1) * (log(n+1))^2)
// for n = 1,2,... .
// [[Rcpp::export]]
List cop_path_1_new_alpha(int N, arma::vec grid, double rho){
  
  // N = size of the population built with predictive resampling and n° of densities/cdf's we are taking
  // grid = grid of points over which we compute the predictive densities/cdf's
  // rho = copula correlation parameter
  
  int m = grid.size();
  
  arma::mat p(N,m);
  arma::mat P(N,m);
  
  // compute p_0/P_0 and then sample from it
  p.row(0) = arma::normpdf(grid).t();
  P.row(0) = arma::normcdf(grid).t();
  double x = arma::randn();
  
  for (int i = 1; i < N; ++i){
    
    // update
    double q_v = R::qnorm5(P(i-1,arma::index_min(abs(grid-x))),0,1,true,false);
    
    for (int c = 0; c < m; ++c){
      
      double q_c = R::qnorm5(P(i-1,c),0,1,true,false);
      
      double numerator = (1/(2*M_PI*sqrt(1-rho*rho)))*exp(-0.5/(1-pow(rho,2))*(pow(q_c,2)+pow(q_v,2)-2*q_c*q_v*rho));
      
      double denominator = R::dnorm4(q_c,0,1,false)*R::dnorm4(q_v,0,1,false); 
      
      p(i,c) = (1-(2 - 1.0/(i+1))/((i+2)*pow(log(i+2),2)))*p(i-1,c)+(2 - 1.0/(i+1))/((i+2)*pow(log(i+2),2))*(numerator/denominator)*p(i-1,c);
      
      P(i,c) = (1-(2 - 1.0/(i+1))/((i+2)*pow(log(i+2),2)))*P(i-1,c)+(2 - 1.0/(i+1))/((i+2)*pow(log(i+2),2))*R::pnorm5((q_c - rho*q_v)/sqrt(1-pow(rho,2)),0,1,true,false);
      
    }
    
    // sampling
    double v = arma::randu();
    x = grid(arma::index_min(abs(P.row(i)-v)));
    
  }
  
  Rcpp::List result;
  result["p"] = p;
  result["P"] = P;
  
  return result;}


