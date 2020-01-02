#include <Rcpp.h>
using namespace Rcpp;

//' @title A random walk Metropolis sampler using RCPP
//' @description A random walk Metropolis sampler for the standard Laplace distribution using RCPP
//' @param sigma the standard deviation of normal distribution for the increment
//' @param x0 the origial position of the chain
//' @param N the number of generating times
//' @return a random sample for the standard Laplace distribution \code{n}
//' @examples
//' \dontrun{
//' N <- 2000
//' x0 <- 20
//' sigma <- 2
//' rwmC(sigma, x0, N)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwmC(double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0;
  NumericVector u = as<NumericVector>(runif(N));
  
  for(int i = 1; i < N; i++){
    double y = as<double>(rnorm(1, x[i-1], sigma));
    if (u[i] <= (exp(-abs(y))/2) / (exp(-abs(x[i-1])) / 2)){
      x[i] = y;
    }else{
      x[i] = x[i-1];
    }
  }
  return(x);
}
