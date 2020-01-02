#' @title kernel density estimation for 2d using R.
#' @description kernden2d function is a kernel density estimation for two dimensional, estimating method by product kernel.
#' @param x the sample data of random variable x, it should be a squence.(vector)
#' @param y the sample data of random variable y, it should be a squence.(vector)
#' @param hx the bandwiths of x (numeric)
#' @param hy the bandwiths of y (numeric)
#' @param kernel chose the kernel function, supported kernel function: "gaussian", "uniform", "epanechnikov", "biweight", "triweight". (character)
#' @return a list of x , y and f.hat, (x,y) is coordinate point , f.hat is the estimation of density corresponding to (x,y)
#' @examples
#' \dontrun{
#' data(faithful)
#' res <- kernden2d_pr(faithful$eruptions, faithful$waiting, 0.23, 33, kernel="gaussian")
#' persp(res$x, res$y, res$f.hat, theta = -30, phi = 30)
#' }
#' @export
kernden2d_pr <- function(x, y, hx , hy, kernel) {
  if (length(x) != length(y))
    stop("the length of x, y must be equal")
  
  # kernel function
  if (kernel=="gaussian") 
  {kern <- function(x) dnorm(x)}
  else if (kernel=="uniform")
  {kern <- function(x) 0.5 * (x>=-1 & x<=1)}
  else if (kernel=="epanechnikov")
  {kern <- function(x) (3/4)*(1-x^2)*(x>=-1 & x<=1)}
  else if (kernel=="biweight")     
  {kern <- function(x) (15/16)*(1-x^2)^2*(x>=-1 & x<=1)}
  else if (kernel=="triweight")    
  {kern <- function(x) (35/32)*(1-x^2)^3*(x>=-1 & x<=1)}
 
  # select the value of x, y corresponding to f.hat
  xx <- seq( min(x)-3*hx, max(x)+3*hx, l = 30)
  yy <- seq( min(y)-3*hy, max(y)+3*hy, l = 30)
  
  ux <- outer(xx, x, "-") / hx
  uy <- outer(yy, y, "-") / hy
  
  # kernel weights
  Kx <- kern(ux)
  Ky <- kern(uy)
  
  fx <- (Kx %*% t(Ky)) / (hx*hy*length(x))
  
  list(x=xx, y=yy, f.hat = fx)
}