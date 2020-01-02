#' @title kernel density estimation for 2d using R.
#' @description 2d kernel density function estimation by spherical-symmetric kernel function.
#' @param x the sample data of random variable x, it should be a squence.(vector)
#' @param y the sample data of random variable y, it should be a squence.(vector)
#' @param hx the bandwiths of x (numeric)
#' @param hy the bandwiths of y (numeric)
#' @param kernel the kernel function (character)
#' @return a list of x , y and f.hat, (x,y) is coordinate point , f.hat is the estimation of density corresponding to (x,y)
#' @examples
#' \dontrun{
#' data(faithful)
#' res <- kernden2d_ss(faithful$eruptions, faithful$waiting, 0.23, 33, kernel="gaussian")
#' persp(res$x, res$y, res$f.hat, theta = -30, phi = 30)
#' }
#' @export
kernden2d_ss <- function(x, y, hx, hy, kernel) {
  if (length(x) != length(y))
    stop("the length of x, y must be equal")
  
  # kernel function
  if (kernel=="gaussian") 
  {kern <- function(x,y) exp(-(x^2+y^2)/2) / (2*pi)}
  else if (kernel=="epanechnikov")
  {kern <- function(x,y) 9/16*(1-(x^2+y^2)) * (x^2+y^2 <= 1)}
  
  # select the value of x, y corresponding to f.hat
  xx <- seq( min(x)-3*hx, max(x)+3*hx, l = 30)
  yy <- seq( min(y)-3*hy, max(y)+3*hy, l = 30)
  
  ux <- outer(xx, x, "-") / hx
  uy <- outer(yy, y, "-") / hy
  
  m <- matrix(0, nrow = 30, ncol = 30)
  z <- numeric( length(x) )
  for (i in seq_along(xx)) {
    for (j in seq_along(xx)) {
      for (d in seq_along(x)) {
        z[d] <- kern( ux[i,d], uy[j,d])
      }
      m[i,j] <- mean(z) / (hx*hy)
    }
  }
  
  list(x = xx, y= yy, f.hat = m)
}