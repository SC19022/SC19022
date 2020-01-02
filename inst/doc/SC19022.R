## ------------------------------------------------------------------------
kernden2d_pr <- function(x, y, hx , hy, kernel) {
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

## ------------------------------------------------------------------------
data("faithful")
res <- kernden2d_pr(faithful$eruptions, faithful$waiting, 0.23, 33, kernel="gaussian")
persp(res$x, res$y, res$f.hat, theta = -30, phi = 30, 
      col = "lightblue", xlab = "eruptions", zlab = "kde",
      ylab = "waiting", main = "product kernel")

## ------------------------------------------------------------------------
kernden2d_ss <- function(x, y, hx, hy, kernel) {
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

## ------------------------------------------------------------------------
data("faithful")
res <- kernden2d_ss(faithful$eruptions, faithful$waiting, 0.23, 33, kernel="gaussian")
persp(res$x, res$y, res$f.hat, theta = -30, phi = 30, 
      col = "red", xlab = "eruptions", zlab = "kde",
      ylab = "waiting", main = "product kernel")

## ---- eval=FALSE---------------------------------------------------------
#  NumericVector rwmC(double sigma, double x0, int N) {
#    NumericVector x(N);
#    x[0] = x0;
#    NumericVector u = as<NumericVector>(runif(N));
#  
#    for(int i = 1; i < N; i++){
#      double y = as<double>(rnorm(1, x[i-1], sigma));
#      if (u[i] <= (exp(-abs(y))/2) / (exp(-abs(x[i-1])) / 2)){
#        x[i] = y;
#      }else{
#        x[i] = x[i-1];
#      }
#    }
#    return(x);
#  }

## ------------------------------------------------------------------------
x<-c(16,19,18,16,20)
p<-c(0.2,0.15,0.15,0.3,0.2)
EX<-sum(x*p)
VarX<-sum(x^2*p)-EX^2

## ------------------------------------------------------------------------
height<-c(160,175,172)
weight<-c(55,66,62)
plot(weight,height)

## ------------------------------------------------------------------------

n1<-c("jack","peter","fisher")
n2<-c("age","height/cm","weight/kg")
a<-c(18,17,20,160,175,172,55,66,62)
m<-matrix(a, byrow=F, nrow=3,dimnames = list(n1,n2))
knitr::kable(m,format = "markdown")


## ------------------------------------------------------------------------

set.seed(100)
n<-1000; sigma<-c(0.25,1,9,25);
#different value of σ
u<-matrix(NA,nrow = 4,ncol = 1000)
y<-runif(n)
#since the CDF F(x) obey to U[0,1]

for (i in 1:4) {
  u[i,]<-sigma[i]*(-2*log(1-y))^(0.5)
  #the inverse function of F(u)
}
for (i in 1:4) { 
  xlength<-c(2,6,50,150);interval<-c(.01,.05,.3,.8)
  x<-seq(0,xlength[i],interval[i])
  hist(u[i,],breaks = x,freq = FALSE,main= c("σ=",paste(sigma[i])),sub= expression(f(x)==frac(x,sigma^2)*e^-frac(x^2,2*sigma^2)),xlab=NULL,col = "aquamarine",border = "aquamarine2")
  lines(x,(x/sigma[i]^2)*exp(-x^2/(2*sigma[i]^2)),col ="blue")
}

## ------------------------------------------------------------------------
set.seed(100)
n<-1000;
x1<-rnorm(n,0,1);x2<-rnorm(n,3,1)
z<-matrix(NA,nrow = 9,ncol = n)
for (i in 1:9) {
  k<-seq(0,1,.125)
  p1<-sample(c(0,1),size = n,replace = TRUE,prob = c(1-k[i],k[i]))
  z[i,]<-p1*x1+(1-p1)*x2#generate the random number from mixture distribution
  d<-seq(-5,8,.1)
  hist(z[i,],breaks = d,prob=TRUE,col = "skyblue",border = "skyblue",main= c("p1=",paste(k[i])),xlab = NULL)
  lines(d,(1/sqrt(2*pi))*(k[i]*exp(-(d^2)/2)+(1-k[i])*exp(-((d-3)^2)/2)),col="blueviolet")#the pdf of mixture distribution.
}

## ------------------------------------------------------------------------
wishart<-function(n,d,sigma){
  L<-t(chol(sigma))#cholesky decomposition
  T<-matrix(0,d,d)
  for (i in 1:d) {
    T[i,i]<-sqrt(rchisq(1,n-i+1))
    j<-1
    while(i>j){
      T[i,j]<-rnorm(1,0,1)
      j<-j+1
    }
  }
  y<-L%*%T%*%t(T)%*%t(L)#the expression of r.v.
 y
}

#example
sigma<-matrix(rbind(c(5,2,1),c(2,5,2),c(1,2,4)),3,3)
x<-wishart(5,3,sigma)
x

## ------------------------------------------------------------------------
# We know that X~U(a,b), so to find the integration of sin(x) in [a,b] equals to the expectation of (b-a)*sin(x)
set.seed(100)
n<-1e5;a<-0;b<-pi/3
x<-runif(n,a,b)
mc.val<-mean((b-a)*sin(x))
ex.val<--cos(b)+cos(a)
print(c(mc.val,ex.val))

## ------------------------------------------------------------------------
set.seed(100)
m<-1e4
#simple monte carlo integration
x<-runif(m)
g1<-exp(-x)/(1+x^2)
mc1<-mean(g1)
#antithetic variables
u<-runif(m/2);v<-1-u
g2<-(exp(-u)/(1+(u)^2)+exp(-v)/(1+(v)^2))/2
mc2<-mean(g2)
#reduction in variance as a percentage
perc<-paste(round(100*(var(g1)-var(g2))/var(g1),2),"%")
print(c(mc1,mc2,perc))

## ------------------------------------------------------------------------
n<-10000
esti<-sd<-numeric(5)
g<-function(x){
  exp(-x-log(1+x^2))
}
for (i in 1:5) {
    x<-runif(n/5,(i-1)/5,i/5)
    fg<-g(x)/(exp(-x)/(1-exp(-1)))
    esti[i]<-mean(fg)
    sd[i]<-sd(fg)
  }
estimates<-mean(esti)
var<-mean(sd)
print(c(estimates,var))

## ------------------------------------------------------------------------
# we use UCL to express the upper confidence limit
# LCL expresses the lower confidence limit.
# n is the sample size.
# m is the the numbers of replicates.
# 1-alpha is the confidence level.
# prob is the probability of confidence interval cover
# mean between normal sample and chisquare sample.

n     <- 20
m     <- 1000
alpha <- .05
prob  <- numeric(2)
LCL.chisq <- UCL.chisq <-
LCL.norm <- UCL.norm <- numeric(m)

f.u <- function(x){
  mean(x) + qt(1-alpha/2, df=n-1) * sqrt(var(x)/n)
}
f.l <- function(x){
  mean(x) - qt(1-alpha/2, df=n-1) * sqrt(var(x)/n)
}

for (i in 1:m) {
  x1 <- rchisq(n, df = 2)
  UCL.chisq[i] <- f.u(x1)
  LCL.chisq[i] <- f.l(x1)
  
  x2 <- rnorm(n, mean = 2, sd =2)
  UCL.norm[i] <- f.u(x2)
  LCL.norm[i] <- f.l(x2)
}
prob[1] <- mean(LCL.chisq < 2 & 2 < UCL.chisq)
prob[2] <- mean(LCL.norm < 2 & 2 < UCL.norm)
print(prob)

## ------------------------------------------------------------------------
# m is the number of replicates
# n is the sample size of random variable from normal.
# k is the number of quantiles
# we define "quant" function to calculate the quantile
# of skewness.
# "se" is the function to calculate the stand error.

m <- 1000
n <- 100
k <- 100
y <- numeric(k)

quant <- function(q){
  y <- numeric(k)
  for (i in 1:k) {
    SK <- replicate(n, expr = {
      x <- rnorm(m, mean = 0, sd = 1)
      mean(((x-mean(x))/sd(x))^3)
    })
    y[i] <- as.numeric(quantile(SK, q))
  }
  mean(y)
}

se <- function(q){
  sqrt(q*(1-q) / n*(dnorm(0,1))^2)
}

rname <- c("large.q","est.q","se")
cname <- c("q=0.025", "q=0.05", "q=0.95", "q=0.975")
results <- matrix(NA, 3, 4, dimnames = list(rname,cname))
results[1, 1:4] <- c(qnorm(c(.025, .05, .95, .975), 0, 6/n))
results[2, 1:4] <- c(quant(.025), quant(.05), quant(.95), quant(.975))
results[3, 1:4] <- c(se(.025), se(.05), se(.95), se(.975))
print(results)

## ------------------------------------------------------------------------
skew <- function(x) {
  # caculates the sample skewness coefficient.
  y <- mean(x)
  return( mean((x - y)^3) / mean((x - y)^2)^1.5 )
}

beta  <- seq(1, 100, 5)  # this is parameter of beta distribution.
m     <- 10000           # this is number of replicates.
n     <- 50              # this is the sample size.
alpha <- .05             # this is the significance level.
cv    <- qnorm(.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
                         # this is the critical value in test.
pwr   <- numeric( length(beta) )
                         # pwr is the value of value.
for (i in 1:length(beta)) {
  reject <- replicate(m, expr = {
    # simulate x from beta(beta, beta).
    x <- rbeta(n, beta[i], beta[i])
    as.integer( abs(skew(x)>cv) )
  })
  pwr[i] <- mean(reject)
}
plot(beta, pwr, type = "b", xlab = "beta of beta(beta, beta)", 
     ylab = "power", ylim = c(0, .05))
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(beta, pwr+se, lty = 3)
lines(beta, pwr-se, lty = 3)

# for t(theta)
theta <- c(seq(1, 100,5))
pwr   <- numeric( length(theta) )
for (i in 1:length(theta)) {
  reject <- replicate(m, expr = {
    # simulate x from t distribution.
    x <- rt(n, df = theta)
    as.integer( abs(skew(x)>cv) )
  })
  pwr[i] <- mean(reject)
}
plot(theta, pwr, type = "b", xlab = "degree of freedom in t", 
     ylab = "power", ylim = c(0, .5))
abline(h = .05, lty = 2)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(theta, pwr+se, lty = 3)
lines(theta, pwr-se, lty = 3)

## ------------------------------------------------------------------------
alpha <- .05  # this is the significance level.
n <- 100      # this is the sample size.
m <- 10000    # this is the number of replicates.
p <- matrix(0, nrow = 3, ncol = m)
m.data <- matrix(0, nrow = 3, ncol = n)
for (i in 1:m) {      # m replicates
  m.data[1, ] <- rchisq(n, df = 1)
  m.data[2, ] <- runif(n, min = 0, max = 2)
  m.data[3, ] <- rexp(n, rate = 1)
  for (j in 1:3) {    # t.test for chisquare, uniform, 
    # and exponential distributions.
    ttest <- t.test(m.data[j, ], alternative = "two.sided", mu = 1)
    p[j,i] <- as.numeric(ttest$p.value < alpha)
  }
}
t1e <- apply(p, 1, mean)
print(t(t1e))

## ------------------------------------------------------------------------
library("bootstrap")
library("boot")

data(scor)
data <- as.matrix(scor)
color <- c("red", "blue", "yellow", "pink", "green")



plot(data[ ,1], col = color[1], pch = 20, ylim = c(0,100),
     xlab = "student", ylab = "score", main = "scatter plots")
for (i in 2:5) {     # draws the scatter plots.
  points(data[ ,i], col = color[i], pch = 20)
}
legend(80,100, c("mec","vec","alg","ana","sta"),
       col = color, pch =20)



subject <- data.frame(scor)
cor.est <- cor(subject, subject)

r <- function(x, i){    # computers the correlation coefficient.
  cor(x[i, 1], x[i,2])
}
r12 <- boot(data = data[ ,1:2], statistic = r, R = 2000)
r34 <- boot(data = data[ ,3:4], statistic = r, R = 2000)
r35 <- boot(data = data[ ,c(3,5)], statistic = r, R = 2000)
r45 <- boot(data = data[ ,4:5], statistic = r, R = 2000)


knitr::kable(round(cor.est, 4))


se.hat <- t(c(sd(r12$t), sd(r34$t), sd(r35$t), sd(r45$t)))
colnames(se.hat) <- c("r12se.hat", "r34se.hat", "r35se.hat", "r45se.hat")
knitr::kable(se.hat, digits = 4)

## ------------------------------------------------------------------------
library("boot")


skew <- function(x, i) {
  # caculates the sample skewness coefficient.
  y <- mean(x[i])
  return( mean((x[i] - y)^3) / mean((x[i] - y)^2)^1.5 )
}



cr <- function(distribution, alpha, r){
  # this function to compute the coverage rates of standard normal
  # bootstrap confidence interval (SNBCI), basic bootstrap confidence
  # interval (BBCI), and percentile confidence interval (PCI) for 
  # sample skewness statistic.
  
  # alpha is the confidence level
  # r is the times of bootstrap resampling.
  
  SNBCI.ucl <- SNBCI.lcl <- BBCI.ucl <- 
  BBCI.lcl <- PCI.ucl <- PCI.lcl <- numeric(r)
  
  for (i in 1:r) {
    if (distribution == "normal")     {x = rnorm(100)}
    else if (distribution == "chisq") {x = rchisq(100, df = 5)}
    
    # SB is the estimated statistic of the bootstrap resampling
    # here, we define the three interval's confidence limits.
    
    SB <- boot(data = x, statistic = skew, R=1000)
    SNBCI.ucl[i] <- SB$t0 - qnorm(alpha/2)*sd(SB$t)
    SNBCI.lcl[i] <- SB$t0 + qnorm(alpha/2)*sd(SB$t)
    BBCI.ucl[i] <- 2*SB$t0 - as.numeric(quantile(SB$t, alpha/2))
    BBCI.lcl[i] <- 2*SB$t0 - as.numeric(quantile(SB$t, 1-alpha/2))
    PCI.ucl[i] <- as.numeric(quantile(SB$t, 1-alpha/2)) 
    PCI.lcl[i] <- as.numeric(quantile(SB$t, alpha/2))
  }
  
  
  
  if (distribution == "normal")     {sktrue = 0}
  else if (distribution == "chisq") {sktrue = sqrt(8/5)}
  l1 <- mean( SNBCI.lcl >= sktrue )
  r1 <- mean( SNBCI.ucl <= sktrue )
  l2 <- mean( BBCI.lcl  >= sktrue )
  r2 <- mean( BBCI.ucl  <= sktrue )
  l3 <- mean( PCI.lcl   >= sktrue )
  r3 <- mean( PCI.ucl   <= sktrue )
  cr1 <- 1 - l1 - r1
  cr2 <- 1 - l2 - r2
  cr3 <- 1 - l3 - r3
  list(c(SNBCI.cr = cr1, SNBCI.leftmiss = l1,
         SNBCI.rightmiss = r1),
       c(BBCI.cr = cr2, BBCI.leftmiss = l2,
         BBCI.rightmiss = r2),
       c(PCI.cr = cr3, PCI.leftmiss = l3,
         PCI.rightmiss = r3))
}
normal <- cr("normal", .05, 1000)
knitr::kable(rbind(unlist(normal)))


chisq <- cr("chisq", .05, 1000)
knitr::kable(rbind(unlist(chisq)))



## ------------------------------------------------------------------------
library(bootstrap)
m <- data.matrix(scor)
n <- nrow(m)

lambda.hat <- eigen(t(m)%*%m/n)$values
theta.hat <- lambda.hat[1] / sum(lambda.hat)

### jackknife
theta.j <- numeric(n)
for (i in 1:n) { # computes the estimator of theta by jackknife
  m.j <- m[-i,]
  lambda.j <- eigen(t(m.j)%*%m.j/(n-1))$values
  theta.j[i] <- lambda.j[1] / sum(lambda.j)
}
bias.est <- (n-1) * ( mean(theta.j)-theta.hat )
se.est <- sqrt( ((n-1)/n) * sum((theta.j-mean(theta.j))^2) )
print(c(bias.est = bias.est, se.est = se.est))

## ------------------------------------------------------------------------
library(DAAG); attach(ironslag)

n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)
# for n-fold cross validation
# fit models on leave-one-out samples
for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  z1 <- lm(y ~ x)      # linear model
  yhat1 <- z1$coef[1] + z1$coef[2] * chemical[k]
  e1[k] <- magnetic[k] - yhat1
  
  z2 <- lm(y ~ x + I(x^2)) # quadratic model
  yhat2 <- z2$coef[1] + z2$coef[2] * chemical[k] +
    z2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k] - yhat2
  
  z3 <- lm(log(y) ~ x)    # exponential model
  logyhat3 <- z3$coef[1] + z3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k] - yhat3
  
  z4 <- lm(y ~ x + I(x^2) + I(x^3))
                          # cubic polynomial model
  yhat4 <- z4$coef[1] + z4$coef[2] * chemical[k] +
    z4$coef[3] * chemical[k]^2 + z4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}

round(c(mean(e1^2), mean(e2^2), mean(e3^2), mean(e4^2)), 4)

round(c(summary(z1)$adj.r.squared, summary(z2)$adj.r.squared,
  summary(z3)$adj.r.squared, summary(z4)$adj.r.squared), 4)

## ------------------------------------------------------------------------
## p.test is a permutation test for the different size
## we suppose x, y are drawn from normal distribution,
## and with the same mean, but different variances and
## different sample sizes

p.test <- function(n1, n2, mu1 ,mu2 ,sigma1, sigma2) {
  x <- rnorm(n1, mu1, sigma1)
  y <- rnorm(n2, mu2, sigma2)    # get the normal sample.
  count <- function(x, y) {      # count the number 
                                 # of extreme points.
    
    X <- x - mean(x)
    Y <- y - mean(y)             # centralization
    countx <- sum(X > max(Y)) + sum(X < min(Y))
    county <- sum(Y > max(X)) + sum(Y < min(X))
    return( max(c(countx,county)) )
  }
  
  
  T.hat <- count(x, y)           # this is the statistic for
                                 # the original sample.
  
  
  r <- 99                      # resampling times.
  T <- numeric(r)
  z <- c(x, y)
  
  for (i in 1:r) {               # compute the statistic 
                                 # for r times resampling.
    k <- sample(1:(n1 + n2), size = n1, replace = FALSE)
    x <- z[k]
    y <- z[-k]
    T[i] <- count(x, y)
  }
  
  mean(c(T.hat, T) >= T.hat)
}



set.seed(1)
p1 <- p.test(50, 70, 0, 0, 1, 1)
p2 <- p.test(50, 70, 0, 0, 1, 2)
print(c(equal.asl = p1, diff.asl = p2))

## ------------------------------------------------------------------------
library(MASS)
library(Ball)
library(energy)


n <- seq(20, 100, l=5)         # sample size
mu <- rep(0, 2)                # this is the mean
sigma <- diag(2)               # this is the variance
alpha <- .1                    # the significance level


p1 <- p2 <- numeric(length(n))

## here, we computer the model 1 firstly.
for (i in 1:length(n)) {
  
  pv1 <- replicate(100, expr = {      # get pvalue
    x <- mvrnorm(n[i], mu, sigma)     # sampling x
    e <- mvrnorm(n[i], mu, sigma)     # sampling e
    y <- x / 4 + e
    dcov.test(x, y, R = 99)$p.value})
  
  p1[i] <- mean( pv1 < alpha )        # this is the power
  
  
  
  pv2 <- replicate(100, expr = {
    x <- mvrnorm(n[i], mu, sigma)     # sampling x
    e <- mvrnorm(n[i], mu, sigma)     # sampling e
    y <- x / 4 + e
    as.numeric(bcov.test(x, y, R = 99, seed = NULL)$p.value)})
                     
  p2[i] <- mean( pv2 < alpha)
}


plot(n, p1, type = "b", xlab = "n", ylab = "power",
     main = expression(Y == X/4 + e), col = "blue", lty = 3)
lines(n, p2, type = "b", col = "red", lty = 3)
legend(5, 1, paste(c("dcov","Ball")), col = c("blue","red"),
       lty = c(3,3), cex = .7)




## then, repeat the process to computer the model 2

n <- seq(10, 210, l=5)
p_1 <- p_2 <- numeric(length(n))

for (i in 1:length(n)) {
  
  pv_1 <- replicate(100, expr = {
    x <- mvrnorm(n[i], mu, sigma)     # sampling x
    e <- mvrnorm(n[i], mu, sigma)     # sampling e
    y <- (x / 4) * e
    dcov.test(x, y, R = 99)$p.value})
  
  p_1[i] <- mean( pv_1 < alpha )
  
  
  pv_2 <- replicate(100, expr = {
    x <- mvrnorm(n[i], mu, sigma)     # sampling x
    e <- mvrnorm(n[i], mu, sigma)     # sampling e
    y <- (x / 4) * e
    as.numeric(bcov.test(x, y, R = 99, seed = NULL)$p.value)})
  
  p_2[i] <- mean( pv_2 < alpha)
}

plot(n, p_1, type = "b", xlab = "n", ylab = "power",
     main = expression(Y == X/4 %*% e), col = "blue", lty = 3)
lines(n, p_2, type = "b", col = "red", lty = 3)
legend(150, .5, paste(c("dcov","Ball")), col = c("blue","red"),
       lty = c(3,3), cex = 1)

## ------------------------------------------------------------------------
set.seed(1)

df <- function(x) {    # df is the density function
                       # of standard laplace distribution.
  f <- exp(-abs(x)) / 2
  return(f)
}

quan <- function(q) {  # quan is the func to find the quantile
                       # of standard laplace distribution.
  
  df <- function(x) {
    f <- exp(-abs(x)) / 2
    return(f)
  }
  
  cumprob <- function(y) {
    f <- integrate(df, -Inf, y)$value - q
    return(f)
  }
  uniroot(cumprob, c(-10, 10))$root
}




rw.Metropolis <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (df(y) / df(x[i-1])))
      x[i] <- y  
    else {
      x[i] <- x[i-1]
      k <- k + 1
    }
  }
  return(list(x=x, k=k))
}





N <- 2000
sigma <- c(.05, .5, 2, 16)
x0 <- 20

rw1 <- rw.Metropolis(sigma[1], x0, N)
rw2 <- rw.Metropolis(sigma[2], x0, N)
rw3 <- rw.Metropolis(sigma[3], x0, N)
rw4 <- rw.Metropolis(sigma[4], x0, N)

rw <- cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
for (j in 1:4) {
  plot(rw[,j], type="l",
       xlab=bquote(sigma == .(round(sigma[j],3))),
       ylab="X", ylim=range(rw[,j]))
  abline(h = c(quan(.025), quan(.975)))
}

a <- c(.05, seq(.1, .9, .1), .95)

Q <- numeric(length(a))
for (i in 1:length(a)) {
  Q[i] <- quan(a[i])
}

mc <- rw[501:N, ]

Qrw <- apply(mc, 2, function(x) quantile(x, a))
qq <- data.frame(round(cbind(Q, Qrw), 3))
names(qq) <- c('True','sigma=0.05',
               'sigma=0.5','sigma=2','sigma=16')

knitr::kable(qq)

rate <- (N - c(rw1$k, rw2$k, rw3$k, rw4$k)) / N
acceptance.rate <- data.frame(sigma = sigma, 
                              acceptance.rate = rate)
knitr::kable(acceptance.rate)


## ------------------------------------------------------------------------
x <- 0.1
ex1 <- exp(log(x))
ex2 <- log(exp(x))
isTRUE(ex1 == ex2)
isTRUE(all.equal(ex1, ex2))

## ------------------------------------------------------------------------
b <- c(4:25, 100, 500, 1000) # this is the freedom degree.
m <- matrix(NA, nrow = 25, ncol = 3)

f <- function(a, k) {
  pt(sqrt(a^2*k/(k+1-a^2)), df = k, log.p = TRUE) -
    pt(sqrt(a^2*(k-1)/(k-a^2)), df = k-1, log.p = TRUE)}


for (i in 1:25) {
  k <- b[i]
  xx <- seq(1e-5, sqrt(b[i])-1e-5, .01)
  g <- function(a) f(a, b[i])
  yy <- apply(as.matrix(xx), 1, g)
  plot(xx, yy, type = "l", main = paste("k=",b[i]),
       xlab = "a", ylab = "f")
  abline(h = 0, col = "red",lty=3)
}

## ------------------------------------------------------------------------
ak1 <- function(k) {
 # to computer the intersection point in exercise 11.4
  res <- uniroot(f, c(1e-5, 2), k=k)
  res$root
}

m[ ,1] <- apply(as.matrix(b), 1, ak1)



F <- function(a, k) {
  # this function is to calculate s(a, k) - s(a, k+1)
  ck <- function(a, k) sqrt( a^2*k/(k+1-a^2))
  g <- function(u) (1+u^2/(k-1))^(-k/2)
  int <- function(a, k) integrate(g, lower = 0, upper = ck(a, k))
  
  s <- function(a, k) {
    exp( log(2/sqrt(pi*(k-1))) + lgamma(k/2)
         - lgamma((k-1)/2) +log(int(a,k)$value))
  }
  
  s(a, k) - s(a, k+1)
}


for (i in 1:25) {
  k <- b[i]
  xx <- seq(1e-5, sqrt(b[i])-1e-5, 0.1)
  G <- function(a) F(a, b[i])
  yy <- apply(as.matrix(xx), 1, G)
   plot(xx, yy, type = "l", main = paste("k=",b[i]),
       xlab = "a", ylab = "F")
  abline(h = 0, col = "red")
}

## ------------------------------------------------------------------------
ak2 <- function(k) {
  # to computer the intersection point in (0,2)
  res <- uniroot(F, c(1e-5, 2), k=k)
  res$root
}
m[ ,2] <- apply(as.matrix(b), 1, ak2)

ak3 <- function(k) {
  # to computer the intersection point
  res <- uniroot(F, c(1.8, 2.5), k=k)
  res$root
}

m[3:25,3] <- apply(as.matrix(b[3:25]), 1, ak3)
m[2,3] <- uniroot(F, c(1.8, sqrt(5)), k=5)$root
colnames(m) <- c("intersection point in 11.4", 
                 "1st intersection point in 11.5",
                 "2nd intersection point in 11.5")

knitr::kable(m)

## ------------------------------------------------------------------------
N <-100
tol <- .Machine$double.eps^0.5   # the digit of converge
na <- 24
nb <- 28
noo <- 41
nab <- 70
n <- sum(na, nb, noo, nab)

p <- q <- r <- Elogl <- rep(0, N)
p[1] <- q[1] <- 1/3              # the original p, q ,r
k <- 0
for (i in 2:N) {
  p.old <- p[i-1]
  q.old <- q[i-1]
  r.old <- 1-p.old-q.old
  
  naa <- na * p.old^2 / (p.old^2 + 2*p.old*r.old)
  nbb <- nb * q.old^2 / (q.old^2 + 2*q.old*r.old)
  nao <- na * 2 * p.old * r.old / (p.old^2 + 2*p.old*r.old)
  nbo <- nb * 2 * q.old * r.old / (q.old^2 + 2*q.old*r.old)
  
  p[i] <- sum(2*naa, nab, nao) / (2*n)
  q[i] <- sum(2*nbb, nab, nbo) / (2*n)
  r[i] <- sum(2*noo, nao, nbo) / (2*n)
  
  if (sum( abs(p[i]-p.old)/p.old < tol))
      break
  
  Elogl[i-1] <- naa*2*log(p[i]) + nbb*2*log(q[i]) + 
    noo*2*log(r[i]) + nab*log(2*p[i]*q[i]) +
    nao*log(2*p[i]*r[i]) +
    nbo*log(2*q[i]*r[i])
  prob.hat <- c(p[i], q[i], r[i])
  k <- k + 1
}

print(round(c(p.hat=prob.hat[1], q.hat=prob.hat[2]),3))
plot(1:k, Elogl[1:k], type = "b",
     xlab = "the recursive times", ylab = "Elogl")

## ------------------------------------------------------------------------
disp <- mtcars$disp
mpg <- mtcars$mpg
wt <- mtcars$wt


formulas <- list(
  model1 = mpg ~ disp,
  model2 = mpg ~ I(1 / disp),
  model3 = mpg ~ disp + wt,
  model4 = mpg ~ I(1 / disp) + wt
  )

# by using the lapply function.
lapply(formulas, function(x) lm(x)$coefficients)

# by using the for loop.
# we creat a list to save the results firstly.

forloops <- coefficient <- list(model1="1", model2="2",
                            model3="3", model4="4")
for (i in 1:4) {
  forloops[[i]] <- lm(formulas[[i]])
  coefficient[[i]] <- forloops[[i]]$coefficients
}
coefficient

## ------------------------------------------------------------------------

bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})


# by using lapply.
temp <- lapply(bootstraps, function(x) lm(x$mpg ~ x$disp)$coefficients)
m1 <- vapply(temp, as.numeric, matrix(0, nrow = 2))
colnames(m1) <- "fit"; rownames(m1) <- c("intercept", "slope")
knitr::kable(round(m1, 4))


# by using a for loop.
forloop <- coefficient2 <- list(fit1="1", fit2="2", fit3="3",
                                fit4="4", fit5="5",fit6="6", 
                                fit7="7",fit8="8",fit9="9", fit10="10")
for (i in 1:10) {
  x <- bootstraps[[i]]
  forloop[[i]] <- lm(x$mpg ~ x$disp)
  coefficient2[[i]] <- forloop[[i]]$coefficients
}
temp <- coefficient2
m2 <- vapply(temp, as.numeric, matrix(0, nrow = 2))
colnames(m2) <- "model"; rownames(m2) <- c("intercept", "slope")
knitr::kable(round(m2, 4))

# without using an anonymous function.
fitfunc <- function(n){
  resmodel <- rescoefficient <- forloop
  for (i in 1:n) {
    rows <- sample(1:nrow(mtcars), rep = TRUE)
    data <- mtcars[rows, ]
    resmodel[[i]] <- lm(data$mpg ~ data$disp)
    rescoefficient[[i]] <- resmodel[[i]]$coefficients
  }
  list(model = resmodel, coefficient = rescoefficient)
}

res <- fitfunc(10)
temp <- res$coefficient
m3 <- vapply(temp, as.numeric, matrix(0, nrow = 2))
colnames(m3) <- "model"; rownames(m3) <- c("intercept", "slope")
knitr::kable(round(m3, 4))

## ------------------------------------------------------------------------
## for exercise 11.1.3
rsq <- function(mod) summary(mod)$r.squared
r1 <- unlist(lapply(forloops, rsq))
knitr::kable(round(t(r1), 4))

## for exercise 11.1.4 with using anonymous function.
r2 <- unlist(lapply(forloop,rsq))
knitr::kable(round(t(r2), 4))

## for exercise 11.1.4 without anonymous function.
r3 <- unlist(lapply(res$model, rsq))
knitr::kable(round(t(r3), 4))

## ------------------------------------------------------------------------
set.seed(1)
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)

# by using sapply and anonymous function.
sapply(trials, function(x) x$p.value)

# by using sapply and '[[' directly.
sapply(trials, '[[', 3)

## ------------------------------------------------------------------------
## recall the additional homework in 11/22, the solution is really 
## taking much time. For this reason, we attempt to write a code 
## with parallel computing by parSapply function. We suppose five
## sampling size (5,10,15,20,25) for the model Y = X/4+e, computing
## the p-value by sapply and parSapply respectively, and compare 
## with the results, by the way, it also contains the consuming 
## time for the consequence.


library(parallel)
library(MASS)
library(energy)

p.value <- function(n) {
  mu <- rep(0, 2)
  sigma <- diag(2) 
  alpha <- .1 
  
  pv1 <- replicate(1000, expr = {
    x <- MASS::mvrnorm(n, mu, sigma)
    e <- MASS::mvrnorm(n, mu, sigma)
    y <- x / 4 + e
    energy::dcov.test(x, y, R = 199)$p.value})
  
  p <- mean( pv1 < alpha )
  return(p)
}

n <- seq(5, 25, 5)
sapply(n, p.value)
cl <- makeCluster(4)
parSapply(cl, n, p.value)
system.time(sapply(n, p.value))
system.time(parSapply(cl, n, p.value))
stopCluster(cl)

## ------------------------------------------------------------------------
library(Rcpp)
library(microbenchmark)
set.seed(1)

# df is the density function
# of standard laplace distribution.
df <- function(x) {    
  f <- exp(-abs(x)) / 2
  return(f)
}

rwmR<- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= (df(y) / df(x[i-1])))
      x[i] <- y  
    else {
      x[i] <- x[i-1]
    }
  }
  return(x)
}


cppFunction(
'NumericVector rwmC(double sigma, double x0, int N){
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
}')


N <- 2000
x0 <- 20
sigma <- c(.05, .5, 2, 16)

rwmR1 <- rwmR(sigma[1], x0, N)
rwmR2 <- rwmR(sigma[2], x0, N)
rwmR3 <- rwmR(sigma[3], x0, N)
rwmR4 <- rwmR(sigma[4], x0, N)

rwmC1 <- rwmC(sigma[1], x0, N)
rwmC2 <- rwmC(sigma[2], x0, N)
rwmC3 <- rwmC(sigma[3], x0, N)
rwmC4 <- rwmC(sigma[4], x0, N)

rwmR <- cbind(rwmR1, rwmR2, rwmR3,  rwmR4)
rwmC <- cbind(rwmC1, rwmC2, rwmC3,  rwmC4)


for (j in 1:4) {
  par(mfrow=c(1,2))
  qqnorm(rwmR[ ,j])
  qqline(rwmR[ ,j])
  qqnorm(rwmC[ ,j])
  qqline(rwmC[ ,j])
}

dir <- 'C:/'
source(paste0(dir, 'rwmR.R'))
sourceCpp(paste0(dir, 'rwmC.cpp'))
microbenchmark(rwmR(2, x0, N), rwmC(2, x0, N))

