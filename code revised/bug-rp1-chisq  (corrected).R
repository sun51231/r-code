

rm(list=ls())
#p=90, chisquare from https://github.com/sun51231/r-code/blob/main/chisq90.R

#### 

loadings <- c(.7, .7, .7, .7, .7)
lambda <- cbind(
  c(loadings, rep(0, 85)),
  c(rep(0, 5), loadings, rep(0, 80)),
  c(rep(0, 10), loadings, rep(0, 75)),
  c(rep(0, 15), loadings, rep(0, 70)),
  c(rep(0, 20), loadings, rep(0, 65)),
  c(rep(0, 25), loadings, rep(0, 60)),
  c(rep(0, 30), loadings, rep(0, 55)),
  c(rep(0, 35), loadings, rep(0, 50)),
  c(rep(0, 40), loadings, rep(0, 45)),
  c(rep(0, 45), loadings, rep(0, 40)),
  c(rep(0, 50), loadings, rep(0, 35)),
  c(rep(0, 55), loadings, rep(0, 30)),
  c(rep(0, 60), loadings, rep(0, 25)),
  c(rep(0, 65), loadings, rep(0, 20)),
  c(rep(0, 70), loadings, rep(0, 15)),
  c(rep(0, 75), loadings, rep(0, 10)),
  c(rep(0, 80), loadings, rep(0, 5)),
  c(rep(0, 85), loadings)
)
phi <- matrix(c(1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,
                .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1
), ncol = 18)
psi <- 0.51*diag(90)
Sigma <- lambda %*% phi %*% t(lambda) + psi


model <- model90 <-  "
F1 =~ x1 + x2 + x3 + x4 + x5
F2 =~ x6 + x7 + x8 + x9 + x10
F3 =~ x11 + x12 + x13 + x14 + x15
F4 =~ x16 + x17 + x18 + x19 + x20
F5 =~ x21 + x22 + x23 + x24 + x25
F6 =~ x26 + x27 + x28 + x29 + x30
F7 =~ x31 + x32 + x33 + x34 + x35
F8 =~ x36 + x37 + x38 + x39 + x40
F9 =~ x41 + x42 + x43 + x44 + x45
F10 =~ x46 + x47 + x48 + x49 + x50
F11 =~ x51 + x52 + x53 + x54 + x55
F12 =~ x56 + x57 + x58 + x59 + x60
F13 =~ x61 + x62 + x63 + x64 + x65
F14 =~ x66 + x67 + x68 + x69 + x70
F15 =~ x71 + x72 + x73 + x74 + x75
F16 =~ x76 + x77 + x78 + x79 + x80
F17 =~ x81 + x82 + x83 + x84 + x85
F18 =~ x86 + x87 + x88 + x89 + x90
"

m<-18


Sigma.sqrt <- chol(Sigma)  # Cholesky

#get_pvalues <- function(n) {

n <- 10^4
  p <- ncol(Sigma.sqrt)
  k <- 2 * p
  A <- rbind(diag(p), diag(p)) / sqrt(2)    
  w <- (matrix(rchisq(n * k, df = 4), n, k) - 4) / sqrt(8)          
  r <- sqrt(6 / rchisq(n, df = 8))          
  epsilon <- (w %*% A) * r                         
  dat <- epsilon %*% t(Sigma.sqrt) 
  
  colnames(dat) <- paste0("x", 1:90)
  
  #check if model is misspecified
  
  residuals <- lavaan::lav_matrix_vech(cov(dat)- Sigma)# should be close to zero
  
  
  plot(residuals)#something is wrong the sample covariance matrix does not converge to Sigma
  

  
  
  
  
  
###
###

###  The corrected code is as follows:

  
library(Matrix)  
library(expm)

loadings <- c(.7, .7, .7, .7, .7)
lambda <- cbind(
    c(loadings, rep(0, 85)),
    c(rep(0, 5), loadings, rep(0, 80)),
    c(rep(0, 10), loadings, rep(0, 75)),
    c(rep(0, 15), loadings, rep(0, 70)),
    c(rep(0, 20), loadings, rep(0, 65)),
    c(rep(0, 25), loadings, rep(0, 60)),
    c(rep(0, 30), loadings, rep(0, 55)),
    c(rep(0, 35), loadings, rep(0, 50)),
    c(rep(0, 40), loadings, rep(0, 45)),
    c(rep(0, 45), loadings, rep(0, 40)),
    c(rep(0, 50), loadings, rep(0, 35)),
    c(rep(0, 55), loadings, rep(0, 30)),
    c(rep(0, 60), loadings, rep(0, 25)),
    c(rep(0, 65), loadings, rep(0, 20)),
    c(rep(0, 70), loadings, rep(0, 15)),
    c(rep(0, 75), loadings, rep(0, 10)),
    c(rep(0, 80), loadings, rep(0, 5)),
    c(rep(0, 85), loadings)
)
phi <- matrix(c(1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1,.3,
                  .3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,.3,1
), ncol = 18)
psi <- 0.51*diag(90)
Sigma <- lambda %*% phi %*% t(lambda) + psi
  
  
model <- model90 <-  "
F1 =~ x1 + x2 + x3 + x4 + x5
F2 =~ x6 + x7 + x8 + x9 + x10
F3 =~ x11 + x12 + x13 + x14 + x15
F4 =~ x16 + x17 + x18 + x19 + x20
F5 =~ x21 + x22 + x23 + x24 + x25
F6 =~ x26 + x27 + x28 + x29 + x30
F7 =~ x31 + x32 + x33 + x34 + x35
F8 =~ x36 + x37 + x38 + x39 + x40
F9 =~ x41 + x42 + x43 + x44 + x45
F10 =~ x46 + x47 + x48 + x49 + x50
F11 =~ x51 + x52 + x53 + x54 + x55
F12 =~ x56 + x57 + x58 + x59 + x60
F13 =~ x61 + x62 + x63 + x64 + x65
F14 =~ x66 + x67 + x68 + x69 + x70
F15 =~ x71 + x72 + x73 + x74 + x75
F16 =~ x76 + x77 + x78 + x79 + x80
F17 =~ x81 + x82 + x83 + x84 + x85
F18 =~ x86 + x87 + x88 + x89 + x90
"
  
m<-18




Sigma.sqrt <- sqrtm(Sigma)       #  Revised   ## Square Root,  Non-Cholesky

####  #  Revised 
generate_A0 <- function(p) {
  ones_p <- matrix(1, nrow = p, ncol = 1)
  identity_p <- diag(p)
  matrix_sum <- identity_p + ones_p %*% t(ones_p)
  inv_sqrt <- solve(sqrtm(matrix_sum))
  A0_top <- identity_p
  A0_bottom <- t(ones_p)
  A0 <- rbind(A0_top, A0_bottom) %*% inv_sqrt
  return(A0)
}

###      Chi-square distribution

n <- 10^4

p <- ncol(Sigma.sqrt)
w_raw <- matrix(rchisq(n * (p+1), df = 4), nrow = n, ncol = p+1)   #####  Revised  
w <- (w_raw - 4) / (2 * sqrt(2))
chi_sq_r <- rchisq(n, df = 8)
r <- sqrt(6 / chi_sq_r)
A0 <- generate_A0(p)
epsilon <- r * (w %*% A0)
dat <- epsilon %*% t(Sigma.sqrt) #
colnames(dat) <- paste0("x", 1:p)



residuals <- lavaan::lav_matrix_vech(cov(dat)- Sigma)#  to zero
plot(residuals)




###    t-distribution


Sigma.sqrt <- sqrtm(Sigma)       #  Revised   ## Square Root,  Non-Cholesky


n <- 10^4
p <- ncol(Sigma.sqrt)
w <- matrix(rnorm(n * p), nrow = n)       # n x p 
r <- sqrt(6 / rchisq(n, df = 8))  
A <- diag(p)
epsilon <- r * (w %*% A)                         
dat <- epsilon %*% t(Sigma.sqrt)   ###  

colnames(dat) <- paste0("x", 1:p)





residuals <- lavaan::lav_matrix_vech(cov(dat)- Sigma)#  to zero
plot(residuals)



###  


