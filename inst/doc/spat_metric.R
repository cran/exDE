## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(exDE)

## -----------------------------------------------------------------------------
l <- 5
p <- 3
n <- 4

calN <- matrix(
  data = c(1,1,1,0,0,
           0,0,0,1,1,
           0,0,0,0,0), 
  nrow = p, ncol = l, byrow = TRUE
)


xi <- matrix(c(.7, .2, .1, .8, .2), 5, 1) 
calU <- t(calN %*% diag(as.vector(xi)))

g <- rep(1/12, p) 
sigma <- rep(1/12/2, p) 
calK <- t(matrix(
  c(c(0, .6, .3), 
    c(.4, 0, .7), 
    c(.6, .4, 0)), 3, 3))
f <- rep(1/3, p) 
q <- rep(0.9, p) 
tau <- 12 
M <- rep(100, p)

Omega <- diag(g) + (diag(1,p)-calK) %*% diag(sigma)
Lambda <- Omega %*% M

H <- matrix(c(10,90, 100, 900), 4, 1)
X <- as.vector(0.2*H)
r <- rep(1/200, n)
b <- rep(0.55, n)
c <- c(0.1, .02, .1, .02)

calJ <- t(matrix(
  c(c(0,0,0,0),
    c(1,1,0,0),
    c(0,0,1,1)), n, p
))

Theta <- t(matrix(
  c(c(0.01,0.01,0.001,0.001),
    c(.95,.92,.04,.02),
    c(.04,.02,.959,.929)), n, p
))

W <- Theta %*% H
beta <- t(Theta) %*% diag_inverse(W)

D <- c/r

calV <- metric_calV(f = f, q = q, Omega = Omega, tau = tau, M = M, W = W)
calD <- metric_calD(W = W, beta = beta, b = b, D = D, H = H)
calR <- metric_calR(b = b, beta = beta, calV = calV, W = W, D = D, H = H)
calZ <- metric_calZ(Omega = Omega, tau = tau, f = f, q = q, M = M, W = W, calD = calD)

