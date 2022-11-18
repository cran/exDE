## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(exDE)
library(expm)
library(data.table)
library(ggplot2)

## -----------------------------------------------------------------------------
params <- list(
  nHabitats = 5,
  nPatches = 3,
  nStrata = 4
)
params <- list2env(params)

## -----------------------------------------------------------------------------
calN <- matrix(
  data = c(1,1,1,0,0,
           0,0,0,1,1,
           0,0,0,0,0), 
  nrow = params$nPatches, ncol = params$nHabitats, byrow = TRUE
)
params$calN <- calN

## -----------------------------------------------------------------------------
xi <- matrix(c(.7, .2, .1, .8, .2), 5, 1) 
params$calU <- t(calN %*% diag(as.vector(xi)))

## -----------------------------------------------------------------------------
L0 <- rep(1, params$nHabitats)
psi <- rep(1/8, params$nHabitats)
phi <- rep(1/8, params$nHabitats)
theta <- c(1/10, 1/20, 1/40, 1/100, 1/10) 

make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta, L0 = L0)

## -----------------------------------------------------------------------------
g <- rep(1/12, params$nPatches) 
sigma <- rep(1/12/2, params$nPatches) 
calK <- t(matrix(
  c(c(0, .6, .3), 
    c(.4, 0, .7), 
    c(.6, .4, 0)), 3, 3))
f <- rep(1/3, params$nPatches) 
q <- rep(0.9, params$nPatches) 
nu <- c(1/3,1/3,0)   
eggsPerBatch <- 30 
tau <- 12 
M0 <- rep(100, params$nPatches)
G0 <- rep(10, params$nPatches)
Y0 <- rep(1, params$nPatches)
Z0 <- rep(0, params$nPatches)

Omega <- make_Omega(g, sigma, calK, params$nPatches)
Upsilon <- expm::expm(-Omega*tau)

make_parameters_MYZ_RM_ode(pars = params, g = g, sigma = sigma, calK = calK, tau = tau, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, M0 = M0, G0 = G0, Y0 = Y0, Z0 = Z0)

## -----------------------------------------------------------------------------
calJ <- t(matrix(
  c(c(0,0,0,0),
    c(1,1,0,0),
    c(0,0,1,1)), 4, 3
))

## -----------------------------------------------------------------------------
Psi <- t(matrix(
  c(c(0.01,0.01,0.001,0.001),
    c(.95,.92,.04,.02),
    c(.04,.02,.959,.929)), 4, 3
))

## -----------------------------------------------------------------------------
H <- matrix(c(10,90, 100, 900), 4, 1)
X0 <- as.vector(0.2*H)
r <- rep(1/200, params$nStrata)
b <- rep(0.55, params$nStrata)
c <- c(0.1, .02, .1, .02)

make_parameters_X_SIS(pars = params, b = b, c = c, r = r, Psi = Psi, X0 = X0, H = as.vector(H))

## -----------------------------------------------------------------------------
make_parameters_exogenous_null(params)
make_parameters_vc_null(params)

## -----------------------------------------------------------------------------
make_indices(params)

## -----------------------------------------------------------------------------
y <- rep(NaN, max(params$X_ix))
y[params$L_ix] <- as.vector(L0)
y[params$M_ix] <- as.vector(M0)
y[params$G_ix] <- as.vector(G0)
y[params$Y_ix] <- as.vector(Y0)
y[params$Z_ix] <- as.vector(Z0)
y[params$Upsilon_ix] <- as.vector(Upsilon)
y[params$X_ix] <- as.vector(X0)

## -----------------------------------------------------------------------------
out <- deSolve::ode(y = y, times = 0:365, func = xDE_diffeqn, parms = params, method = "lsoda")

## ---- out.width = "100%"------------------------------------------------------
colnames(out)[params$L_ix+1] <- paste0('L_', 1:params$nHabitats)
colnames(out)[params$M_ix+1] <- paste0('M_', 1:params$nPatches)
colnames(out)[params$G_ix+1] <- paste0('G_', 1:params$nPatches)
colnames(out)[params$Y_ix+1] <- paste0('Y_', 1:params$nPatches)
colnames(out)[params$Z_ix+1] <- paste0('Z_', 1:params$nPatches)
colnames(out)[params$X_ix+1] <- paste0('X_', 1:params$nStrata)
out <- out[, -c(params$Upsilon_ix+1)]

out <- as.data.table(out)
out <- melt(out, id.vars = 'time')
out[, c("Component", "Stratification") := tstrsplit(variable, '_', fixed = TRUE)]
out[, variable := NULL]

ggplot(data = out, mapping = aes(x = time, y = value, color = Stratification)) +
  geom_line() +
  facet_wrap(. ~ Component, scales = 'free') +
  theme_bw()

