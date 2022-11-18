## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(exDE)
library(expm)
library(deSolve)
library(data.table)
library(ggplot2)

## -----------------------------------------------------------------------------
nPatches <- 3
f <- 0.3
q <- 0.9
g <- 1/20
sigma <- 1/10
tau <- 11
nu <- 1/2
eggsPerBatch <- 30

calK <- matrix(0, nPatches, nPatches)
calK[1, 2:3] <- c(0.2, 0.8)
calK[2, c(1,3)] <- c(0.5, 0.5)
calK[3, 1:2] <- c(0.7, 0.3)
calK <- t(calK)

Omega <- make_Omega(g, sigma, calK, nPatches)
Upsilon <- expm::expm(-Omega * tau)

kappa <- c(0.1, 0.075, 0.025)
Lambda <- c(5, 10, 8)

params <- list(
  nPatches = nPatches
)
params <- list2env(params)

## -----------------------------------------------------------------------------
Omega_inv <- solve(Omega)

M_eq <- as.vector(Omega_inv %*% Lambda)
G_eq <- as.vector(solve(diag(nu+f, nPatches) + Omega) %*% diag(f, nPatches) %*% M_eq)
Y_eq <- as.vector(solve(diag(f*q*kappa) + Omega) %*% diag(f*q*kappa) %*% M_eq)
Z_eq <- as.vector(Omega_inv %*% Upsilon %*% diag(f*q*kappa) %*% (M_eq - Y_eq))

## ---- out.width = "100%"------------------------------------------------------
 make_parameters_MYZ_RM_ode(pars = params, g = g, sigma = sigma, calK = calK, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = rep(0, nPatches), G0 = rep(0, nPatches), Y0 = rep(0, nPatches), Z0 = rep(0, nPatches))
  make_indices(params)

y0 <- rep(NaN, 12)
y0[params$M_ix] <- M_eq
y0[params$G_ix] <- G_eq
y0[params$Y_ix] <- Y_eq
y0[params$Z_ix] <- Z_eq
y0[params$Upsilon_ix] <- as.vector(Upsilon)

MosyBehavior <- MosquitoBehavior(0, y, params)

out <- deSolve::ode(y = y0, times = 0:50, func = function(t, y, pars, Lambda, kappa, MosyBehavior) {
  list(dMYZdt(t, y, pars, Lambda, kappa, MosyBehavior))
}, parms = params, method = 'lsoda', Lambda = Lambda, kappa = kappa, MosyBehavior = MosyBehavior)

colnames(out)[params$M_ix+1] <- paste0('M_', 1:params$nPatches)
colnames(out)[params$G_ix+1] <- paste0('G_', 1:params$nPatches)
colnames(out)[params$Y_ix+1] <- paste0('Y_', 1:params$nPatches)
colnames(out)[params$Z_ix+1] <- paste0('Z_', 1:params$nPatches)
out <- out[, -c(params$Upsilon_ix+1)]

out <- as.data.table(out)
out <- melt(out, id.vars = 'time')
out[, c("Component", "Patch") := tstrsplit(variable, '_', fixed = TRUE)]
out[, variable := NULL]

ggplot(data = out, mapping = aes(x = time, y = value, color = Patch)) +
  geom_line() +
  facet_wrap(. ~ Component, scales = 'free') +
  theme_bw()

