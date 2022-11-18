## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(exDE)
library(MASS)
library(expm)
library(deSolve)
library(data.table)
library(ggplot2)

## -----------------------------------------------------------------------------
pars <- new.env()
pars$nPatches <- 3
pars$nStrata <- 3
pars$nHabitats <- 3

# human parameters
b <- 0.55
c <- 0.15
r <- 1/200
wf <- rep(1, pars$nStrata)

pfpr <- runif(n = pars$nStrata, min = 0.25, max = 0.35)
H <- rpois(n = pars$nStrata, lambda = 1000)
X <- rbinom(n = pars$nStrata, size = H, prob = pfpr)

Psi <- matrix(
  data = c(
    0.9, 0.05, 0.05,
    0.05, 0.9, 0.05,
    0.05, 0.05, 0.9
  ), nrow = pars$nStrata, ncol = pars$nPatches, byrow = T
)
Psi <- t(Psi)

# adult mosquito parameters
f <- 0.3
q <- 0.9
g <- 1/10
sigma <- 1/100
nu <- 1/2
eggsPerBatch <- 30
tau <- 11

# mosquito movement calK
calK <- matrix(0, pars$nPatches, pars$nPatches)
calK[upper.tri(calK)] <- 1/(pars$nPatches-1)
calK[lower.tri(calK)] <- 1/(pars$nPatches-1)
calK <- calK/rowSums(calK)
calK <- t(calK)

Omega <- make_Omega(g = g, sigma = sigma, K = calK, nPatches = pars$nPatches)
Omega_inv <- solve(Omega)
Upsilon <- expm::expm(-Omega * tau)
Upsilon_inv <- expm::expm(Omega * tau)

## -----------------------------------------------------------------------------
# derived EIR to sustain equilibrium pfpr
EIR <- diag(1/b, pars$nStrata) %*% ((r*X) / (H - X))

# ambient pop
W <- Psi %*% H

# biting distribution matrix
beta <- diag(wf) %*% t(Psi) %*% diag(1/as.vector(W), pars$nPatches)

# kappa
kappa <- t(beta) %*% (X*c)

# equilibrium solutions
Z <- diag(1/(f*q), pars$nPatches) %*% ginv(beta) %*% EIR
MY <- diag(1/as.vector(f*q*kappa), pars$nPatches) %*% Upsilon_inv %*% Omega %*% Z
Y <- Omega_inv %*% (diag(as.vector(f*q*kappa), pars$nPatches) %*% MY)
M <- MY + Y
G <- solve(diag(nu+f, pars$nPatches) + Omega) %*% diag(f, pars$nPatches) %*% M
Lambda <- Omega %*% M

## -----------------------------------------------------------------------------
# aquatic habitat membership matrix (assume 1 habitat per patch)
calN <- matrix(0, pars$nPatches, pars$nHabitats)
diag(calN) <- 1

# egg dispersal matrix (assume 1 habitat per patch)
calU <- matrix(0, pars$nHabitats, pars$nPatches)
diag(calU) <- 1

alpha <- as.vector(solve(calN) %*% Lambda)

psi <- 1/10
phi <- 1/12
eta <- as.vector(calU %*% G * nu * eggsPerBatch)

L <- alpha/psi
theta <- (eta - psi*L - phi*L)/(L^2)

## -----------------------------------------------------------------------------
pars$calU <- calU
pars$calN <- calN

# ITN coverage
ITN_cov <- function(t) {(sin(2*pi*t / 365) + 1) / 2}

# set parameters
make_parameters_exogenous_null(pars)
make_parameters_vc_lemenach(pars, phi = ITN_cov)
make_parameters_MYZ_RM_dde(pars = pars, g = g, sigma = sigma, calK = calK, f = f, q = q, nu = nu, eggsPerBatch = eggsPerBatch, tau = tau, M0 = M, G0 = G, Y0 = Y, Z0 = Z)
make_parameters_L_basic(pars = pars, psi = psi, phi = phi, theta = theta, L0 = L)
make_parameters_X_SIS(pars = pars, b = b, c = c, r = r, Psi = Psi, wf = wf, X0 = X, H = H)
make_indices(pars = pars)

## -----------------------------------------------------------------------------
# ICs
y <- rep(NaN, max(pars$X_ix))
y[pars$L_ix] <- pars$Lpar$L0
y[pars$X_ix] <- pars$Xpar$X0
y[pars$M_ix] <- pars$MYZpar$M0
y[pars$G_ix] <- pars$MYZpar$G0
y[pars$Y_ix] <- pars$MYZpar$Y0
y[pars$Z_ix] <- pars$MYZpar$Z0
y[pars$Upsilon_ix] <- as.vector(Upsilon)

# solve the model
out <- dede(y = y, times = 0:(365*5), func = xDE_diffeqn, parms = pars)

## ---- out.width = "100%"------------------------------------------------------
colnames(out)[pars$L_ix+1] <- paste0('L_', 1:pars$nHabitats)
colnames(out)[pars$M_ix+1] <- paste0('M_', 1:pars$nPatches)
colnames(out)[pars$G_ix+1] <- paste0('G_', 1:pars$nPatches)
colnames(out)[pars$Y_ix+1] <- paste0('Y_', 1:pars$nPatches)
colnames(out)[pars$Z_ix+1] <- paste0('Z_', 1:pars$nPatches)
colnames(out)[pars$X_ix+1] <- paste0('X_', 1:pars$nStrata)
out <- out[, -c(pars$Upsilon_ix+1)]

out <- as.data.table(out)
out <- melt(out, id.vars = 'time')
out[, c("Component", "Stratification") := tstrsplit(variable, '_', fixed = TRUE)]
out[, variable := NULL]

ggplot(data = out, mapping = aes(x = time, y = value, color = Stratification)) +
  geom_line() +
  facet_wrap(. ~ Component, scales = 'free') +
  theme_bw()

