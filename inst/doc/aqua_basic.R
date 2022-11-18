## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(exDE)
library(deSolve)
library(data.table)
library(ggplot2)

## ---- out.width = "100%"------------------------------------------------------
nHabitats <- 3
alpha <- c(10, 50, 20)
eta <- c(250, 500, 170)
psi <- 1/10
phi <- 1/12

L <- alpha/psi
theta <- (eta - psi*L - phi*L)/(L^2)

params <- list(
  nHabitats = nHabitats
)
params <- list2env(params)

make_parameters_L_basic(pars = params, psi = psi, phi = phi, theta = theta, L0 = L)
make_indices(params)

y0 <- L

out <- deSolve::ode(y = y0, times = 0:50, func = function(t, y, pars, eta) {
  list(dLdt(t, y, pars, eta))
}, parms = params, method = 'lsoda', eta = eta)

colnames(out)[params$L_ix+1] <- paste0('L_', 1:params$nHabitats)

out <- as.data.table(out)
out <- melt(out, id.vars = 'time')
out[, c("Component", "Patch") := tstrsplit(variable, '_', fixed = TRUE)]
out[, variable := NULL]

ggplot(data = out, mapping = aes(x = time, y = value, color = Patch)) +
  geom_line() +
  theme_bw()

