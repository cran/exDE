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
nStrata <- 3
H <- c(100, 500, 250)
X <- c(20, 120, 80)
b <- 0.55
c <- 0.15
r <- 1/200
Psi <- matrix(data = 1,nrow = 1, ncol = nStrata)

EIR <- diag(1/b, nStrata) %*% ((r*X)/(H-X))

params <- list(
  nStrata = nStrata
)
params <- list2env(params)

make_parameters_X_SIS(pars = params, b = b, c = c, r = r, Psi = Psi, X0 = X, H = H)
make_indices(params)

y0 <- rep(0, 3)
y0[params$X_ix] <- X

out <- deSolve::ode(y = y0, times = c(0, 365), func = function(t, y, pars, EIR) {
  list(dXdt(t, y, pars, EIR))
}, parms = params, method = 'lsoda', EIR = as.vector(EIR))

colnames(out)[params$X_ix+1] <- paste0('X_', 1:params$nStrata)

out <- as.data.table(out)
out <- melt(out, id.vars = 'time')
out[, c("Component", "Strata") := tstrsplit(variable, '_', fixed = TRUE)]
out[, variable := NULL]

ggplot(data = out, mapping = aes(x = time, y = value, color = Strata)) +
  geom_line() +
  theme_bw()

