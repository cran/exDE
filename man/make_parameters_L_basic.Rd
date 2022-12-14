% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aquatic-basic.R
\name{make_parameters_L_basic}
\alias{make_parameters_L_basic}
\title{Make parameters for basic competition aquatic mosquito model}
\usage{
make_parameters_L_basic(pars, psi, phi, theta, L0)
}
\arguments{
\item{pars}{an \link{environment}}

\item{psi}{maturation rates for each aquatic habitat}

\item{phi}{density-independent mortality rates for each aquatic habitat}

\item{theta}{density-dependent mortality terms for each aquatic habitat}

\item{L0}{initial conditions}
}
\value{
a \link{list} with class \code{basic}.
}
\description{
Make parameters for basic competition aquatic mosquito model
}
