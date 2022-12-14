% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adult-RM.R
\name{F_eggs.RM}
\alias{F_eggs.RM}
\title{Number of eggs laid by adult mosquitoes}
\usage{
\method{F_eggs}{RM}(t, y, pars)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \link{environment}}
}
\value{
a \link{numeric} vector of length \code{nPatches}
}
\description{
Implements \link{F_eggs} for the generalized RM model.
}
