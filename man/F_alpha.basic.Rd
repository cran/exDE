% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aquatic-basic.R
\name{F_alpha.basic}
\alias{F_alpha.basic}
\title{Number of newly emerging adults from each larval habitat}
\usage{
\method{F_alpha}{basic}(t, y, pars)
}
\arguments{
\item{t}{current simulation time}

\item{y}{state vector}

\item{pars}{an \link{environment}}
}
\value{
a \link{numeric} vector of length \code{nHabitats}
}
\description{
Implements \link{F_alpha} for the basic competition model.
}
