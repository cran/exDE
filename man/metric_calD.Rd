% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{metric_calD}
\alias{metric_calD}
\title{Parasite dispersal by humans}
\usage{
metric_calD(W, beta, b, D, H)
}
\arguments{
\item{W}{ambient human population at each patch}

\item{beta}{the biting distribution matrix}

\item{b}{transmission efficiency from mosquitoes to humans}

\item{D}{human transmitting capacity}

\item{H}{human population size of each strata}
}
\value{
a numeric \link{matrix}
}
\description{
Compute the \code{p} by \code{p} matrix \eqn{\mathcal{D}} whose columns describe
how potentially infectious person time from persons in that patch are dispersed
across other patches.
\deqn{\mathcal{D} = \mbox{diag}\left(W \right) \cdot \beta^T \cdot \mbox{diag}\left(bDH\right) \cdot \beta}
}
