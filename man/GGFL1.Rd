% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GGFL_cda1.R
\name{GGFL1}
\alias{GGFL1}
\title{Sub-function for GGFL}
\usage{
GGFL1(M, c, Lambda, B, Beta0, tol = 1e-05, convC = NULL, maxit = 500)
}
\arguments{
\item{M}{k * k positive definite matrix}

\item{c}{k-dimensional vector}

\item{Lambda}{r-dimensional vector which all elements are positive}

\item{B}{r * k matrix}

\item{Beta0}{k-dimensional initial vector}

\item{tol}{tolerance for convergence}

\item{convC}{function which calculates a convergence criterion}

\item{maxit}{iteration limit}
}
\value{
a list object which has the following elements:
\item{solution}{vector of the minimizer}

\item{type}{\code{"join"} or \code{"disjoin"}}

\item{idx}{index of parameter when \code{type=join}}
}
\description{
\code{GGFL1} This is required to execute "GGFL" or "pGGFL" (v0.4.1)
}
\examples{
#GGFL1(M, c, Lambda, B, Beta0)
}
