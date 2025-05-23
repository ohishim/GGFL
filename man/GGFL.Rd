% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GGFL_cda.R
\name{GGFL}
\alias{GGFL}
\title{Coordinate optimization for GGFL}
\usage{
GGFL(
  y,
  X,
  group,
  adj,
  Lambda = "default",
  lambda.type = NULL,
  adaptive = TRUE,
  weight = NULL,
  alpha = NULL,
  standardize = TRUE,
  MPinv = FALSE,
  tol = 1e-05,
  maxit = 500,
  progress = FALSE,
  out.all = FALSE,
  cwu.control = NULL
)
}
\arguments{
\item{y}{vector of an objective variable}

\item{X}{matrix of explanatory variables without intercept}

\item{group}{numeric vector of group labels for each individual}

\item{adj}{data.frame or matrix with two columns of adjacent information}

\item{Lambda}{candidates of tuning parameter
if \code{"default"}, the candidates are defined following the paper;
if \code{"uniform"}, the candidates are defined by uniformly dividing;
if vector/scalar, the candidates are defined by following \code{lambda.type}}

\item{lambda.type}{option when \code{Lambda} is numeric;
if "value", \code{Lambda} is searching points;
if "rate", \code{lam.max*Lambda} is searching points}

\item{adaptive}{option only when \code{weight=NULL}; if \code{TRUE}, adaptive penalty is adopted}

\item{weight}{list of penalty weight (default is \code{NULL})}

\item{alpha}{a value expressing penalty strength for EGCV criterion;
default is \code{NULL} which corresponds to \code{log(n)}}

\item{standardize}{if \code{TRUE}, \code{y} and \code{X} are standardized in the sense of norm}

\item{MPinv}{if \code{TRUE}, the ordinary least squares estimator is calculated by the Moore-Penrose inverse matrix}

\item{tol}{tolerance for convergence}

\item{maxit}{iteration limit}

\item{progress}{If \code{TRUE}, progress is displayed of the form
"lambda index_iteration number_convergence criterion_objective function"}

\item{out.all}{if \code{TRUE}, results for all tuning parameters are output;
if \code{FALSE}, results for only the optimal tuning parameter are output}

\item{cwu.control}{a list object of parameters for controlling the coordinate wise update
which can has the two elements \code{conv.check} and \code{tol};
\code{conv.check} is one of "coef" (default) and "grad" and solution convergence is
determined based on the coefficients or gradient, respectively;
\code{tol} is a tolerance for convergence of which the default is 1e-5}
}
\value{
a list object which has the following elements:
\item{coefficients}{a list object which has
\code{GGFL}: matrix of GGFL estimates; \code{OLS}: matrix of OLS estimates;
\code{MAX} matrix of estimates when all groups are equal
}

\item{fitted.values}{a matrix of fitted.values corresponding to
each element of \code{coefficients}
}

\item{summary}{a data.frame with results of the optimal solution}

\item{weight}{a list of penalty weight}

\item{cluster}{a list of cluster}

\item{cluster.labels}{a vector expressing which cluster each group is in}

\item{cwu.control}{\code{cwu.control}}

\item{all.coef}{a list with all candidates of GGFL estimates}

\item{all.cluster}{a matrix of \code{cluster}s for all candidates}

\item{log}{results for all candidates}
}
\description{
\code{GGFL} solving GGFL optimization problem via coordinate descent algorithm (v1.0.2)
}
\examples{
#GGFL(y, X, group, adj)
}
