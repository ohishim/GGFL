#' @title Sub-function for GGFL
#' @description \code{GGFL1} This is required to execute "GGFL" or "pGGFL" (v0.4.1)
#'
#' @importFrom magrittr %>% equals
#'
#' @param M k * k positive definite matrix
#' @param c k-dimensional vector
#' @param Lambda r-dimensional vector which all elements are positive
#' @param B r * k matrix
#' @param Beta0 k-dimensional initial vector
#' @param tol tolerance for convergence
#' @param convC function which calculates a convergence criterion
#' @param maxit iteration limit
#'
#' @return a list object which has the following elements:
#' \item{solution}{vector of the minimizer}
#'
#' \item{type}{`"join"` or `"disjoin"`}
#'
#' \item{idx}{index of parameter when `type=join`}
#'
#' @examples
#' #GGFL1(M, c, Lambda, B, Beta0)

GGFL1 <- function(M, c, Lambda, B, Beta0, tol=1e-5, convC=NULL, maxit=500){

  type <- "-"

  norm2 <- function(x){sqrt(sum(x^2))}

  if(is.null(convC)){
    convC <- function(Beta, M, c, Lambda, B, r, Beta0){
      max((Beta - Beta0)^2 / Beta0^2)
    }
  }

  r <- length(Lambda)

  #---   checking fusion   -----------------------------------------------------

  V <- lapply(1:r, function(s){
    bs <- B[s,]
    out1 <- 2*(M%*%bs - c) %>% drop

    if(r == 1)
    {
      out2 <- 0
    } else
    {
      out2 <- lapply((1:r)[-s], function(l){
        a <- bs - B[l,]
        Lambda[l] * a / norm2(a)
      }) %>% Reduce("+", .)
    }

    return(out1 + out2)
  })

  s <- which(Lambda >= sapply(V, norm2))

  if(length(s) == 1)
  { #---   fusion case   -------------------------------------------------------
    Beta.hat <- B[s,]
    type <- "join"

  } else if(length(s) == 0)
  { #---   searching solution   ------------------------------------------------
    k <- length(Beta0)

    Q <- function(Beta){
      a <- sapply(1:r, function(l){Lambda[l] / norm2(Beta - B[l,])}) %>% sum
      return(2*M + diag(a, k, k))
    }

    h <- function(Beta){
      a <- lapply(1:r, function(l){
        Lambda[l]*B[l,] / norm2(Beta - B[l,])
      }) %>% Reduce("+", .)
      return(2*c + a)
    }

    updateB <- function(Beta){
      check <- scale(B, center=Beta, scale=F) %>% abs %>% rowSums %>% equals(0) %>% which
      checkN <- length(check)

      if(checkN == 0)
      {
        .Beta <- Beta
      } else if(checkN == 1)
      {
        v <- V[[check]]

        if(r == 1)
        {
          d <- 1
        } else
        {
          d <- Lambda[check]*min(
            scale(B[-check,,drop=F], center=B[check,], scale=F) %>% apply(1, norm2)
          )/sum(Lambda)
        }

        .Beta <- Beta - (d/5)*v/norm2(v)
      } else
      {
        stop("The error 1 occurs in `GGFL1`")
      }

      return(chol_solve(Q(.Beta), h(.Beta)) %>% drop)
    }

    Beta.aft <- Beta0
    dif <- 1
    iter <- 0

    while((dif > tol) & (iter < maxit))
    {
      iter <- iter + 1

      Beta.bef <- Beta.aft
      Beta.aft <- updateB(Beta.bef)

      dif <- convC(Beta.aft, M, c, Lambda, B, r, Beta.bef)
    } #end while dif

    obj <- function(Beta){
      sum(Beta * drop(M%*%Beta)) - 2*sum(c*Beta) + sum(
        sapply(1:r, function(j){
          Lambda[j]*norm2(Beta - B[j,])
        })
      )
    }

    idx <- apply(B, 1, obj) %>% c(obj(Beta.aft)) %>% which.min

    if(idx == r+1)
    {
      Beta.hat <- Beta.aft

      if(scale(B, center=Beta0, scale=F) %>% abs %>% rowSums %>% equals(0) %>% any){
        type <- "disjoin"
      }
    } else
    {
      Beta.hat <- B[idx,]
      type <- "join"
      s <- idx
    }

  } else
  {
    stop("The error 2 occurs in `GGFL1`")
  }

  return(
    list(
      solution = Beta.hat,
      type = type,
      idx = s
    )
  )
}
