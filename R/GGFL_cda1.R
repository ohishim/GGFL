#' @title Sub-function for GGFL
#' @description \code{GGFL1} This is required to execute "GGFL" or "pGGFL" (v0.2.0)
#'
#' @importFrom magrittr %>%
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

  v <- sapply(1:r, function(s){
    bs <- B[s,]
    out1 <- 2*(M%*%bs - c) %>% drop
    out2 <- lapply((1:r)[-s], function(l){
      a <- bs - B[l,]
      Lambda[l] * a / norm2(a)
    }) %>% Reduce("+", .)

    return(norm2(out1 + out2))
  })

  s <- which(Lambda >= v)

  if(length(s) == 1)
  {
    Beta.hat <- B[s,]
    type <- "join"

  } else if(length(s) == 0)
  {
    check <- apply(B, 1, function(b){all(b == Beta0)}) %>% which
    if(length(check) == 1)
    {
      v <- lapply(check, function(s){
        bs <- B[s,]
        out1 <- 2*(M%*%bs - c) %>% drop
        out2 <- lapply((1:r)[-s], function(l){
          a <- bs - B[l,]
          Lambda[l] * a / norm2(a)
        }) %>% Reduce("+", .)

        return(out1 + out2)
      })[[1]]

      Beta0 <- Beta0 - (v/norm2(v))

      type <- "disjoin"
    } else if(length(check) > 1)
    {
      stop("The error 1 occurs in GGFL1")
    }

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

    Beta.aft <- Beta0
    dif <- 1
    iter <- 0

    while((dif > tol) & (iter < maxit))
    {
      iter <- iter + 1

      Beta.bef <- Beta.aft
      Beta.aft <- solve(Q(Beta.bef)) %*% h(Beta.bef) %>% drop

      dif <- convC(Beta.aft, M, c, Lambda, B, r, Beta.bef)
    } #end while dif

    Beta.hat <- Beta.aft

  } else
  {
    stop("The error 2 occurs in GGFL1")
  }

  return(
    list(
      solution = Beta.hat,
      type = type,
      idx = s
    )
  )
}
