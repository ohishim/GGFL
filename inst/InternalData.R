
svd_solve <- function(M, b, tol=sqrt(.Machine$double.eps)){

  if(missing(b)){b <- diag(ncol(M))}

  svdM <- svd(M); d <- svdM$d; r <- sum(d > tol)
  l <- c(1/(d[1:r]), numeric(ncol(M) - r))

  return(drop(
    svdM$v %*% (l*crossprod(svdM$u, b))
  ))
}

chol_solve <- function(M, b){

  if(missing(b)){b <- diag(ncol(M))}

  R <- chol(M)
  y <- backsolve(R, b, transpose=TRUE)
  return(backsolve(R, y))
}

usethis::use_data(svd_solve, chol_solve, internal=T, overwrite = T)
