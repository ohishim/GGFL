#' @title Coordinate optimization for GGFL
#' @description \code{pGGFL} More general version of Coordinate optimization for GGFL.
#'   Only partial explanatory variables are penalized. (v1.0.0)
#'
#' @importFrom magrittr %>% add set_rownames
#' @importFrom MASS ginv
#' @importFrom purrr transpose
#'
#' @param y vector of an objective variable
#' @param X matrix of explanatory variables without intercept for local estimation
#'   (for GGFL penalty)
#' @param Z matrix of explanatory variables without intercept for global estimation
#' @param group numeric vector of group labels for each individual
#' @param adj data.frame or matrix with two columns of adjacent information
#' @param tol tolerance for convergence
#' @param Lambda candidates of tuning parameter
#'   if `"default"`, the candidates are defined following the paper;
#'   if `"uniform"`, the candidates are defined by uniformly dividing
#' @param alpha a value expressing penalty strength for EGCV criterion;
#'   default is `NULL` which corresponds to `log(n)`
#' @param standardize if `TRUE`, `X` is standardized in the sense of norm
#' @param MPinv if TRUE, the ordinary least squares estimator is calculated by the Moore-Penrose inverse matrix
#' @param maxit iteration limit
#' @param progress If `TRUE`, progress is displayed of the form
#'   "lambda index_iteration number_convergence criterion_objective function"
#' @param out.all if `TRUE`, results for all tuning parameters are output;
#'   if `FALSE`, results for only the optimal tuning parameter are output
#' @param cwu.control a list object of parameters for controlling the coordinate wise update
#'   which can has the two elements `conv.check` and `tol`;
#'   `conv.check` is one of "coef" (default) and "grad" and solution convergence is
#'     determined based on the coefficients or gradient, respectively;
#'   `tol` is a tolerance for convergence of which the default is 1e-5
#'
#' @return a list object which has the following elements:
#' \item{coefficients.local}{a list object for `X` which has
#'   `GGFL`: matrix of GGFL estimates; `OLS`: matrix of OLS estimates;
#'   `MAX` matrix of estimates when all groups are equal
#' }
#'
#' \item{coefficients.global}{a list object for `Z` which has
#'   `GGFL`: matrix of GGFL estimates; `OLS`: matrix of OLS estimates;
#'   `MAX` matrix of estimates when all groups are equal
#' }
#'
#' \item{fitted.values}{a matrix of fitted.values corresponding to
#'   each element of `coefficients`
#' }
#'
#' \item{summary}{a data.frame with results of the optimal solution}
#'
#' \item{weight}{a list of penalty weight}
#'
#' \item{cluster}{a vector expressing which cluster each group is in}
#'
#' \item{cwu.control}{`cwu.control`}
#'
#' \item{all.coef}{a list with all candidates of GGFL estimates}
#'
#' \item{all.cluster}{a matrix of `cluster`s for all candidates}
#'
#' \item{log}{results for all candidates}
#'
#' @export
#' @examples
#' #pGGFL(y, X, Z, group, adj)

pGGFL <- function(
  y, X, Z, group, adj, tol=1e-5, Lambda="default", alpha=NULL, standardize=TRUE,
  MPinv=FALSE, maxit=500, progress=FALSE, out.all=FALSE, cwu.control=NULL
){

  ##############################################################################
  ###   preparation
  ##############################################################################

  t1 <- proc.time()[3]

  norm2 <- function(x){sqrt(sum(x^2))}

  if(standardize)
  {
    normX <- c(1, colSums(X^2) %>% sqrt)
    X <- scale(cbind(1, X), center=F, scale=normX) %>% `attr<-`("scaled:scale", NULL)
    X_ <- as.data.frame(X) %>% split(group) %>% lapply(data.matrix)
  } else
  {
    X_ <- cbind(1, X) %>% as.data.frame %>% split(group) %>% lapply(data.matrix)
  }

  y_ <- split(y, group)
  Z_ <- as.data.frame(Z) %>% split(group) %>% lapply(data.matrix)
  D <- split(adj[,2], adj[,1])

  .y <- unlist(y_); .X <- do.call(rbind, X_); .Z <- do.call(rbind, Z_)

  n <- length(y); k <- ncol(.X); q <- ncol(Z); m <- length(y_); r <- sapply(D, length)
  id_ <- split(1:n, group)

  if(is.null(alpha)){alpha <- log(n)}
  s.t <- mean((y-mean(y))^2)

  if(is.null(cwu.control))
  {
    cwu.control <- list(conv.check = "coef", tol = 1e-5)
  } else
  {
    cwu.control <- list(
      conv.check = ifelse(is.null(cwu.control$conv.check), "coef", cwu.control$conv.check),
      tol = ifelse(is.null(cwu.control$tol), 1e-5, cwu.control$tol)
    )

    if(!(cwu.control$conv.check %in% c("coef", "grad")))
    {
      cwu.control$conv.check <- "coef"
      cwu.control$warning <- "`conv.check` is replaced by `coef`"
      warning("`conv.check` is replaced by `coef` because it is wrong")
    }
  }

  tol1 <- cwu.control$tol

  if(cwu.control$conv.check == "coef")
  {
    convC <- function(Beta, M, c, Lambda, B, r, Beta0){
      max((Beta - Beta0)^2 / Beta0^2, na.rm = TRUE)
    }
  } else
  {
    convC <- function(Beta, M, c, Lambda, B, r, Beta0){

      Grad <- 2*(drop(M%*%Beta) - c) + (
        lapply(1:r, function(.x){
          a <- Beta - B[.x,]
          return(Lambda[.x]*a/norm2(a))
        }) %>% Reduce("+", .)
      )

      return(abs(Grad) %>% max)
    }
  }

  fobj <- function(yli, Xli, Zli, D, BETA, Gamma, lambda, W, r){
    sapply(1:m, function(j){
      wj <- W[[j]]
      Dj <- D[[j]]
      Bj <- BETA[Dj,,drop=F]
      rj <- r[j]
      Betaj <- BETA[j,]

      rss <- sum((yli[[j]] - Xli[[j]]%*%BETA[j,] - Zli[[j]]%*%Gamma)^2)
      pen <- scale(Bj, center=Betaj, scale=F) %>% apply(1, norm2) %>%
        multiply_by(wj) %>% sum

      return(rss + lambda*pen)
    }) %>% sum
  }

  if(MPinv)
  {
    SOLVE <- svd_solve
  } else
  {
    SOLVE <- chol_solve
  }

  #=============================================================================
  ###   estimates for lambda.max
  #=============================================================================

  .XZ <- cbind(.X, .Z)

  Coef.max <- SOLVE(crossprod(.XZ), crossprod(.XZ, .y)) %>% drop
  Beta.max <- Coef.max[1:k]
  BETA.max <- matrix(Beta.max, m, k, byrow=T)
  Gamma.max <- Coef.max[(k+1):(k+q)]

  rm(.XZ); invisible(gc())

  #=============================================================================
  ###   LSE
  #=============================================================================

  GWcomp <- lapply(1:m, function(j){
    Xj <- X_[[j]]; X.X <- crossprod(Xj); X.y <- crossprod(Xj, y_[[j]])
    Z.X <- crossprod(Z_[[j]], Xj)
    Minv <- SOLVE(X.X)
    return(list(
      Z.X = Z.X, Beta0 = Minv%*%X.y %>% drop,
      MinvX.Z = tcrossprod(Minv, Z.X), X.X = X.X, X.y = X.y
    ))
  }) %>% purrr::transpose()
  Z.X_ <- GWcomp$Z.X
  BETA0 <- GWcomp$Beta0 %>% do.call(rbind, .)
  X.XinvX.Z_ <- GWcomp$MinvX.Z
  X.X_ <- GWcomp$X.X; X.y_ <- GWcomp$X.y
  rm(GWcomp); invisible(gc())

  .Z.Zinv <- SOLVE(crossprod(.Z))
  Gamma0 <- .Z.Zinv%*%crossprod(.Z, .y) %>% drop

  BETA.aft <- BETA.max
  Gamma.aft <- Gamma.max
  dif <- 1
  iter <- 0

  while(dif > tol)
  {
    iter <- iter + 1

    BETA.bef <- BETA.aft
    Gamma.bef <- Gamma.aft

    Betacomp <- lapply(1:m, function(j){
      Betaj <- BETA0[j,] - drop(X.XinvX.Z_[[j]]%*%Gamma.bef)
      return(list(
        Beta = Betaj, Z.XBeta = drop(Z.X_[[j]]%*%Betaj)
      ))
    }) %>% purrr::transpose()
    BETA.aft <- Betacomp$Beta %>% do.call(rbind, .)
    Z.XBeta <- Betacomp$Z.XBeta %>% Reduce(`+`, .)

    Gamma.aft <- Gamma0 - drop(.Z.Zinv%*%Z.XBeta)

    dif <- max(
      max(abs(BETA.aft - BETA.bef) / abs(BETA.bef)),
      max(abs(Gamma.aft - Gamma.bef) / abs(Gamma.bef))
    )
  } #end for while

  BETA.LSE <- BETA.aft
  Gamma.LSE <- Gamma.aft

  #=============================================================================
  ###   penalty weight
  #=============================================================================

  W <- lapply(1:m, function(j){
    Betaj <- BETA.LSE[j,]
    out <- sapply(D[[j]], function(l){1 / norm2(Betaj - BETA.LSE[l,])})
    return(out)
  })

  #=============================================================================
  ###   lambda
  #=============================================================================

  lam.max <- sapply(1:m, function(j){
    norm2(X.X_[[j]]%*%Beta.max - X.y_[[j]] + crossprod(Z.X_[[j]], Gamma.max)) /
      sum(W[[j]])
  }) %>% max

  if(Lambda == "default")
  {
    lamR <- sqrt(.Machine$double.eps)/lam.max

    Lambda <- lam.max * exp(
      seq(0, log(lamR), by = log(lamR)/99)
    ) %>% rev
  } else if(Lambda == "uniform")
  {
    Lambda <- seq(0, lam.max, length=101)[-1]
  }

  lam.n <- length(Lambda)

  ##############################################################################
  ###   CDA
  ##############################################################################

  gr.labs <- 1:m
  BETA.aft <- BETA.LSE
  Gamma.aft <- Gamma.LSE

  BETAs <- Gr.labs <- list()
  GAMMAs <- matrix(0, lam.n, q)
  LOG <- matrix(0, lam.n, 8)

  for(lam.i in 1:lam.n)
  {
    lambda <- Lambda[lam.i]

    dif <- 1
    iter <- 0
    while(dif > tol)
    {
      iter <- iter + 1
      BETA.bef <- BETA.aft
      Gamma.bef <- Gamma.aft

      ##########################################################################
      ###   update beta
      ##########################################################################

      X.yt <- lapply(1:m, function(j){
        X.y_[[j]] - drop(crossprod(Z.X_[[j]], Gamma.bef))
      })

      #=========================================================================
      ###   descent cycle
      #=========================================================================

      for(j in 1:m)
      {
        Dj <- D[[j]]
        labj <- gr.labs[Dj]

        Bj <- BETA.aft[Dj,,drop=F]
        wj <- W[[j]]
        rj <- r[j]

        if(length(unique(labj)) < rj)
        {
          idxj <- split(1:rj, labj)

          Bj0 <- Bj
          Bj <- lapply(idxj, function(x){Bj0[x[1],]}) %>% do.call(rbind, .)
          wj <- split(wj, labj) %>% sapply(., sum)
          labj <- sapply(idxj, function(x){labj[x[1]]})
        }

        resj <- GGFL1(X.X_[[j]], X.yt[[j]], 2*lambda*wj, Bj, BETA.aft[j,], tol1, convC, maxit)
        BETA.aft[j,] <- resj$solution

        if(resj$type == "join")
        {
          gr.labs[j] <- labj[resj$idx]
        }
        if(resj$type == "disjoin")
        {
          gr.labs[j] <- setdiff(1:m, gr.labs) %>% min
        }
      } #end for j

      #=========================================================================
      ###   fusion cycle
      #=========================================================================

      t <- unique(gr.labs) %>% length

      if(t == 1)
      {
        BETA.aft <- BETA.max
      } else if(t < m)
      {
        E <- split(1:m, gr.labs)
        idx.f <- (sapply(E, length) > 1) %>% which

        XI <- lapply(E, function(Ej){BETA.aft[Ej[1],]}) %>% do.call(rbind, .)
        grf.labs <- sapply(E, function(Ej){gr.labs[Ej[1]]})

        for(l in idx.f)
        {
          El <- E[[l]]
          Fl <- lapply(El, function(j){
            Dj <- D[[j]]
            return(setdiff(Dj, El))
          }) %>% unlist %>% unique
          D1 <- lapply((1:t)[-l], function(s){
            if(length(intersect(E[[s]], Fl)) > 0){return(s)}
          }) %>% unlist

          Jl <- lapply(D1, function(i){
            Ei <- E[[i]]

            Out <- lapply(El, function(j){
              Dj <- D[[j]]
              if(length(intersect(Ei, Dj)) > 0)
              {
                return(cbind(j, intersect(Ei, Dj)))
              }
            }) %>% do.call(rbind, .)
          })

          w1 <- lapply(Jl, function(x){
            lapply(1:nrow(x), function(i){
              W[[x[i,1]]][which(x[i,2]==D[[x[i,1]]])]
            }) %>% unlist %>% sum
          }) %>% unlist

          Bl <- XI[D1,,drop=F]
          rl <- length(D1)

          labl <- grf.labs[D1]

          if(length(unique(labl)) < rl)
          {
            idxl <- split(1:rl, labl)

            Bl0 <- Bl
            Bl <- lapply(idxl, function(x){Bl0[x[1],]}) %>% do.call(rbind, .)
            w1 <- split(w1, labl) %>% sapply(., sum)
            labl <- sapply(idxl, function(x){labl[x[1]]})
          }

          resl <- GGFL1(
            X.X_[El] %>% Reduce("+", .),
            X.yt[El] %>% Reduce("+", .),
            2*lambda*w1,
            Bl, XI[l,], tol1, convC, maxit
          )
          XI[l,] <- resl$solution
          BETA.aft[El,] <- matrix(XI[l,], length(El), k, byrow=T)

          if(resl$type == "join")
          {
            gr.labs[El] <- grf.labs[l] <- labl[resl$idx]
          }
          if(resl$type == "disjoin")
          {
            gr.labs[El] <- grf.labs[l] <- setdiff(1:m, grf.labs) %>% min
          }
        } #end for l
      } #end if t < m

      ##########################################################################
      ###   update Gamma
      ##########################################################################

      Z.XBeta <- lapply(1:m, function(j){
        drop(Z.X_[[j]]%*%BETA.aft[j,])
      }) %>% Reduce(`+`, .)

      Gamma.aft <- Gamma0 - drop(.Z.Zinv%*%Z.XBeta)

      ##########################################################################
      ###   convergence
      ##########################################################################

      dif <- max(
        max((BETA.aft - BETA.bef)^2 / BETA.bef^2, na.rm = TRUE),
        max((Gamma.aft - Gamma.bef)^2 / Gamma.bef^2, na.rm = TRUE)
      )

      if(progress)
      {
        obj <- fobj(y_, X_, Z_, D, BETA.aft, Gamma.aft, lambda, W, r)
        print(paste(lam.i, iter, dif, obj, sep="_"))
      }
    } #end while dif

    BETAs[[lam.i]] <- BETA.aft
    GAMMAs[lam.i,] <- Gamma.aft
    Gr.labs[[lam.i]] <- gr.labs

    y.hat <- (lapply(1:m, function(j){X_[[j]] %*% BETA.aft[j,] %>% drop}) %>% unlist) + .Z%*%Gamma.aft
    rss <- mean((.y-y.hat)^2)
    df <- length(unique(gr.labs))*k + q

    if(!progress)
    {
      obj <- fobj(y_, X_, Z_, D, BETA.aft, Gamma.aft, lambda, W, r)
    }

    LOG[lam.i,] <- c(
      lambda, rss, df, rss/((1 - (df/n))^alpha), 1 - (rss/s.t),
      median(abs(.y-y.hat)/abs(.y)), obj, iter
    )

  } #end for lam.i
  t2 <- proc.time()[3]

  LOG <- as.data.frame(LOG) %>% set_names(c(
    "lambda", "rss", "df", "MSC", "R2", "MER", "obje", "iter"
  ))

  opt <- which.min(LOG$MSC)
  BETA.hat <- BETAs[[opt]]
  Gamma.hat <- GAMMAs[opt,]

  FIT <- lapply(1:m, function(j){
    X_[[j]] %*% cbind(BETA.hat[j,], BETA.LSE[j,], Beta.max)
  }) %>% do.call(rbind, .) %>% add(.Z%*%cbind(Gamma.hat, Gamma.LSE, Gamma.max)) %>%
    set_colnames(c("GGFL", "OLS", "max")) %>%
    as.data.frame %>% dplyr::mutate(id = unlist(id_)) %>% dplyr::arrange(id) %>%
    dplyr::select(-id)

  if(standardize)
  {
    BETA.hat <- scale(BETA.hat, center=F, scale=normX) %>%
      `attr<-`("scaled:scale", NULL)
    BETA.LSE <- scale(BETA.LSE, center=F, scale=normX) %>%
      `attr<-`("scaled:scale", NULL)
    BETA.max <- matrix(Beta.max / normX, m, k, byrow=T)
  }

  out <- list(
    coefficients.local = list(
      GGFL = BETA.hat, OLS = BETA.LSE, MAX = BETA.max
    ),
    coefficients.global = list(
      GGFL = Gamma.hat, OLS = Gamma.LSE, MAX = Gamma.max
    ),
    fitted.values = FIT,
    summary = LOG[opt,] %>% dplyr::mutate(alpha=alpha, time=(t2-t1) %>% set_names(NULL)),
    weight = W,
    cluster = Gr.labs[[opt]],
    cwu.control = cwu.control
  )

  if(out.all)
  {
    if(standardize)
    {
      BETAs <- lapply(BETAs, function(BETA){
        scale(BETA, center=F, scale=normX) %>%
          `attr<-`("scaled:scale", NULL)
      })
    }

    out$all.local.coef = BETAs
    out$all.global.coef = GAMMAs
    out$all.cluster = do.call(rbind, Gr.labs) %>%
      set_colnames(paste0("g", 1:m)) %>% set_rownames(paste0("lam", 1:lam.n))
    out$log = LOG
  }

  return(out)
}
