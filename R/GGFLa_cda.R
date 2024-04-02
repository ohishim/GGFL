#' @title Coordinate optimization for GGFL
#' @description \code{pGGFL} More general version of Coordinate optimization for GGFL.
#'   Only partial explanatory variables are penalized. (v0.2.1)
#'
#' @importFrom magrittr %>%
#' @importFrom MASS ginv
#'
#' @param yli list of group-wise vectors of an objective variable
#' @param Xli list of group-wise matrixes of explanatory variables which are took GGFL penalty
#' @param Zli list of group-wise matrixes of explanatory variables which are not took GGFL penalty
#' @param D list of adjacency relations
#' @param tol tolerance for convergence
#' @param Lambda candidates of tuning parameter
#'   if `"default"`, the candidates are defined following the paper;
#'   if `"uniform"`, the candidates are defined by uniformly dividing
#' @param alpha a value expressing penalty strength for EGCV criterion;
#'   default is `NULL` which corresponds to `log(n)`
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
#' \item{coef1}{list of GGFL estimates for penalized variables}
#'
#' \item{coef2}{list of GGFL estimates for non-penalized variables}
#'
#' \item{rss}{vector or scalar of residual sum of squares}
#'
#' \item{df}{vector or scalar of degrees of freedom}
#'
#' \item{msc}{vector or scalar of values of model selection criterion}
#'
#' \item{r2}{vector or scalar of coefficient of determination}
#'
#' \item{mer}{vector or scalar of median error rate}
#'
#' \item{lambda}{vector or scalar of tuning parameters}
#'
#' \item{weight}{list of penalty weight}
#'
#' \item{gr.labs}{list of labels for groups}
#'
#' \item{time}{runtime (s)}
#'
#' \item{coef1GGFL}{matrix of GGFL estimates for penalized variables}
#'
#' \item{coef1OLS}{matrix of OLS estimates for penalized variables}
#'
#' \item{coef1MAX}{matrix of estimates for penalized variables when all groups are equal}
#'
#' \item{coef2GGFL}{matrix of GGFL estimates for non-penalized variables}
#'
#' \item{coef2OLS}{matrix of OLS estimates for non-penalized variables}
#'
#' \item{coef2MAX}{matrix of estimates for non-penalized variables when all groups are equal}
#'
#' \item{fitGGFL}{vector of fitted values for `coef1GGFL` and `coef2GGFL`}
#'
#' \item{fitOLS}{vector of fitted values for `coef1OLS` and `coef2OLS`}
#'
#' \item{fitMAX}{vector of fitted values for `coef1MAX` and `coef2MAX`}
#'
#' \item{cwu.control}{`cwu.control`}
#'
#' @export
#' @examples
#' #pGGFL(yli, Xli, Zli, D)

pGGFL <- function(
  yli, Xli, Zli, D, tol=1e-5, Lambda="default", alpha=NULL,
  MPinv=FALSE, maxit=500, progress=FALSE, out.all=FALSE, cwu.control=NULL
){

  ##############################################################################
  ###   preparation
  ##############################################################################

  norm2 <- function(x){sqrt(sum(x^2))}

  y <- unlist(yli)

  n <- length(y)
  k <- Xli[[1]] %>% ncol
  q <- Zli[[1]] %>% ncol
  m <- length(yli)
  r <- sapply(D, length)

  Z <- do.call(rbind, Zli)
  Z. <- t(Z)
  Z.Zinv <- solve(Z. %*% Z)

  X.li <- lapply(Xli, function(X){t(X)})
  X.Xli <- lapply(1:m, function(j){X.li[[j]] %*% Xli[[j]]})

  ns <- sapply(yli, length)
  gr.idx <- mapply(rep, 1:m, ns) %>% unlist

  if(is.null(alpha)){alpha <- log(n)}
  s.t <- sum((y-mean(y))^2)/n

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

  if(progress)
  {
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
  }

  #=============================================================================
  ###   estimates for lambda.max
  #=============================================================================

  XZ <- cbind(do.call(rbind, Xli), Z); ZX <- t(XZ)

  Coef.max <- solve(ZX %*% XZ) %*% ZX %*% y %>% drop
  Beta.max <- Coef.max[1:k]
  BETA.max <- lapply(1:m, function(j){Beta.max}) %>% do.call(rbind, .)
  Gamma.max <- Coef.max[(k+1):(k+q)]

  rm(XZ, ZX, Coef.max)

  #=============================================================================
  ###   LSE
  #=============================================================================

  Z.ZinvZ. <- Z.Zinv %*% Z.

  if(MPinv)
  {
    X.XinvX.li <- lapply(1:m, function(j){
      ginv(X.Xli[[j]]) %*% X.li[[j]]
    })
  } else
  {
    X.XinvX.li <- lapply(1:m, function(j){
      solve(X.Xli[[j]], tol=10^(-100)) %*% X.li[[j]]
    })
  }

  BETA.aft <- BETA.max
  Gamma.aft <- Gamma.max
  dif <- 1
  iter <- 0

  while(dif > tol)
  {
    iter <- iter + 1

    BETA.bef <- BETA.aft
    Gamma.bef <- Gamma.aft

    BETA.aft <- lapply(1:m, function(j){
      X.XinvX.li[[j]] %*% (yli[[j]] - Zli[[j]]%*%Gamma.bef) %>% drop
    }) %>% do.call(rbind, .)

    XBeta <- lapply(1:m, function(j){Xli[[j]] %*% BETA.aft[j,]}) %>% unlist
    Gamma.aft <- Z.ZinvZ. %*% (y - XBeta) %>% drop

    dif <- max(
      max(abs(BETA.aft - BETA.bef) / abs(BETA.bef)),
      max(abs(Gamma.aft - Gamma.bef) / abs(Gamma.bef))
    )
  } #end for while

  BETA.LSE <- BETA.aft
  Gamma.LSE <- Gamma.aft

  rm(Z.ZinvZ., X.XinvX.li)

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

  ytli <- split(y - Z%*%Gamma.max, gr.idx)

  lam.max <- sapply(1:m, function(j){norm2(X.Xli[[j]]%*%Beta.max - X.li[[j]]%*%ytli[[j]]) / sum(W[[j]])}) %>% max

  if(Lambda == "default")
  {
    Lambda0 <- rev(lam.max * ((3/4)^((1:100)-1)))

    lam.check <- sapply(Lambda0, function(lambda){
      identical(all.equal(lambda, 0), TRUE)
    })
    Lambda <- Lambda0[!lam.check]
  } else if(Lambda == "uniform")
  {
    Lambda <- seq(0, lam.max, length=101)[-1]
  }

  lam.n <- length(Lambda)

  rm(ytli)

  ##############################################################################
  ###   CDA
  ##############################################################################

  gr.labs <- 1:m
  BETA.aft <- BETA.LSE
  Gamma.aft <- Gamma.LSE

  BETAs <- Gr.labs <- list()
  GAMMAs <- matrix(0, lam.n, q)
  RSS <- DF <- MSC <- R2 <- MER <- MER <- numeric(lam.n)

  t1 <- proc.time()[3]
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

      ytli <- (y - Z%*%Gamma.bef) %>% split(., gr.idx)
      X.yt <- lapply(1:m, function(j){
        X.li[[j]] %*% ytli[[j]] %>% drop
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

        resj <- GGFL1(X.Xli[[j]], X.yt[[j]], 2*lambda*wj, Bj, BETA.aft[j,], tol1, convC, maxit)
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
            X.Xli[El] %>% Reduce("+", .),
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

      yt <- lapply(1:m, function(j){
        yli[[j]] - Xli[[j]]%*%BETA.aft[j,]
      }) %>% unlist
      Gamma.aft <- Z.Zinv %*% Z. %*% yt %>% drop

      ##########################################################################
      ###   convergence
      ##########################################################################

      dif <- max(
        max((BETA.aft - BETA.bef)^2 / BETA.bef^2, na.rm = TRUE),
        max((Gamma.aft - Gamma.bef)^2 / Gamma.bef^2, na.rm = TRUE)
      )

      if(progress)
      {
        print(paste(lam.i, iter, dif, fobj(yli, Xli, Zli, D, BETA.aft, Gamma.aft, lambda, W, r), sep="_"))
      }
    } #end while dif

    BETAs[[lam.i]] <- BETA.aft
    GAMMAs[lam.i,] <- Gamma.aft
    Gr.labs[[lam.i]] <- gr.labs

    y.hat <- (lapply(1:m, function(j){Xli[[j]] %*% BETA.aft[j,] %>% drop}) %>% unlist) + Z%*%Gamma.aft
    RSS[lam.i] <- rss <- sum((y-y.hat)^2) / n
    DF[lam.i] <- df <- length(unique(gr.labs))*k + q

    MSC[lam.i] <- rss / ((1 - (df/n))^alpha)

    R2[lam.i] <- 1 - (rss/s.t)
    MER[lam.i] <- median(abs(y-y.hat)/abs(y))

  } #end for lam.i
  t2 <- proc.time()[3]

  if(out.all)
  {
    out <- list(
      coef1 = BETAs,
      coef2 = GAMMAs,
      rss = RSS,
      df = DF,
      msc = MSC,
      r2 = R2,
      mer = MER,
      lambda = Lambda,
      weight = W,
      gr.labs = Gr.labs,
      time = t2-t1,
      cwu.control = cwu.control
    )
  } else
  {
    opt <- which.min(MSC)
    BETA.hat <- BETAs[[opt]]
    Gamma.hat <- GAMMAs[opt,]

    out <- list(
      coef1GGFL=BETA.hat,
      coef1OLS=BETA.LSE,
      coef1MAX=BETA.max,
      coef2GGFL=Gamma.hat,
      coef2OLS=Gamma.LSE,
      coef2MAX=Gamma.max,
      fitGGFL = (lapply(1:m, function(j){Xli[[j]] %*% BETA.hat[j,] %>% drop}) %>% unlist) +
        Z %*% Gamma.hat,
      fitOLS = (lapply(1:m, function(j){Xli[[j]] %*% BETA.LSE[j,] %>% drop}) %>% unlist) +
        Z %*% Gamma.LSE,
      fitMAX = (lapply(1:m, function(j){Xli[[j]] %*% BETA.max[j,] %>% drop}) %>% unlist) +
        Z %*% Gamma.max,
      rss = RSS[opt],
      df = DF[opt],
      msc = MSC[opt],
      r2 = R2[opt],
      mer = MER[opt],
      lambda = Lambda[opt],
      weight = W,
      gr.labs = Gr.labs[[opt]],
      time = t2-t1,
      cwu.control = cwu.control
    )
  }

  return(out)
}
