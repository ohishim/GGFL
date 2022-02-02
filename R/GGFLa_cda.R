#' @title Coordinate optimization for GGFL
#' @description \code{GGFLa.cda} More general version of Coordinate optimization for GGFL.
#'   This allows to use explanatory variables which are not took GGFL penalty.
#'
#' @importFrom magrittr %>%
#'
#' @param yli list of group-wise vectors of an objective variable
#' @param Xli list of group-wise matrixes of explanatory variables which are took GGFL penalty
#' @param Zli list of group-wise matrixes of explanatory variables which are not took GGFL penalty
#' @param D list of adjacency relations
#' @param thres threshold for convergence judgement
#' @param progress If TRUE, progress is displayed
#' @param out.all if TRUE, results for all tuning parameters are output;
#'   if FALSE, results for only the optimal tuning parameter are output
#'
#' @return coef: list of GGFL estimates
#' @return rss: vector or scalar of residual sum of squares
#' @return df: vector or scalar of degrees of freedom
#' @return msc: vector or scalar of values of model selection criterion
#' @return r2: vector or scalar of coefficient of determination
#' @return mer: vector or scalar of median error rate
#' @return lambda: vector or scalar of tuning parameters
#' @return weight: list of penalty weight
#' @return gr.labs: list of labels for groups
#' @return time: runtime (s)
#'
#' @return coef.opt: matrix of GGFL estimates
#' @return coef.lse: matrix of OLS estimates
#' @return coef.max: matrix of estimates when all groups are equal
#' @return pred: vecter of predictive values for coef.opt
#' @return pred.lse: vecter of predictive values for coef.lse
#' @return pred.max: vecter of predictive values for coef.max
#'
#' @export
#' @examples
#' #GGFLa.cda(yli, Xli, Zli, D)

GGFLa.cda <- function(yli, Xli, Zli, D, thres=1e-5, progress=FALSE, out.all=FALSE){

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

  alpha <- log(n)
  s.t <- sum((y-mean(y))^2)/n

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

  X.XinvX.li <- lapply(1:m, function(j){
    solve(X.Xli[[j]], tol=10^(-100)) %*% X.li[[j]]
  })

  BETA.aft <- BETA.max
  Gamma.aft <- Gamma.max
  dif <- 1
  iter <- 0

  while(dif > thres)
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
  Lambda0 <- rev(lam.max * ((3/4)^((1:100)-1)))

  lam.check <- sapply(Lambda0, function(lambda){
    identical(all.equal(lambda, 0), TRUE)
  })
  Lambda <- Lambda0[!lam.check]
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
    while(dif > thres)
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

        Bj <- BETA.aft[Dj,]
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

        resj <- GGFL.cda1(X.Xli[[j]], X.yt[[j]], 2*lambda*wj, Bj, BETA.aft[j,])
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

          resl <- GGFL.cda1(
            X.Xli[El] %>% Reduce("+", .),
            X.yt[El] %>% Reduce("+", .),
            2*lambda*w1,
            Bl,
            XI[l,]
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
        max((BETA.aft - BETA.bef)^2) / max(BETA.bef^2),
        max((Gamma.aft - Gamma.bef)^2) / max(Gamma.bef^2)
      )

      if(progress){print(paste0(lam.i, "_", iter, "_", dif))}
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
      time = t2-t1
    )
  } else
  {
    opt <- which.min(MSC)
    BETA.hat <- BETAs[[opt]]
    Gamma.hat <- GAMMAs[opt,]

    out <- list(
      coef1.opt= BETA.hat,
      coef1.lse=BETA.LSE,
      coef1.max=BETA.max,
      coef2.opt=Gamma.hat,
      coef2.lse=Gamma.LSE,
      coef2.max=Gamma.max,
      pred = (lapply(1:m, function(j){Xli[[j]] %*% BETA.hat[j,] %>% drop}) %>% unlist) +
        Z %*% Gamma.hat,
      pred.lse = (lapply(1:m, function(j){Xli[[j]] %*% BETA.LSE[j,] %>% drop}) %>% unlist) +
        Z %*% Gamma.LSE,
      pred.max = (lapply(1:m, function(j){Xli[[j]] %*% BETA.max[j,] %>% drop}) %>% unlist) +
        Z %*% Gamma.max,
      rss = RSS[opt],
      df = DF[opt],
      msc = MSC[opt],
      r2 = R2[opt],
      mer = MER[opt],
      lambda = Lambda[opt],
      weight = W,
      gr.labs = Gr.labs[[opt]],
      time = t2-t1
    )
  }

  return(out)
}
