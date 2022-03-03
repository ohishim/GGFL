#' @title Coordinate optimization for GGFL
#' @description \code{GGFL.cda} solving GGFL optimization problem via coordinate descent algorithm
#'
#' @importFrom magrittr %>%
#' @importFrom MASS ginv
#'
#' @param y vector of an objective variable
#' @param X matrix of explanatory variables without intercept
#' @param area vector of area labels for each individual
#' @param adj dataframe of adjacent information of which
#'   the labels for the first and second column are respectively "area" and "adj"
#' @param thres threshold for convergence judgement
#' @param Lambda candidates of tuning parameter
#'   if "default", the candidates are defined following the paper;
#'   if "uniform", the candidates are defined by uniformly dividing;
#'   if vector/scalar, the candidates are defined by following lambda.type
#' @param lambda.type option when Lambda is numeric;
#'   if "value", Lambda is searching points;
#'   if "rate", lam.max*Lambda is searching points
#' @param standardize if TRUE, y and X are standardized in the sense of norm
#' @param MPinv if TRUE, the ordinary least square estimator is calculated by the Moore–Penrose inverse matrix
#' @param progress If TRUE, progress is displayed
#' @param out.all if TRUE, results for all tuning parameters are output;
#'   if FALSE, results for only the optimal tuning parameter are output
#'
#' @return idx: the index of the optimal tuning parameter
#' @return Coef: list of GGFL estimates
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
#' #GGFL.cda(y, X, area, adj)

GGFL.cda <- function(
  y, X, area, adj, thres=1e-5, Lambda="default", lambda.type=NULL,
  standardize=TRUE, MPinv=FALSE, progress=FALSE, out.all=FALSE
){

  ##############################################################################
  ###   preparations
  ##############################################################################

  t1 <- proc.time()[3]

  norm2 <- function(x){sqrt(sum(x^2))}

  if(standardize)
  {
    y0 <- y
    normy <- norm2(y)
    y <- y0 / normy
    yli <- split(y, area)

    X0 <- X
    normX <- c(1, apply(X0, 2, norm2))
    X <- apply(X0, 2, function(x){x / norm2(x)}) %>% cbind(1, .)
    Xli <- as.data.frame(X) %>% split(area) %>% lapply(data.matrix)
  } else
  {
    yli <- split(y, area)
    Xli <- cbind(1, X) %>% as.data.frame %>% split(area) %>% lapply(data.matrix)
  }

  area1 <- unique(area) %>% sort
  D <- split(adj$adj, adj$area) %>% lapply(function(x){which(area1 %in% x)})

  n <- sapply(yli, length) %>% sum
  k <- Xli[[1]] %>% ncol
  m <- length(yli)

  X <- do.call(rbind, Xli)
  X. <- t(X)
  M <- X. %*% X
  Minv <- solve(M)
  MinvX. <- Minv %*% X.

  rm(Minv, X, X.)

  X.Xli <- lapply(Xli, function(A){t(A)%*%A})

  if(MPinv)
  {
    X.Xinvli <- lapply(X.Xli, ginv)
  } else
  {
    X.Xinvli <- lapply(X.Xli, solve)
  }

  r <- sapply(D, length)

  X.y <- lapply(1:m, function(j){t(Xli[[j]]) %*% yli[[j]] %>% drop})

  BETA.LSE <- lapply(1:m, function(j){X.Xinvli[[j]] %*% X.y[[j]] %>% drop}) %>% do.call(rbind, .)
  Beta.max <- MinvX. %*% y %>% drop
  BETA.max <- lapply(1:m, function(j){Beta.max}) %>% do.call(rbind, .)

  rm(MinvX., X.Xinvli)

  W <- lapply(1:m, function(j){
    Betaj <- BETA.LSE[j,]
    out <- sapply(D[[j]], function(l){1 / norm2(Betaj - BETA.LSE[l,])})
    return(out)
  })

  if(is.character(Lambda))
  {
    lam.max <- sapply(1:m, function(j){norm2(X.Xli[[j]]%*%Beta.max - X.y[[j]]) / sum(W[[j]])}) %>% max

    if(Lambda == "default")
    {
      Lambda <- rev(lam.max * ((3/4)^((1:100)-1)))
    } else if(Lambda == "uniform")
    {
      Lambda <- seq(0, lam.max, length=101)[-1]
    }
  } else if(is.null(lambda.type))
  {
    stop("Please select lambda.type")
  } else if(lambda.type == "rate")
  {
    lam.max <- sapply(1:m, function(j){norm2(X.Xli[[j]]%*%Beta.max - X.y[[j]]) / sum(W[[j]])}) %>% max
    Lambda <- lam.max * Lambda
  } else if(lambda.type != "value")
  {
    stop("The lambda.type is missing")
  }

  lam.n <- length(Lambda)

  alpha <- log(n)
  s.t <- sum((y-mean(y))^2)/n

  gr.labs <- 1:m

  ##############################################################################
  ###   CDA
  ##############################################################################

  BETA.aft <- BETA.LSE

  BETAs <- Gr.labs <- list()
  RSS <- DF <- MSC <- R2 <- MER <- numeric(lam.n)

  for(lam.i in 1:lam.n)
  {
    lambda <- Lambda[lam.i]

    dif <- 1
    iter <- 0
    while(dif > thres)
    {
      iter <- iter + 1
      BETA.bef <- BETA.aft

      ##########################################################################
      ###   descent cycle
      ##########################################################################

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

        resj <- GGFL.cda1(X.Xli[[j]], X.y[[j]], 2*lambda*wj, Bj, BETA.aft[j,])
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

      ##########################################################################
      ###   fusion cycle
      ##########################################################################

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
            X.y[El] %>% Reduce("+", .),
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
      ###   convergence
      ##########################################################################

      dif <- max((BETA.aft - BETA.bef)^2) / max(BETA.bef^2)

      if(progress){print(paste0(lam.i, "_", iter, "_", dif))}
    } #end while dif

    BETAs[[lam.i]] <- BETA.aft
    Gr.labs[[lam.i]] <- match(gr.labs, unique(gr.labs))

    yli.hat <- lapply(1:m, function(j){Xli[[j]] %*% BETA.aft[j,] %>% drop})
    RSS[lam.i] <- rss <- (sapply(1:m, function(j){sum((yli[[j]] - yli.hat[[j]])^2)}) %>% sum) / n
    DF[lam.i] <- df <- length(unique(gr.labs))*k

    MSC[lam.i] <- rss / ((1 - (df/n))^alpha)

    R2[lam.i] <- 1 - (rss/s.t)
    MER[lam.i] <- median(abs(y-unlist(yli.hat))/abs(y))

  } #end for lam.i
  t2 <- proc.time()[3]

  opt <- which.min(MSC)
  BETA.hat <- BETAs[[opt]]
  pred <- lapply(1:m, function(j){Xli[[j]] %*% BETA.hat[j,] %>% drop}) %>% unlist
  pred.lse <- lapply(1:m, function(j){Xli[[j]] %*% BETA.LSE[j,] %>% drop}) %>% unlist
  pred.max <- lapply(1:m, function(j){Xli[[j]] %*% BETA.max[j,] %>% drop}) %>% unlist

  if(standardize)
  {
    BETA.hat <- normy * BETA.hat / matrix(normX, m, k, byrow=T)
    BETA.LSE <- normy * BETA.LSE / matrix(normX, m, k, byrow=T)
    BETA.max <- matrix(normy * Beta.max / normX, m, k, byrow=T)

    pred <- pred * normy
    pred.lse <- pred.lse * normy
    pred.max <- pred.max * normy
  }

  if(out.all)
  {
    if(standardize)
    {
      BETAs <- lapply(BETAs, function(BETA){
        normy * BETA / matrix(normX, m, k, byrow=T)
      })
    }

    out <- list(
      idx = opt,
      coef.opt = BETA.hat,
      coef.lse = BETA.LSE,
      coef.max = BETA.max,
      pred = pred,
      pred.lse = pred.lse,
      pred.max = pred.max,
      Coef = BETAs,
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
    out <- list(
      coef.opt = BETA.hat,
      coef.lse = BETA.LSE,
      coef.max = BETA.max,
      pred = pred,
      pred.lse = pred.lse,
      pred.max = pred.max,
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
