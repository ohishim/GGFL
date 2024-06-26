#' @title Coordinate optimization for GGFL
#' @description \code{GGFL} solving GGFL optimization problem via coordinate descent algorithm (v0.2.2)
#'
#' @importFrom magrittr %>% multiply_by
#' @importFrom MASS ginv
#' @importFrom stats median
#'
#' @param y vector of an objective variable
#' @param X matrix of explanatory variables without intercept
#' @param area vector of area labels for each individual
#' @param adj dataframe of adjacent information of which
#'   the labels for the first and second column are respectively "area" and "adj"
#' @param Lambda candidates of tuning parameter
#'   if `"default"`, the candidates are defined following the paper;
#'   if `"uniform"`, the candidates are defined by uniformly dividing;
#'   if vector/scalar, the candidates are defined by following `lambda.type`
#' @param lambda.type option when `Lambda` is numeric;
#'   if "value", `Lambda` is searching points;
#'   if "rate", `lam.max*Lambda` is searching points
#' @param adaptive option only when `weight=NULL`; if `TRUE`, adaptive penalty is adopted
#' @param weight list of penalty weight (default is `NULL`)
#' @param alpha a value expressing penalty strength for EGCV criterion;
#'   default is `NULL` which corresponds to `log(n)`
#' @param standardize if `TRUE`, `y` and `X` are standardized in the sense of norm
#' @param MPinv if `TRUE`, the ordinary least squares estimator is calculated by the Moore-Penrose inverse matrix
#' @param tol tolerance for convergence
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
#' \item{idx}{the index of the optimal tuning parameter}
#'
#' \item{coefGGFL}{matrix of GGFL estimates}
#'
#' \item{coefOLS}{matrix of OLS estimates}
#'
#' \item{coefMAX}{matrix of estimates when all groups are equal}
#'
#' \item{fitGGFL}{vecter of fitted values for `coefGGFL`}
#'
#' \item{fitOLS}{vecter of fitted values for `coefOLS`}
#'
#' \item{fitMAX}{vecter of fitted values for `coefMAX`}
#'
#' \item{Coef}{list of GGFL estimates}
#'
#' \item{obje}{objective function}
#'
#' \item{rss}{residual sum of squares}
#'
#' \item{df}{degrees of freedom}
#'
#' \item{msc}{EGCV criterion}
#'
#' \item{r2}{coefficient of determination}
#'
#' \item{mer}{median error rate}
#'
#' \item{iter}{the number of iteration}
#'
#' \item{lambda}{tuning parameter}
#'
#' \item{weight}{list of penalty weight}
#'
#' \item{gr.labs}{list of labels for groups}
#'
#' \item{time}{runtime (in seconds)}
#'
#' \item{cwu.control}{`cwu.control`}
#'
#' @export
#' @examples
#' #GGFL(y, X, area, adj)

GGFL <- function(
  y, X, area, adj, Lambda="default", lambda.type=NULL, adaptive=TRUE, weight=NULL, alpha=NULL,
  standardize=TRUE, MPinv=FALSE, tol=1e-5, maxit=500, progress=FALSE, out.all=FALSE,
  cwu.control=NULL
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

  if(is.null(weight))
  {
    if(adaptive)
    {
      W <- lapply(1:m, function(j){
        Betaj <- BETA.LSE[j,]
        out <- sapply(D[[j]], function(l){1 / norm2(Betaj - BETA.LSE[l,])})
        return(out)
      })
    } else
    {
      W <- lapply(r, function(x){rep(1, x)})
    }
  } else
  {
    W <- weight
    rm(weight)
  }

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
  s.t <- sum((y-mean(y))^2)/n

  if(is.null(alpha)){alpha <- log(n)}

  gr.labs <- 1:m

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

  fobj <- function(rss, D, BETA, lambda, W){

    pen <- sapply(1:m, function(j){
      wj <- W[[j]]
      Dj <- D[[j]]
      Bj <- BETA[Dj,,drop=F]
      Betaj <- BETA[j,]

      out <- scale(Bj, center=Betaj, scale=F) %>% apply(1, norm2) %>%
        multiply_by(wj) %>% sum

      return(out)
    }) %>% sum

    return(rss + lambda*pen)
  }

  ##############################################################################
  ###   CDA
  ##############################################################################

  BETA.aft <- BETA.LSE

  BETAs <- Gr.labs <- list()
  RSS <- DF <- MSC <- R2 <- MER <- Obje <- IterN <- numeric(lam.n)

  for(lam.i in 1:lam.n)
  {
    lambda <- Lambda[lam.i]

    dif <- 1
    iter <- 0
    while((dif > tol) & (iter < maxit))
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

        resj <- GGFL1(X.Xli[[j]], X.y[[j]], 2*lambda*wj, Bj, BETA.aft[j,], tol1, convC, maxit)
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

          resl <- GGFL1(
            X.Xli[El] %>% Reduce("+", .),
            X.y[El] %>% Reduce("+", .),
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
      ###   convergence
      ##########################################################################

      dif <- max((BETA.aft - BETA.bef)^2 / BETA.bef^2, na.rm = TRUE)

      if(progress)
      {
        yli.hat <- lapply(1:m, function(j){Xli[[j]] %*% BETA.aft[j,] %>% drop})
        .rss <- sapply(1:m, function(j){sum((yli[[j]] - yli.hat[[j]])^2)}) %>% sum
        obj <- fobj(.rss, D, BETA.aft, lambda, W)

        print(paste(lam.i, iter, dif, obj, sep="_"))
      }
    } #end while dif

    BETAs[[lam.i]] <- BETA.aft
    Gr.labs[[lam.i]] <- match(gr.labs, unique(gr.labs))

    if(progress)
    {
      RSS[lam.i] <- rss <- .rss/n
      Obje[lam.i] <- obj
    } else
    {
      yli.hat <- lapply(1:m, function(j){Xli[[j]] %*% BETA.aft[j,] %>% drop})
      RSS[lam.i] <- rss <- (sapply(1:m, function(j){sum((yli[[j]] - yli.hat[[j]])^2)}) %>% sum) / n
      Obje[lam.i] <- fobj(n*rss, D, BETA.aft, lambda, W)
    }

    DF[lam.i] <- df <- length(unique(gr.labs))*k

    MSC[lam.i] <- rss / ((1 - (df/n))^alpha)

    R2[lam.i] <- 1 - (rss/s.t)
    MER[lam.i] <- median(abs(y-unlist(yli.hat))/abs(y))
    IterN[lam.i] <- iter

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
      coefGGFL = BETA.hat,
      coefOLS = BETA.LSE,
      coefMAX = BETA.max,
      fitGGFL = pred,
      fitOLS = pred.lse,
      fitMAX = pred.max,
      Coef = BETAs,
      obje = Obje,
      rss = RSS,
      df = DF,
      msc = MSC,
      r2 = R2,
      mer = MER,
      iter = IterN,
      lambda = Lambda,
      weight = W,
      gr.labs = Gr.labs,
      time = t2-t1,
      cwu.control = cwu.control
    )
  } else
  {
    out <- list(
      coefGGFL = BETA.hat,
      coefOLS = BETA.LSE,
      coefMAX = BETA.max,
      fitGGFL = pred,
      fitOLS = pred.lse,
      fitMAX = pred.max,
      obje = Obje[opt],
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
