#' @title Coordinate optimization for GGFL
#' @description \code{GGFL} solving GGFL optimization problem via coordinate descent algorithm (v1.0.2)
#'
#' @importFrom dplyr arrange mutate select
#' @importFrom magrittr %>% multiply_by set_colnames set_names
#' @importFrom purrr transpose
#' @importFrom stats median
#'
#' @param y vector of an objective variable
#' @param X matrix of explanatory variables without intercept
#' @param group numeric vector of group labels for each individual
#' @param adj data.frame or matrix with two columns of adjacent information
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
#' \item{coefficients}{a list object which has
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
#' \item{cluster}{a list of cluster}
#'
#' \item{cluster.labels}{a vector expressing which cluster each group is in}
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
#' #GGFL(y, X, group, adj)

GGFL <- function(
  y, X, group, adj, Lambda="default", lambda.type=NULL, adaptive=TRUE, weight=NULL, alpha=NULL,
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
    normy <- norm2(y)
    y <- y / normy
    y_ <- split(y, group)

    normX <- c(1, colSums(X^2) %>% sqrt)
    X <- scale(cbind(1, X), center=F, scale=normX) %>% `attr<-`("scaled:scale", NULL)
    X_ <- as.data.frame(X) %>% split(group) %>% lapply(data.matrix)
  } else
  {
    y_ <- split(y, group)
    X_ <- cbind(1, X) %>% as.data.frame %>% split(group) %>% lapply(data.matrix)
  }

  .y <- unlist(y_)
  D <- split(adj[,2], adj[,1]); r <- sapply(D, length)

  n <- length(y); k <- ncol(X_[[1]]); m <- length(y_)

  id_ <- split(1:n, group)

  if(MPinv)
  {
    SOLVE <- svd_solve
  } else
  {
    SOLVE <- chol_solve
  }

  GWcomp <- lapply(1:m, function(j){
    Xj <- X_[[j]]; yj <- y_[[j]]
    X.X <- crossprod(Xj); X.y <- crossprod(Xj, yj) %>% drop

    return(list(
      X.X = X.X, X.y = X.y, LSE = SOLVE(X.X, X.y)
    ))
  }) %>% purrr::transpose()
  X.X_ <- GWcomp$X.X; X.y_ <- GWcomp$X.y
  BETA.LSE <- GWcomp$LSE %>% do.call(rbind, .)

  rm(GWcomp); invisible(gc())

  Beta.max <- SOLVE(Reduce(`+`, X.X_), Reduce(`+`, X.y_))
  BETA.max <- matrix(Beta.max, m, k, byrow=T)

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
    rm(weight); invisible(gc())
  }

  if(is.character(Lambda))
  {
    lam.max <- sapply(1:m, function(j){norm2(X.X_[[j]]%*%Beta.max - X.y_[[j]]) / sum(W[[j]])}) %>% max

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
  } else if(is.null(lambda.type))
  {
    stop("Please select lambda.type")
  } else if(lambda.type == "rate")
  {
    lam.max <- sapply(1:m, function(j){norm2(X.X_[[j]]%*%Beta.max - X.y_[[j]]) / sum(W[[j]])}) %>% max
    Lambda <- lam.max * Lambda
  } else if(lambda.type != "value")
  {
    stop("The lambda.type is missing")
  }

  lam.n <- length(Lambda)
  s.t <- mean((y-mean(y))^2)

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
  # RSS <- DF <- MSC <- R2 <- MER <- Obje <- IterN <- numeric(lam.n)
  LOG <- matrix(0, lam.n, 8)

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

        resj <- GGFL1(X.X_[[j]], X.y_[[j]], 2*lambda*wj, Bj, BETA.aft[j,], tol1, convC, maxit)
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
            X.X_[El] %>% Reduce("+", .),
            X.y_[El] %>% Reduce("+", .),
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
        y.hat_ <- lapply(1:m, function(j){X_[[j]] %*% BETA.aft[j,] %>% drop})
        .rss <- sapply(1:m, function(j){sum((y_[[j]] - y.hat_[[j]])^2)}) %>% sum
        obj <- fobj(.rss, D, BETA.aft, lambda, W)

        print(paste(lam.i, iter, dif, obj, sep="_"))
      }
    } #end while dif

    BETAs[[lam.i]] <- BETA.aft
    Gr.labs[[lam.i]] <- match(gr.labs, unique(gr.labs))

    if(progress)
    {
      rss <- .rss/n
    } else
    {
      y.hat_ <- lapply(1:m, function(j){X_[[j]] %*% BETA.aft[j,] %>% drop})
      rss <- (sapply(1:m, function(j){sum((y_[[j]] - y.hat_[[j]])^2)}) %>% sum) / n
      obj <- fobj(n*rss, D, BETA.aft, lambda, W)
    }

    df <- length(unique(gr.labs))*k
    LOG[lam.i,] <- c(
      lambda, rss, df, rss / ((1 - (df/n))^alpha), 1 - (rss/s.t),
      median(abs(.y-unlist(y.hat_))/abs(.y)), obj, iter
    )
  } #end for lam.i
  t2 <- proc.time()[3]

  LOG <- as.data.frame(LOG) %>% set_names(c(
    "lambda", "rss", "df", "MSC", "R2", "MER", "obje", "iter"
  ))

  opt <- which.min(LOG$MSC)
  BETA.hat <- BETAs[[opt]]

  FIT <- lapply(1:m, function(j){
    X_[[j]]%*%cbind(BETA.hat[j,], BETA.LSE[j,], Beta.max)
  }) %>% do.call(rbind, .) %>% set_colnames(c("GGFL", "OLS", "max")) %>%
    as.data.frame %>% dplyr::mutate(id = unlist(id_)) %>% dplyr::arrange(id) %>%
    dplyr::select(-id)

  if(standardize)
  {
    BETA.hat <- scale(normy*BETA.hat, center=F, scale=normX) %>%
      `attr<-`("scaled:scale", NULL)
    BETA.LSE <- scale(normy*BETA.LSE, center=F, scale=normX) %>%
      `attr<-`("scaled:scale", NULL)
    BETA.max <- matrix(normy * Beta.max / normX, m, k, byrow=T)

    FIT <- normy*FIT
  }

  cluster <- Gr.labs[[opt]]
  m. <- max(cluster)
  rownames(BETA.hat) <- paste0("c", cluster)

  out <- list(
    coefficients = list(
      GGFL = BETA.hat, OLS = BETA.LSE, MAX = BETA.max
    ),
    fitted.values = FIT,
    summary = LOG[opt,] %>% dplyr::mutate(
      cluster = m., alpha = alpha, time = (t2-t1) %>% set_names(NULL)
    ),
    weight = W,
    cluster = split(1:m, cluster) %>% set_names(paste0("c", 1:m.)),
    cluster.labels = cluster %>% set_names(paste0("g", 1:m)),
    cwu.control = cwu.control
  )

  if(out.all)
  {
    if(standardize)
    {
      BETAs <- lapply(BETAs, function(BETA){
        scale(normy*BETA, center=F, scale=normX) %>%
          `attr<-`("scaled:scale", NULL)
      })
    }

    out$all.coef = BETAs
    out$all.cluster = do.call(rbind, Gr.labs) %>%
      set_colnames(paste0("g", 1:m)) %>% set_rownames(paste0("lam", 1:lam.n))
    out$log = LOG
  }

  return(out)
}
