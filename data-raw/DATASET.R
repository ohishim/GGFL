## code to prepare `exdata` dataset goes here

library(magrittr)

rm(list=ls(all.names = T));
gc();gc();

set.seed(22)

n0 <- 50
k <- 20
m <- 10
n <- n0*m
rhox <- 0.5

F.sta <- list(c(1,2,3),c(4,5,6,9,10),c(7,8))
m.sta <- length(F.sta)

PSI <- diag(k)
for(i in 1:(k-1))
{
  for(j in (i+1):(k))
  {
    PSI[i, j] <- PSI[j, i] <- rhox^{j-i}
  }
}

svdPSI <- svd(PSI)
PSI.sq <- svdPSI$u %*% diag(sqrt(svdPSI$d), k, k) %*% t(svdPSI$v)

BETA.sta <- matrix(0, m, k)
for(j in 1:m.sta)
{
  BETA.sta[F.sta[[j]], ] <- j
} #end for j

X0 <- runif(n*k, -1, 1) %>% matrix(., n, k)
X <- (X0%*%PSI.sq) %>% as.data.frame
names(X) <- paste0("X", 1:k)

area <- rep(1:m, each=n0)
X1 <- split(X, area) %>% lapply(data.matrix)

y <- lapply(1:m, function(j){
  Xj <- X1[[j]]
  Betaj <- BETA.sta[j,]
  return(Xj %*% Betaj + rnorm(n0))
}) %>% unlist

exdata <- data.frame(
  y = y,
  X,
  group = area
)[sample(1:n, n),]
rownames(exdata) <- NULL

usethis::use_data(exdata, overwrite = TRUE)

# adj <- read.table("inst/adj.txt", header = T)
# usethis::use_data(adj, overwrite = TRUE)
