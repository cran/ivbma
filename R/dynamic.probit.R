probit.bma.update.M <- function(theta,D)
  {
    ##------ Keep book --------
    n <- length(D$Y)
    X <- D$X
    Z <- D$Z
    p.X <- dim(D$X)[2]
    M <- theta$M
    ##------------------------

    ##----- Flip a variable -----
    w <- sample(1:p.X,1)
    M.new <- theta$M
    M.new[w] <- 1 - M.new[w]
    if(sum(M.new) == 0) return(M)
    ##-----------------------------
    
    ##---- New Score ------------
    p.2 <- sum(M.new)
    X.2 <- X[,(1:p.X)[M.new == 1],drop=FALSE]
    Omega.2 <- diag(p.2) + t(X.2) %*% X.2
    lambda.2 <- t(Z) %*% X.2 %*% solve(Omega.2)
    score.2 <- 0.5 * lambda.2 %*% Omega.2 %*% t(lambda.2) - 0.5 * mydet(Omega.2)
    ##---------------------------

    ##---- Old Score ------------
    p.1 <- sum(M)
    X.1 <- X[,(1:p.X)[M == 1],drop=FALSE]
    Omega.1 <- diag(p.1) + t(X.1) %*% X.1
    lambda.1 <- t(Z) %*% X.1 %*% solve(Omega.1)
    score.1 <- 0.5 * lambda.1 %*% Omega.1 %*% t(lambda.1) - 0.5 * mydet(Omega.1)
    ##---------------------------

    ##---- Decide ---------------
    alpha <- score.2 - score.1
    if(log(runif(1)) < alpha)
      {
        M <- M.new
        Omega.1 <- Omega.2
        lambda.1 <- lambda.2
      }
    ##---------------------------

    lambda <- rep(0,p.X)
    lambda[ (1:p.X)[M==1]] <- rmvnorm.precision(lambda.1,Omega.1)

    return(list(M = M, lambda = lambda))
  }

probit.bma.updateZ <- function (theta,D)
{
  means <- D$X %*% theta$beta
  theta$Z[D$Y == 0] <- qnorm(runif(sum(D$Y == 0), 0, 0.5), means,1)
  theta$Z[D$Y == 1] <- qnorm(runif(sum(D$Y == 1), 0.5,1),means,1)
  return(theta$Z)

}

probit.bma.init <- function(D,full)
  {
    theta <- NULL
    n <- length(D$Y)
    theta$Z <- NULL
    theta$Z[D$Y == 0] <- qnorm(runif(sum(D$Y == 0), 0, 0.5))
    theta$Z[D$Y == 1] <- qnorm(runif(sum(D$Y == 1), 0.5,1))
    theta$beta <- solve(t(D$X) %*% D$X) %*% t(D$X) %*% D$Z
    p.X <- dim(D$X)[2]
    theta$M <- rbinom(p.X,1,.5)
    theta$beta[theta$M == 0] <- 0
    return(theta)
  }

probit.bma.sample <- function(theta,D)
  {
    theta$Z <- probit.bma.updateZ(theta,D)
    l <- probit.bma.update.M(theta,D)
    theta$beta <- l$lambda
    theta$M <- l$M
    return(theta)
  }

probit.bma.results.init <- function(D,odens)
  {
      p.X <- dim(D$X)[2]

      ##-------- Information to be returned ----------
      results <- NULL
      results$beta  <- matrix(0, odens, p.X)
      results$beta.bar <- rep(0,p.X)
      results$M <- matrix(0,odens,p.X)
      results$M.bar <- rep(0, p.X)
      ##----------------------------------------------
      
      return(results)
  }

probit.bma <- function(Y,X,s=1e3,b = round(s/10),
                       full = FALSE,odens = min(c(5e3,s-b)),
                       print.every = round(s/10))
  {
    D <- NULL
    D$Y <- Y
    D$X <- X
    theta <- probit.bma.init(D,full)

    results <- probit.bma.results.init(D,odens)
    which.save <- round(seq(b + 1, s, length = odens))
    save.loc <- 1
    next.save <- which.save[save.loc]

    for(i in 1:s)
      {
        theta <- probit.bma.sample(theta,D)
        ##---------- Record ----------------------------
        if(i == next.save)
          {
            ##--------------- Save -----------------------------------
            results$beta[save.loc,] <- theta$beta
            results$M[save.loc,] <- theta$M
            save.loc <- save.loc + 1
            next.save <- which.save[save.loc]
            ##--------------------------------------------------------
          }
        if(i > b)
          {
            results$beta.bar <- results$beta.bar + theta$beta/ (s - b)
            results$M.bar <- results$M.bar + theta$M / (s - b)
          }
        ##------------------------------------------------------
      }
    
    return(results)
  }
