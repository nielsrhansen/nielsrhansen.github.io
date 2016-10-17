gibbsstep <- function(mu, alpha, ybar, J, 
                      sigmaeps = 2, sigmaalph = 4) {
  p <- length(alpha)
  n <- sum(J)
  ymean <- sum(J * ybar) / n
  V <- 1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  mu <- rnorm(1, ymean - sum(J * alpha) / n, sigmaeps / sqrt(n))
  alpha <- rnorm(p, J * V * (ybar - mu) / sigmaeps^2, sqrt(V))
  list(mu = mu, alpha = alpha)
}

gibbs <- function(mu0, alpha0, data, t = 1000, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  ybar <- tapply(data$y, data$lab, mean)
  for(i in 2:t) {
    tmp <- gibbsstep(mu[i - 1], alpha[, i - 1], ybar, J, ...)
    mu[i] <- tmp$mu
    alpha[, i] <- tmp$alpha
  }
  list(mu = mu, alpha = alpha) 
}

gibbsfast <- function(mu0, alpha0, data, t = 1000, sigmaeps = 2, 
                      sigmaalph = 4, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  p <- length(alpha0)
  n <- sum(J)
  ybar <- tapply(data$y, data$lab, mean)
  ymean <- sum(J * ybar) / n
  V  <-  1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  sigmaepsn <-  sigmaeps / sqrt(n)
  sqrtV <- sqrt(V)
  JVsig <- J * V / sigmaeps^2
  for(i in 2:t) {
    mu[i] <- rnorm(1, ymean - sum(J * alpha[, i - 1]) / n, sigmaepsn)
    alpha[, i] <- rnorm(p, JVsig * (ybar - mu[i]), sqrtV)
  }
  list(mu = mu, alpha = alpha) 
}

gibbsfastest <- function(mu0, alpha0, data, t = 1000, sigmaeps = 2, 
                         sigmaalph = 4, ...) {
  mu <- numeric(t)
  alpha <- matrix(0, length(alpha0), t)
  mu[1] <- mu0
  alpha[, 1] <- alpha0
  J <- table(data$lab)
  p <- length(alpha0)
  n <- sum(J)
  ybar <- tapply(data$y, data$lab, mean)
  ymean <- sum(J * ybar) / n
  V  <-  1 / (J / sigmaeps^2 + 1 / sigmaalph^2)
  sigmaepsn <-  sigmaeps / sqrt(n)
  sqrtV <- sqrt(V)
  JVsig <- J * V / sigmaeps^2
  J <- J / n
  W <- rnorm(t, ymean, sigmaepsn)
  U <- matrix(rnorm(t * p, JVsig * ybar, sqrtV), p, t)
  for(i in 2:t) {
    mu[i] <- W[i] - sum(J * alpha[, i - 1])
    alpha[, i] <- U[, i] - JVsig * mu[i]
  }
  list(mu = mu, alpha = alpha) 
}



