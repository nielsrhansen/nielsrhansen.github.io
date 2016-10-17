library(ggplot2)
library(gridExtra)
library(Matrix)

mixedmodel <- function(data, par, sigma = c(2, 4)) {
  J <- table(data$lab)
  p <- length(par)
  n <- sum(J)
  ybar <- tapply(data$y, data$lab, mean)
  ymean <- sum(J * ybar) / n
  V  <-  1 / (J / sigma[1]^2 + 1 / sigma[2]^2)
  sqrtV <- sqrt(V)
  JVsig <- J * V / sigma[1]^2
  J <- J / n
  step <- function(z, u) {
    mu <- u[1] - sum(J * z[-1])
    alpha <- u[-1] - JVsig * mu
    c(mu, alpha)
  }
  structure(
    list(par = par,
         p = p,
         stepmu = c(ymean, JVsig * ybar),
         stepsigma = c(sigma[1] / sqrt(n), sqrtV),
         step = step
    ), class = "mixedmodel")
}

gibbs <- function(x, t, ...) 
  UseMethod("gibbs")

gibbs.mixedmodel <- function(x, t = 1000, ...) {
  par <- matrix(0, length(x$par), t)
  par[, 1] <- x$par
  U <- matrix(rnorm(x$p * t, x$stepmu, x$stepsigma), x$p, t)
  for(i in 2:t)
    par[, i] <- x$step(par[, i - 1], U[, i])
  x$gibbspar <- par
  x
}

plot.mixedmodel <- function(x, ...) {
  t <- 1:ncol(x$gibbspar)
  temp <- data.frame(t(rbind(t, x$gibbspar)))
  temp <- reshape2::melt(temp, id.var = "t")
  qplot(t, value, data = temp, color = variable, geom = "line", size = I(2)) +
    ylab("") + geom_point(size = 4)
}

summary.mixedmodel <- function(object, burnin = 100, plot = TRUE, showcor = TRUE, ...) {
  temp <- t(object$gibbspar)
  temp <- temp[, -1] + temp[, 1] 
  muhat <- colMeans(temp)
  varhat <- var(temp)
  sigma <- sqrt(diag(varhat))
  t <- 1:ncol(object$gibbspar)
  temp <- data.frame(cbind(t, temp))
  varnames <- colnames(temp)[-1]
  if (plot) {
    p1 <- qplot(variable, value, 
                data = reshape2::melt(temp[-seq(1, burnin), ], id.var = "t"), 
                geom = "violin", fill = I("lightblue"), alpha = I(0.5)) + 
      geom_point(data = data.frame(x = varnames, y = muhat), 
                 aes(x = x, y = y), color = "red", size = 4) + 
      geom_errorbar(data = data.frame(x = varnames, 
                                      min = muhat - 1.96 * sigma, 
                                      max = muhat + 1.96 * sigma, 
                                      y = muhat),
                    aes(x = x, min = min, max = max, y = y)) +
      ylab("posterior density") + xlab("mean")
    if (showcor) {
      p2 <- Matrix::image(Matrix(varhat), useAbs = FALSE)
      gridExtra::grid.arrange(p1, p2, ncol = 2)
    } else {
      p1
    }
  } 
  out <- cbind(muhat, sigma)
  rownames(out) <- varnames
  corind <- c()
  if (showcor) { 
    out <- cbind(out, varhat / (sigma %*% t(sigma)))
    dimnames(out) <- list(varnames, c("Mean est.", "Stand. err.", varnames))
    corind <- seq(3, ncol(varhat))
  }
  printCoefmat(out, digits = 3, tst.ind = corind, dig.tst = 2)
}
  



