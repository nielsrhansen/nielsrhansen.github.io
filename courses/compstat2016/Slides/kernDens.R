kernDens <- function (x, h, n = 512, 
                      kern = function(x, h) dnorm(x, sd = h)) {
  rg <- range(x)
  ## xx equivalent to spec. in density()
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = n)
  y <- numeric(n)
  for (i in seq_along(xx))
    y[i] <- mean(kern(xx[i] - x, h))
  list(x = xx, y = y, h = h)
}

kernDens2 <- function (x, h, n = 512, 
                       kern = function(x, h) dnorm(x, sd = h)) {
  ## Silverman's rule if h is missing
  if(missing(h))
    h <- 1.06 * length(x)^(-0.2) * sd(x)
  rg <- range(x)
  xx <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = n)
  dif <- outer(x, xx, "-")
  kerndif <- kern(dif, h)
  y <- colMeans(kerndif)
  list(x = xx, y = y, h = h)
}

kernbin <- function(x, lo, hi, n) {
  w <- numeric(n)
  delta <- (hi - lo) / (n - 1)
  for(i in seq_along(x)) {
    ii <- floor((x[i] - lo) / delta + 0.5) + 1
    w[ii] <- w[ii] + 1
  }
  w / sum(w)
}

## This implementation assumes a symmetric kernel! 
## It is possible to make an implementation that does
## not rely on symmetry, but it is a little more 
## complicated.
kernDens3 <- function(x, h, n = 512) {
  kern = function(x, h) dnorm(x, sd = h)
  ## Silverman's rule if h is missing
  if(missing(h))
    h <- 1.06 * length(x)^(-0.2) * sd(x)
  rg <- range(x) + c(- 3* h, 3 * h)
  xx <- seq(rg[1], rg[2], length.out = n)
  weights <- kernbin(x, rg[1], rg[2], n)
  kerneval <- kern(xx - xx[1], h)
  kerndif <- matrix(kerneval[toeplitz(1:n)], n, n)
  y <- colSums(weights * kerndif)
  list(x = xx, y = y, h = h)
}