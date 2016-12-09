## R course, Dong 2016 
## Solution of Exercise 4


library("tidyverse")
easy0 <- read_csv("http://nielsrhansen.github.io/Dong/easysmooth.txt")
easy <- easy0

## Plotting with ggplot

ggplot(easy, aes(x = x, y = y)) + 
  geom_point()

## Alternatively with 'qplot'

qplot(x, y, data = easy)

## First implementation

n <- nrow(easy)
s <- rep(NA, length(easy$y))
m <- 10
for(i in (m + 1):(n - m)) 
  s[i] <- mean(easy$y[(i - m):(i + m)])

easy$smooth <- s

## This could have been done with 'mutate' as well
## mutate(easy, smooth = s)

qplot(x, y, data = easy) +
  geom_line(aes(y = smooth), 
            color = "red", 
            size = 1)

##
## Reimplementation using a function
##

runMean <- function(y, x, m = 10, fun = mean) {
  if(!missing(x))
    y <- y[order(x)]
  n <- length(y)
  s <- rep(NA, length(y))
  for(i in (m + 1):(n - m)) 
    s[i] <- fun(y[(i - m):(i + m)])
  s
}

easy <- easy0
easy <- mutate(easy, smooth = runMean(y),
  smoothMed = runMean(y, fun = median))

qplot(x, y, data = easy) +
  geom_line(aes(y = smooth), 
              color = "red", 
              size = 1) + 
  geom_line(aes(y = smoothMed), 
              color = "blue", 
              size = 1)




