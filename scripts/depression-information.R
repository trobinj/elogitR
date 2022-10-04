library(numDeriv)
library(dplyr)
library(ggplot2)
library(tikzDevice)

setwd("/home/trjohns/Dropbox/Research/elogit/manuscript")

n <- 1000
m <- 5
k <- 4
v <- 3

x <- cbind(1, -diag(m))

theta <- c(alph, gamm, beta)

elprob <- function(z, theta, yvec, xvec, Kmax) {
  exp(elprofile(z, theta, yvec, xvec, Kmax))
}

elinfo <- function(z, theta, yvec, xvec, Kmax) {
  p0 <- elprob(z, theta, yvec, xvec, Kmax)
  p1 <- numDeriv::grad(elprob, x = z, yvec = yvec, xvec = xvec, theta = theta, Kmax = Kmax)
  p2 <- numDeriv::hessian(elprob, x = z, yvec = yvec, xvec = xvec, theta = theta, Kmax = Kmax)
  as.numeric(p1^2/p0 - p2)
}

flprob <- function(z, theta, yvec, xvec) {
  exp(flprofile(z, theta, yvec, xvec))
}

flinfo <- function(z, theta, yvec, xvec) {
  p0 <- flprob(z, theta, yvec, xvec)
  p1 <- numDeriv::grad(flprob, x = z, yvec = yvec, xvec = xvec, theta = theta)
  p2 <- numDeriv::hessian(flprob, x = z, yvec = yvec, xvec = xvec, theta = theta)
  as.numeric(p1^2/p0 - p2)
}

zeta <- seq(-3, 3, by = 0.1/4)

##########################

theta <- estimates[,1]

tau <- c(rep(0, k-2), theta[1:(k-2)])
sigma <- theta[k-1]
delta <- theta[k:length(theta)]

alph <- -cumsum(tau)
gamm <- 0:(k-1)
beta <- c(sigma, delta)

theta <- c(alph, gamm, beta)

Kmax <- 1
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[v,], Kmax)
  }
}
inf1 <- apply(inf, 1, sum)

##########################

theta <- estimates[,2]

tau <- c(rep(0, k-2), theta[1:(k-2)])
sigma <- theta[k-1]
delta <- theta[k:length(theta)]

alph <- -cumsum(tau)
gamm <- 0:(k-1)
beta <- c(sigma, delta)

theta <- c(alph, gamm, beta)

Kmax <- 2
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[v,], Kmax)
  }
}
inf2 <- apply(inf, 1, sum)

##########################

theta <- estimates[,3]

tau <- c(rep(0, k-2), theta[1:(k-2)])
sigma <- theta[k-1]
delta <- theta[k:length(theta)]

alph <- -cumsum(tau)
gamm <- 0:(k-1)
beta <- c(sigma, delta)

theta <- c(alph, gamm, beta)

Kmax <- 4
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[v,], Kmax)
  }
}
infK <- apply(inf, 1, sum)

##########################

theta <- estimates[,4]

tau <- c(rep(0, k-2), theta[1:(k-2)])
sigma <- theta[k-1]
delta <- theta[k:length(theta)]

alph <- -cumsum(tau)
gamm <- 0:(k-1)
beta <- c(sigma, delta)

theta <- c(alph, gamm, beta)

y <- flpatterns(K)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- flinfo(zeta[i], theta, y[j,], x[v,])
  }
}
inf0 <- apply(inf, 1, sum)

d2 <- data.frame(inf = c(inf0, inf1, inf2, infK),
  zeta = rep(zeta, 4), K = rep(c(0,1,2,4), each = length(zeta)))

prnt <- FALSE

if (prnt) tikz("depressionscale-information.tex", width = 6, height = 4)

p <- ggplot(d2, aes(x = zeta, y = inf, group = K)) + theme_classic() +
  geom_line() + scale_x_continuous(breaks = -3:3) + 
  labs(x = "$\\zeta$", y = "Information") + 
  theme(axis.text = element_text(color = "black")) + 
  theme(axis.title.y = element_text(size = 10))

p <- p + annotate("segment", x = -1.5, xend = -0.62, y = 1.75, yend = 1.5,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = -1.5, y = 1.75, hjust = "right", label = "B/W")

p <- p + annotate("segment", x = 1.5, xend = 0.42, y = 2.25, yend = 2,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = 1.5, y = 2.25, hjust = "left", label = "$K' = 3$")

p <- p + annotate("segment", x = 1.5, xend = 0.4, y = 1.75, yend = 1.51,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = 1.5, y = 1.75, hjust = "left", label = "$K' = 2$")

p <- p + annotate("segment", x = 1.5, xend = 0.75, y = 1.25, yend = 0.85,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = 1.5, y = 1.25, hjust = "left", label = "$K' = 1$")

plot(p)

if (prnt) graphics.off()
