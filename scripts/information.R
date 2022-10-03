library(numDeriv)
library(dplyr)
library(ggplot2)
library(tikzDevice)

setwd(this.path::this.dir())

n <- 1000

delta <- c(-0.5, -1.5, -2.5)
tau <- c(0, 0, 1, 2, 3) 
sigma <- 1

m <- length(delta)
K <- length(tau)

alph <- -cumsum(tau)
gamm <- 0:(K-1)
beta <- c(sigma, delta)

# d <- expand.grid(x = seq(-3, 3, length = n), k = 1:K, j = 1:m) %>%
#   mutate(p = exp(alph[k] + gamm[k] * (sigma*x - delta[j]))) %>% group_by(x, j) %>% 
#   mutate(p = p / sum(p)) 
# 
# p <- ggplot(d, aes(x = x, y = p, color = factor(k))) + theme_minimal() + 
#   geom_line() + scale_x_continuous(breaks = -6:6) + facet_wrap(~ j, nrow = 1)
# plot(p)

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

Kmax <- K
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[2,], Kmax)
  }
}
infK <- apply(inf, 1, sum)

Kmax <- 1
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[2,], Kmax)
  }
}
inf1 <- apply(inf, 1, sum)

Kmax <- 2
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[2,], Kmax)
  }
}
inf2 <- apply(inf, 1, sum)

Kmax <- 3
y <- elpatterns(K, Kmax)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[2,], Kmax)
  }
}
inf3 <- apply(inf, 1, sum)

y <- flpatterns(K)
inf <- matrix(NA, length(zeta), nrow(y))
for (i in 1:length(zeta)) {
  for (j in 1:nrow(y)) {
    inf[i,j] <- flinfo(zeta[i], theta, y[j,], x[2,])
  }
}
inf0 <- apply(inf, 1, sum)

d2 <- data.frame(inf = c(inf0, inf1, inf2, inf3, infK), 
  zeta = rep(zeta, 5), K = rep(c(0,1,2,3,5), each = length(zeta)))

prnt <- TRUE

if (prnt) tikz("../manuscript/information.tex", width = 6, height = 4)

p <- ggplot(d2, aes(x = zeta, y = inf, group = K)) + theme_classic() +
  geom_line() + scale_x_continuous(breaks = -3:3) + 
  labs(x = "$\\zeta$", y = "Information") + 
  theme(axis.text = element_text(color = "black")) + 
  theme(axis.title.y = element_text(size = 10))

p <- p + annotate("segment", x = 1.5, xend = 0.51, y = 4, yend = 3.6,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + 
  annotate("text", x = 1.5, y = 4, hjust = "left", label = "$K' = 4$")

p <- p + annotate("segment", x = 1.5, xend = 0.65, y = 3.5, yend = 2.63,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + 
  annotate("text", x = 1.5, y = 3.5, hjust = "left", label = "$K' = 3$")

p <- p + annotate("segment", x = 1.5, xend = 0.5, y = 3, yend = 1.51,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + 
  annotate("text", x = 1.5, y = 3, hjust = "left", label = "$K' = 2$")

p <- p + annotate("segment", x = 1.5, xend = 0.75, y = 2.5, yend = 0.88,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + 
  annotate("text", x = 1.5, y = 2.5, hjust = "left", label = "$K' = 1$")

p <- p + annotate("segment", x = -1, xend = -0.42, y = 3.5, yend = 3,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) + 
  annotate("text", x = -1, y = 3.5, hjust = "right", label = "B/W")

plot(p)

if (prnt) graphics.off()

# Kmax <- K
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[1,], Kmax)
#   }
# }
# infK <- apply(inf, 1, sum)
# 
# Kmax <- 1
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[1,], Kmax)
#   }
# }
# inf1 <- apply(inf, 1, sum)
# 
# Kmax <- 2
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[1,], Kmax)
#   }
# }
# inf2 <- apply(inf, 1, sum)
# 
# Kmax <- 3
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[1,], Kmax)
#   }
# }
# inf3 <- apply(inf, 1, sum)
# 
# y <- flpatterns(K)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- flinfo(zeta[i], theta, y[j,], x[1,])
#   }
# }
# inf0 <- apply(inf, 1, sum)
# 
# d1 <- data.frame(inf = c(inf0, inf1, inf2, inf3, infK), 
#   zeta = rep(zeta, 5), K = rep(c(0,1,2,3,5), each = length(zeta)))
# 
# Kmax <- K
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[3,], Kmax)
#   }
# }
# infK <- apply(inf, 1, sum)
# 
# Kmax <- 1
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[3,], Kmax)
#   }
# }
# inf1 <- apply(inf, 1, sum)
# 
# Kmax <- 2
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[3,], Kmax)
#   }
# }
# inf2 <- apply(inf, 1, sum)
# 
# Kmax <- 3
# y <- elpatterns(K, Kmax)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[3,], Kmax)
#   }
# }
# inf3 <- apply(inf, 1, sum)
# 
# y <- flpatterns(K)
# inf <- matrix(NA, length(zeta), nrow(y))
# for (i in 1:length(zeta)) {
#   for (j in 1:nrow(y)) {
#     inf[i,j] <- flinfo(zeta[i], theta, y[j,], x[3,])
#   }
# }
# inf0 <- apply(inf, 1, sum)
# 
# d3 <- data.frame(inf = c(inf0, inf1, inf2, inf3, infK), 
#   zeta = rep(zeta, 5), K = rep(c(0,1,2,3,5), each = length(zeta)))
# 
# d <- rbind(d1,d2,d3) %>% group_by(zeta,K) %>% summarize(inf = sum(inf))
# 
# p <- ggplot(d, aes(x = zeta, y = inf, group = K, color = factor(K))) + theme_classic() +
#   geom_line() + scale_x_continuous(breaks = -3:3) + 
#   labs(x = "$\\zeta$", y = "Information") + 
#   theme(axis.text = element_text(color = "black"))
# plot(p)
# 
