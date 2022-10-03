library(mosaicData)
library(dplyr)
library(tidyr)
library(eRm)
library(ggplot2)
library(numDeriv)
library(dplyr)
library(tikzDevice)

setwd(this.path::this.dir())

# Items:
# 1. Bothered by things not generally bothered by
# 2. Had trouble keeping mind on what doing
# 3. I felt everything I did was an effort
# 4. My sleep was restless
# 5. I could not get going

d <- HELPfull %>% select(F1A,F1E,F1G,F1K,F1T) %>% na.omit()

m <- ncol(d)
n <- nrow(d)

d <- d %>% mutate(id = 1:n()) %>%
  pivot_longer(cols = -id, names_to = "item", values_to = "y")

p <- 21
k <- length(unique(d$y))
y <- matrix(0, nrow(d), k)
y[,1] <- d$y

if (min(y[,1]) == 0) {
  y[,1] <- y[,1] + 1
}

theta <- c(1:(k-2), 1, rep(-1, m))

theta <- elrsm(theta, y, m, p, 1, TRUE)$estimate

tau <- c(rep(0, k-2), theta[1:(k-2)])
sigma <- theta[k-1]
delta <- theta[k:length(theta)]

alph <- -cumsum(tau)
gamm <- 0:(k-1)
beta <- c(sigma, theta[k:length(theta)])

set.seed(123)
x <- cbind(rep(rnorm(n, 0, 1), each = m), kronecker(rep(1, n), -diag(m)))
y <- stm_samp(c(alph, gamm, beta), x, k)

tmp1 <- elrsm(theta, y, m, p, 1, TRUE)
tmp2 <- elrsm(theta, y, m, p, 2, TRUE)
tmpK <- elrsm(theta, y, m, p, 3, TRUE)
tmp0 <- flrsm(theta, y, m, p, TRUE)

estimates <- cbind(tmp1$estimate, tmp2$estimate, tmpK$estimate, tmp0$estimate)
stderrors <- cbind(tmp1$stderror, tmp2$stderror, tmpK$stderror, tmp0$stderror)

# compute percent reduction in standard errors
tmp <- stderrors
for (i in 1:ncol(tmp)) {
  tmp[,i] <- (tmp[,i] - stderrors[,1])/stderrors[,1] * 100
}
print(tmp)

tbl <- data.frame(
  e1 = estimates[,1], s1 = stderrors[,1],
  e2 = estimates[,2], s2 = stderrors[,2],
  e3 = estimates[,3], s3 = stderrors[,3],
  e4 = estimates[,4], s4 = stderrors[,4])

for (i in 1:ncol(tbl)) {
  tbl[,i] <- format(round(tbl[,i], 3), nsmall = 3)
}
for (i in c(2,4,6,8)) {
  tbl[,i] <- paste("(", tbl[,i], ")", sep = "")
}

tbl <- tbl[c(1,2,4:8,3),]

param <- c("$\\tau_3$", "$\\tau_4$", "$\\delta_1$",
  "$\\delta_2$", "$\\delta_3$", "$\\delta_4$", "$\\delta_5$", "$\\sigma$")

tbl <- cbind(param, tbl)

textable <- function(x, file) {
  for (i in 1:nrow(x)) {
    if (i == nrow(x)) {
      eolstring <- " \\tabularnewline\\bottomrule"
    } else {
      eolstring <- " \\\\\n"
    }
    write.table(x[i,], file, append = i > 1, row.names = FALSE,
      col.names = FALSE, quote = FALSE, sep = " & ", na = "", eol = eolstring)
  }
}
textable(tbl, "../manuscript/depression-table.tex")

n <- 1000
m <- 5
k <- 4
K <- k
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
inf <- matrix(0, length(zeta), nrow(y))
inf1 <- rep(0, length(zeta))
for (v in 1:5) {
  for (i in 1:length(zeta)) {
    for (j in 1:nrow(y)) {
      inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[v,], Kmax)
    }
  }
  inf1 <- inf1 + apply(inf, 1, sum)
}

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
inf <- matrix(0, length(zeta), nrow(y))
inf2 <- rep(0, length(zeta))
for (v in 1:5) {
  for (i in 1:length(zeta)) {
    for (j in 1:nrow(y)) {
      inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[v,], Kmax)
    }
  }
  inf2 <- inf2 + apply(inf, 1, sum)
}

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
inf <- matrix(0, length(zeta), nrow(y))
infK <- rep(0, length(zeta))
for (v in 1:5) {
  for (i in 1:length(zeta)) {
    for (j in 1:nrow(y)) {
      inf[i,j] <- elinfo(zeta[i], theta, y[j,], x[v,], Kmax)
    }
  }
  infK <- infK + apply(inf, 1, sum)
}

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
inf <- matrix(0, length(zeta), nrow(y))
inf0 <- rep(0, length(zeta))
for (v in 1:5) {
  for (i in 1:length(zeta)) {
    for (j in 1:nrow(y)) {
      inf[i,j] <- flinfo(zeta[i], theta, y[j,], x[v,])
    }
  }
  inf0 <- inf0 + apply(inf, 1, sum)
}

d2 <- data.frame(inf = c(inf0, inf1, inf2, infK),
  zeta = rep(zeta, 4), K = rep(c(0,1,2,4), each = length(zeta)))

prnt <- TRUE

if (prnt) tikz("../manuscript/depression-information.tex", width = 6, height = 4)

p <- ggplot(d2, aes(x = zeta, y = inf, group = K)) + theme_classic() +
  geom_line() + scale_x_continuous(breaks = -3:3) + 
  labs(x = "$\\zeta$", y = "Information") + 
  theme(axis.text = element_text(color = "black")) + 
  theme(axis.title.y = element_text(size = 10))

y <- 8
x <- -1.5
p <- p + annotate("segment", x = x, xend = -0.65, y = y, yend = 6,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = x, y = y, hjust = "right", label = "B/W")

y <- 10
x <- 1.5
p <- p + annotate("segment", x = x, xend = 0.78, y = y, yend = 8,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = x, y = y, hjust = "left", label = "$K' = 3$")

y <- 8
x <- 1.5
p <- p + annotate("segment", x = x, xend = 0.75, y = y, yend = 6.5,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = x, y = y, hjust = "left", label = "$K' = 2$")

y <- 6
x <- 1.5
p <- p + annotate("segment", x = x, xend = 0.97, y = y, yend = 4,
  arrow = arrow(length = unit(0.1, "cm"), type = "closed")) +
  annotate("text", x = x, y = y, hjust = "left", label = "$K' = 1$")

plot(p)

if (prnt) graphics.off()
