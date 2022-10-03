library(ggplot2)
library(dplyr)
library(tidyr)
library(progress)

setwd(this.path::this.dir())

n <- 200

delta <- c(-2.5, -1.5, -0.5)
tau <- c(0, 0, 1, 2, 3)
sigma <- 1

m <- length(delta)
K <- length(tau)

alph <- -cumsum(tau)
gamm <- 0:(K-1)
beta <- c(sigma, delta)

d <- expand.grid(x = seq(-3, 3, length = n), k = 1:K, j = 1:m) %>%
  mutate(p = exp(alph[k] + gamm[k] * (sigma*x - delta[j]))) %>% 
  group_by(x, j) %>% mutate(p = p / sum(p)) 

# p <- ggplot(d, aes(x = x, y = p, color = factor(k))) + theme_minimal() + 
#    geom_line() + scale_x_continuous(breaks = -6:6) + facet_wrap(~ j, nrow = 1)
# plot(p)

p <- 21 # quadrature points

reps <- 10000

theta <- c(tau[-c(1:2)], beta)

d1 <- matrix(NA, reps, length(theta))
d2 <- d1
d3 <- d1
dK <- d1
d0 <- d1

set.seed(123)

pb <- progress_bar$new(
  format = "[:bar] :percent eta: :eta",
  total = reps, clear = FALSE)

for (i in 1:reps) {

  pb$tick()

  x <- cbind(rep(rnorm(n, 0, 1), each = m), kronecker(rep(1, n), -diag(m)))  
  y <- dat_samp(c(alph, gamm, beta), x, K)
  
  d1[i,] <- elrsm(theta, y, m, p, 1) - theta
  d2[i,] <- elrsm(theta, y, m, p, 2) - theta
  d3[i,] <- elrsm(theta, y, m, p, 3) - theta
  dK[i,] <- elrsm(theta, y, m, p, K) - theta
  d0[i,] <- flrsm(theta, y, m, p) - theta
}

d <- as.data.frame(rbind(d1, d2, d3, dK, d0))
names(d) <- c("t3","t4","t5","sigma","d1","d2","d3")

write.csv(d, file = "elrsm-simdata.csv", row.names = FALSE)




