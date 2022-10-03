library(ggplot2)
library(dplyr)
library(tidyr)

setwd(this.path::this.dir())

K <- 4
n <- 200

alph <- -cumsum(c(0, c(-2, 0, 2)))
gamm <-  seq(0, K - 1, by = 1) / (K - 1)
beta <- K - 1

d <- expand.grid(x = seq(-3, 3, length = n), k = 1:K) %>%
   mutate(p = exp(alph[k] + gamm[k] * beta * x)) %>% group_by(x) %>% 
   mutate(p = p / sum(p)) 

p <- ggplot(d, aes(x = x, y = p, color = factor(k))) + theme_minimal() + 
   geom_line() + scale_x_continuous(breaks = -6:6)
plot(p)

theta <- c(alph[-1], gamm[-c(1,K)], beta)

x <- matrix(seq(-3, 3, length = n), n, 1)

reps <- 100

d1 <- matrix(NA, reps, length(theta))
d2 <- d1
dK <- d1
d0 <- d1

set.seed(123)

pb <- progress_bar$new(total = reps)

for (i in 1:reps) {
   
   pb$tick()

   y <- dat_samp(c(alph, gamm, beta), x, K)
   
   theta <- c(alph[-1], gamm[-c(1,K)], beta)
   
   d1[i,] <- elreg(theta, y, x, 1) - theta
   d2[i,] <- elreg(theta, y, x, 2) - theta
   dK[i,] <- elreg(theta, y, x, K) - theta
   d0[i,] <- flreg(theta, y, x) - theta
}

d <- as.data.frame(rbind(d1, d2, dK, d0))
names(d) <- c("a2","a3","a4","g2","g3","beta")

#write.csv(d, file = "elreg-4-simdata.csv", row.names = FALSE)


   

