library(ggplot2)
library(dplyr)
library(tidyr)
library(progress)

setwd(this.path::this.dir())

K <- 3
n <- 200

alph <- -cumsum(c(0, c(-1.5, 1.5)))
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

reps <- 10000

d1 <- matrix(NA, reps, length(theta))
dK <- d1

set.seed(123)

pb <- pb <- progress_bar$new(
   format = "[:bar] :percent eta: :eta",
   total = reps, clear = FALSE)


for (i in 1:reps) {
   
   pb$tick()

   y <- dat_samp(c(alph, gamm, beta), x, K)
   
   theta <- c(alph[-1], gamm[-c(1,K)], beta)
   
   d1[i,] <- elreg(theta, y, x, 1) - theta
   dK[i,] <- elreg(theta, y, x, K) - theta
}

d <- as.data.frame(rbind(d1, dK))
names(d) <- c("a2","a3","g2","beta")

write.csv(d, file = "elreg-3-simdata.csv", row.names = FALSE)


   

