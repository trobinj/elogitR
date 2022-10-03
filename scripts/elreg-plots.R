library(ggplot2)
library(dplyr)
library(tikzDevice)

setwd(this.path::this.dir())

K <- 3
n <- 200

alph <- -cumsum(c(0, c(-1.5, 1.5)))
gamm <-  seq(0, K - 1, by = 1) / (K - 1)
beta <- K - 1

d1 <- expand.grid(x = seq(-3, 3, length = n), k = 1:K) %>%
   mutate(p = exp(alph[k] + gamm[k] * beta * x)) %>% group_by(x) %>% 
   mutate(p = p / sum(p)) %>% mutate(K = 3)

K <- 4
n <- 200

alph <- -cumsum(c(0, c(-2, 0, 2)))
gamm <-  seq(0, K - 1, by = 1) / (K - 1)
beta <- K - 1

d2 <- expand.grid(x = seq(-3, 3, length = n), k = 1:K) %>%
   mutate(p = exp(alph[k] + gamm[k] * beta * x)) %>% group_by(x) %>% 
   mutate(p = p / sum(p)) %>% mutate(K = 4)

K <- 5
n <- 200

alph <- -cumsum(c(0, c(-2.25, -0.75, 0.75, 2.25)))
gamm <-  seq(0, K - 1, by = 1) / (K - 1)
beta <- K - 1

d3 <- expand.grid(x = seq(-3, 3, length = n), k = 1:K) %>%
   mutate(p = exp(alph[k] + gamm[k] * beta * x)) %>% group_by(x) %>% 
   mutate(p = p / sum(p)) %>% mutate(K = 5)

K <- 7
n <- 200

alph <- -cumsum(c(0, c(-2.5, -1.5, -0.5, 0.5, 1.5, 2.5)))
gamm <-  seq(0, K - 1, by = 1) / (K - 1)
beta <- K - 1

d4 <- expand.grid(x = seq(-3, 3, length = n), k = 1:K) %>%
   mutate(p = exp(alph[k] + gamm[k] * beta * x)) %>% group_by(x) %>% 
   mutate(p = p / sum(p)) %>% mutate(K = 7)

d <- rbind(d1, d2, d3, d4) %>%
   mutate(K = factor(K, levels = c(3,4,5,7), 
      labels = paste("$K = ", c(3,4,5,7), "$", sep = "")))

prnt <- FALSE

if (prnt) tikz("../manuscript/elreg-plot.tex", width = 6, height = 3)

p <- ggplot(d, aes(x = x, y = p)) + theme_minimal() + geom_line(aes(group = k)) + 
   facet_wrap(~ K, ncol = 4) + labs(x = "$x$", y = "Probability") + ylim(0,1) + 
   theme(axis.text = element_text(color = "black", size = 6)) + 
   scale_x_continuous(breaks = -3:3) + theme(axis.title.y = element_text(size = 10))
plot(p)

if (prnt) graphics.off()



