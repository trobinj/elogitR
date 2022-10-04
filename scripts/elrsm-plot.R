library(ggplot2)
library(dplyr)
library(tidyr)
library(tikzDevice)

setwd("/home/trjohns/Dropbox/Research/elogit/manuscript")

n <- 200

delta <- c(-2.5, -1.5, -0.5)
tau <- c(0, 0, 1, 2, 3)
sigma <- 1

m <- length(delta)
K <- length(tau)

alph <- -cumsum(tau)
gamm <- 0:(K-1)
beta <- c(sigma, delta)

prnt <- TRUE

d <- expand.grid(x = seq(-3, 3, length = n), k = 1:K, j = 1:m) %>%
   mutate(p = exp(alph[k] + gamm[k] * (sigma*x - delta[j]))) %>% group_by(x, j) %>% 
   mutate(p = p / sum(p)) %>%
   mutate(j = factor(j, levels = 1:3, labels = paste("$j = ", 1:3, "$", sep = "")))

if (prnt) tikz("elrsm-plot.tex", width = 6, height = 3)

p <- ggplot(d, aes(x = x, y = p, group = factor(k))) + theme_minimal() + 
   geom_line() + scale_x_continuous(breaks = -6:6) + facet_grid(. ~ j) + 
   labs(x = "$\\zeta$", y = "Probability") + 
   theme(axis.text = element_text(color = "black", size = 6)) + ylim(0,1) + 
  theme(axis.title.y = element_text(size = 10))
plot(p)

if (prnt) graphics.off()
