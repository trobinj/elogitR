library(dplyr)
library(tidyr)

setwd("/home/trjohns/Dropbox/Research/elogit/manuscript")

d <- read.csv("elreg-3-simdata.csv")

K <- 3
reps <- nrow(d)/2

pnames <- c(paste("$\\alpha_", 2:3, "$", sep = ""),
   paste("$\\gamma_", 2:2, "$", sep = ""), "$\\beta$")

d <- d %>% mutate(g2 = g2 * (K-1), beta = beta / (K - 1))

d <- d %>% mutate(sample = rep(1:reps, 2), K = rep(c(1,K), each = reps)) %>% 
   pivot_longer(cols = c("a2","a3","g2","beta"), 
      names_to = "parameter", values_to = "estimate")

d <- d %>% group_by(parameter, K) %>%
   summarize(rmse = sqrt(mean(estimate^2)),
      bias = mean(estimate), sd = sd(estimate)) %>%
   pivot_longer(cols = c(rmse, bias, sd), names_to = "type", values_to = "y") %>%
   mutate(y = format(round(y, 3), nsmall = 3)) %>% 
   arrange(type, parameter, K) %>% 
   pivot_wider(names_from = K, values_from = y) %>%
   arrange(factor(type, levels = c("rmse","bias","sd")),
      factor(parameter, levels = c("a2","a3","g2","beta"))) %>%
   ungroup() %>% select(-c(1,2)) %>%
   mutate(parameter = rep(pnames, 3)) %>%
   mutate(type = c(c("RMSE", rep("", 3)), c("Bias", rep("", 3)), c("SE", rep("", 3)))) %>%
   select(c(4,3,1,2))

textable <- function(x, file, midrule) {
   for (i in 1:nrow(x)) {
      if (i %in% midrule) {
         eolstring <- " \\\\ \\midrule \n"
      } else if (i == nrow(x)) {
         eolstring <- " \\tabularnewline\\bottomrule"
      } else {
         eolstring <- " \\\\\n"
      }
      write.table(x[i,], file, append = i > 1, row.names = FALSE, col.names = FALSE,
         quote = FALSE, sep = " & ", na = "", eol = eolstring)
   }
}

textable(d, "elreg-3-table.tex", midrule = c(4,8))

