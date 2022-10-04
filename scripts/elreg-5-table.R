library(dplyr)

setwd("/home/trjohns/Dropbox/Research/elogit/manuscript")

d <- read.csv("elreg-5-simdata.csv")

K <- 5
reps <- nrow(d)/K

pnames <- c(paste("$\\alpha_", 2:5, "$", sep = ""),
   paste("$\\gamma_", 2:4, "$", sep = ""), "$\\beta$")

d <- d %>% mutate(g2 = g2 * (K-1), g3 = g3 * (K-1), g4 = g4 * (K-1), beta = beta / (K - 1))

d <- d %>% mutate(sample = rep(1:reps, 5), K = rep(c(1,2,3,K,0), each = reps)) %>% 
   pivot_longer(cols = c("a2","a3","a4","a5","g2","g3","g4","beta"), 
      names_to = "parameter", values_to = "estimate")

d <- d %>% group_by(parameter, K) %>%
   summarize(rmse = sqrt(mean(estimate^2)),
      bias = mean(estimate), sd = sd(estimate)) %>%
   pivot_longer(cols = c(rmse, bias, sd), names_to = "type", values_to = "y") %>%
   mutate(y = format(round(y, 3), nsmall = 3)) %>% 
   arrange(type, parameter, K) %>% 
   pivot_wider(names_from = K, values_from = y) %>%
   arrange(factor(type, levels = c("rmse","bias","sd")),
      factor(parameter, levels = c("a2","a3","a4","a5","g2","g3","g4","beta"))) %>%
   select(c(2,1,4,5,6,7,3)) %>% ungroup() %>% select(-c(1,2)) %>%
   mutate(parameter = rep(pnames, 3)) %>%
   mutate(type = c(c("RMSE", rep("", 7)), c("Bias", rep("", 7)), c("SE", rep("", 7)))) %>%
   select(c(7,6,1,2,3,4,5))

textable <- function(x, file, midrule) {
   for (i in 1:nrow(d)) {
      if (i %in% midrule) {
         eolstring <- " \\\\ \\midrule \n"
      } else if (i == nrow(d)) {
         eolstring <- " \\tabularnewline\\bottomrule"
      } else {
         eolstring <- " \\\\\n"
      }
      write.table(d[i,], file, append = i > 1, row.names = FALSE, col.names = FALSE,
                  quote = FALSE, sep = " & ", na = "", eol = eolstring)
   }
}

textable(d, "elreg-5-table.tex", midrule = c(8,16))

