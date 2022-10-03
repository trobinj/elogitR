library(dplyr)

setwd(this.path::this.dir())

d <- read.csv("elreg-4-simdata.csv")

K <- 4
reps <- nrow(d)/K

pnames <- c(paste("$\\alpha_", 2:4, "$", sep = ""),
   paste("$\\gamma_", 2:3, "$", sep = ""), "$\\beta$")

d <- d %>% mutate(g2 = g2 * (K-1), g3 = g3 * (K-1), beta = beta / (K - 1))

d <- d %>% mutate(sample = rep(1:reps, 4), K = rep(c(1,2,K,0), each = reps)) %>% 
   pivot_longer(cols = c("a2","a3","a4","g2","g3","beta"), 
      names_to = "parameter", values_to = "estimate")

d <- d %>% group_by(parameter, K) %>%
   summarize(rmse = sqrt(mean(estimate^2)),
      bias = mean(estimate), sd = sd(estimate)) %>%
   pivot_longer(cols = c(rmse, bias, sd), names_to = "type", values_to = "y") %>%
   mutate(y = format(round(y, 3), nsmall = 3)) %>% 
   arrange(type, parameter, K) %>% 
   pivot_wider(names_from = K, values_from = y) %>%
   arrange(factor(type, levels = c("rmse","bias","sd")),
      factor(parameter, levels = c("a2","a3","a4","g2","g3","beta"))) %>%
   select(c(2,1,4,5,6,3)) %>% ungroup() %>% select(-c(1,2)) %>%
   mutate(parameter = rep(pnames, 3)) %>%
   mutate(type = c(c("RMSE", rep("", 5)), c("Bias", rep("", 5)), c("SE", rep("", 5)))) %>%
   select(c(6,5,1,2,3,4))

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

textable(d, "../manuscript/elreg-4-table.tex", midrule = c(6,12))
