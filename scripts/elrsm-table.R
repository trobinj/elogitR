library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/home/trjohns/Dropbox/Research/elogit/manuscript")

d <- read.csv("elrsm-simdata.csv")
K <- 5

reps <- nrow(d)/5

pnames <- c(paste("$\\tau_", 3:5, "$", sep = ""),
   paste("$\\delta_", 1:3, "$", sep = ""), "$\\sigma$")

d <- d %>% mutate(sample = rep(1:reps, 5), K = rep(c(1,2,3,K,0), each = reps)) %>% 
   pivot_longer(cols = c("t3","t4","t5","sigma","d1","d2","d3"), 
      names_to = "parameter", values_to = "estimate")

p <- ggplot(d, aes(x = estimate)) + theme_minimal() + facet_grid(K ~ parameter) + 
   geom_histogram()
plot(p)

d <- d %>% group_by(parameter, K) %>%
   summarize(rmse = sqrt(mean(estimate^2)),
      bias = mean(estimate), sd = sd(estimate)) %>%
   pivot_longer(cols = c(rmse, bias, sd), names_to = "type", values_to = "y") %>%
   mutate(y = format(round(y, 3), nsmall = 3)) %>%
   arrange(type, parameter, K) %>%
   pivot_wider(names_from = K, values_from = y) %>%
   arrange(factor(type, levels = c("rmse","bias","sd")),
      factor(parameter, levels = c("t3","t4","t5","d1","d2","d3","sigma"))) %>%
   ungroup() %>% select(-c(1,2)) %>%
   mutate(parameter = rep(pnames, 3)) %>%
   mutate(type = c(c("RMSE", rep("", 6)), c("Bias", rep("", 6)), c("SE", rep("", 6)))) %>%
   select(c(7,6,2,3,4,5,1))

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

textable(d, "elrsm-table.tex", midrule = c(7,14))
