library(VGAM)
library(dplyr)
library(tidyr)
library(gssr)
library(haven) # use as_factor to convert to factor
library(forcats)

setwd("/home/trjohns/Dropbox/Research/elogit/manuscript")

if (!exists("gss_all")) data(gss_all)
if (!exists("gss_doc")) data(gss_doc)

set.seed(123)
d <- gss_all %>% 
  select(year, age, sex, race, partyid, polviews) %>%
  mutate(sex = as_factor(sex)) %>% 
  mutate(race = as_factor(race)) %>%
  mutate(partyid = as_factor(partyid)) %>%
  mutate(polviews = as_factor(polviews)) %>% drop_na()

d <- d %>% filter(partyid != "other party") %>%
  mutate(party = case_when(
    grepl("independent", partyid, fixed = TRUE) ~ "I",
    grepl("democrat", partyid, fixed = TRUE) ~ "D",
    grepl("republican", partyid, fixed = TRUE) ~ "R"
  ))

d <- d %>% mutate(polviews = fct_collapse(polviews,
  EL = c("extremely liberal","liberal"),
  SL = "slightly liberal",
  MD = "moderate, middle of the road",
  SC = "slightly conservative",
  EC = c("extremely conservative","conservative")
))

d <- d %>% 
  filter(race != "other") %>%
  mutate(race = droplevels(relevel(race, ref = "white"))) %>% 
  mutate(party = factor(party)) %>%
  mutate(party = relevel(party, ref = "I"))

d <- d %>% filter(year == 2014) %>%
  mutate(y = as.numeric(polviews)) 

rm("gss_all", "gss_doc")
gc()

m <- rrvglm(y ~ party * race - race + age, multinomial, data = d)
AIC(m) # 5709
m <- vglm(y ~ party * race + age, multinomial(refLevel = 1), data = d)
summary(m)
AIC(m) # 5715

x <- model.matrix(y ~ party * race - race + I(age/10), data = d)[,-1]
y <- matrix(0, nrow(d), 5)
y[,1] <- d$y

theta <- as.vector(elstm(c(0,0,0,0,0.5,0.5,0.5,0,0,0,0,0,0), y, x, 1, FALSE)$estimate)
alph <- c(0, theta[1:4])
gamm <- c(0, theta[5:7], 1)
beta <- c(theta[8:13])

set.seed(123)
y <- stm_samp(c(alph, gamm, beta), x, 5)

tmp1 <- elstm(theta, y, x, 1, TRUE)
tmp2 <- elstm(theta, y, x, 2, TRUE)
tmp3 <- elstm(theta, y, x, 3, TRUE)
tmpK <- elstm(theta, y, x, 4, TRUE)
tmp0 <- flstm(theta, y, x, TRUE)

estimates <- cbind(tmp1$estimate, tmp2$estimate, tmp3$estimate, tmpK$estimate, tmp0$estimate)
stderrors <- cbind(tmp1$stderror, tmp2$stderror, tmp3$stderror, tmpK$stderror, tmp0$stderror)

tmp <- stderrors
for (i in 1:ncol(tmp)) {
  tmp[,i] <- (tmp[,i] - stderrors[,1])/stderrors[,1] * 100
}
tmp

tbl <- data.frame(
  e1 = estimates[,1], s1 = stderrors[,1],
  e2 = estimates[,2], s2 = stderrors[,2], 
  e3 = estimates[,3], s3 = stderrors[,3],
  e4 = estimates[,4], s4 = stderrors[,4],
  e5 = estimates[,5], s5 = stderrors[,5])

for (i in 1:ncol(tbl)) {
  tbl[,i] <- format(round(tbl[,i], 2), nsmall = 2)
}
for (i in c(2,4,6,8,10)) {
  tbl[,i] <- paste("(", tbl[,i], ")", sep = "")
}

tbl <- tbl[c(4:7, 1:3, 8:13),]

param <- c("$\\alpha_2$","$\\alpha_3$","$\\alpha_4$","$\\alpha_5$",
  "$\\gamma_2$","$\\gamma_3$","$\\gamma_4$",
  paste("$\\beta_", 1:6, "$", sep = ""))

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
textable(tbl, "polviews-table.tex")