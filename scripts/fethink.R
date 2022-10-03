library(gssr)
library(dplyr)
library(tidyr)
library(haven) # use as_factor to convert to factor
library(forcats)
library(VGAM)

setwd(this.path::this.dir())

if (!exists("gss_all")) data(gss_all)
if (!exists("gss_doc")) data(gss_doc)

set.seed(123)
d <- gss_all %>% 
  select(year, age, sex, race, fethink) %>%
  mutate(sex = as_factor(sex)) %>% 
  mutate(race = as_factor(race)) %>% na.omit()
    
# d <- d %>% filter(partyid != "other party") %>%
#   mutate(party = case_when(
#     grepl("independent", partyid, fixed = TRUE) ~ "I",
#     grepl("democrat", partyid, fixed = TRUE) ~ "D",
#     grepl("republican", partyid, fixed = TRUE) ~ "R"
#   ))

# d <- d %>% mutate(polviews = fct_collapse(polviews,
#   EL = c("extremely liberal","liberal"),
#   SL = "slightly liberal",
#   MD = "moderate, middle of the road",
#   SC = "slightly conservative",
#   EC = c("extremely conservative","conservative")
# )) %>% mutate(polviews = relevel(polviews, ref = "MD"))

d <- d %>% 
  filter(race != "other") %>%
  mutate(race = droplevels(relevel(race, ref = "white")))

d <- d %>% filter(age < 89) %>% mutate(age = as.numeric(age))

m <- rrvglm(fethink ~ I(age/10) + race + sex, multinomial, data = d)
AIC(m)
m <- vglm(fethink ~ I(age/10) + race + sex, multinomial(refLevel = 1), data = d)
AIC(m)

x_mnl <- model.matrix(fethink ~ I(age/10) + race + sex, data = d)
x_stm <- x_mnl[,-1]
y <- matrix(0, nrow(d), 4)
y[,1] <- d$fethink

theta_mnl <- as.vector(elmnl(rep(0, length(coef(m))), y, x_mnl, 1, FALSE)$estimate)
theta_stm <- as.vector(elstm(c(0,0,0,0.25,0.75,0,0,0), y, x_stm, 1, FALSE)$estimate)

set.seed(101)
y <- mnl_samp(theta_mnl, x_mnl, 4)

tmp1 <- elmnl(theta_mnl, y, x_mnl, 1, TRUE)
tmp2 <- elmnl(theta_mnl, y, x_mnl, 2, TRUE)
tmpK <- elmnl(theta_mnl, y, x_mnl, 3, TRUE)
tmp0 <- flmnl(theta_mnl, y, x_mnl, TRUE)

foo1 <- elstm(theta_stm, y, x_stm, 1, TRUE)
foo2 <- elstm(theta_stm, y, x_stm, 2, TRUE)
fooK <- elstm(theta_stm, y, x_stm, 3, TRUE)
foo0 <- flstm(theta_stm, y, x_stm, TRUE)

estimates <- cbind(tmp1$estimate, tmp2$estimate, tmpK$estimate, tmp0$estimate)
stderrors <- cbind(tmp1$stderror, tmp2$stderror, tmpK$stderror, tmp0$stderror)
aic <- c(foo1$aic - tmp1$aic, foo2$aic - tmp2$aic, fooK$aic - tmpK$aic, foo0$aic - tmp0$aic)

tbl <- data.frame(
  e1 = estimates[,1], s1 = stderrors[,1],
  e2 = estimates[,2], s2 = stderrors[,2], 
  e3 = estimates[,3], s3 = stderrors[,3],
  e4 = estimates[,4], s4 = stderrors[,4])

for (i in 1:ncol(tbl)) {
  tbl[,i] <- format(round(tbl[,i], 2), nsmall = 2)
}
for (i in c(2,4,6,8)) {
  tbl[,i] <- paste("(", tbl[,i], ")", sep = "")
}

parms <- c(paste("$\\beta_{0,", 2:4, "}$", sep = ""),
           paste("$\\beta_{1,", 2:4, "}$", sep = ""),
           paste("$\\beta_{2,", 2:4, "}$", sep = ""),
           paste("$\\beta_{3,", 2:4, "}$", sep = ""))

tbl <- cbind(parms, tbl)

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
textable(tbl, "../manuscript/fethink-table.tex")
