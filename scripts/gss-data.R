library(gssr)
library(dplyr)
library(tidyr)
library(haven) # use as_factor to convert to factor
library(forcats)

setwd(this.path::this.dir())

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
