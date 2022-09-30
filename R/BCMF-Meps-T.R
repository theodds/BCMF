## Load (37) ----

library(SoftBart)
library(progress)
library(tidyverse)
library(caret)
source("R/softbart_bcmf.R")

meps <- read_csv("Data/meps2011.csv") %>%
  filter(totexp > 0 & bmi > 0) %>%
  mutate(logY = log(totexp)) %>%
  mutate(smoke = ifelse(smoke == "No", 0, 1)) %>%
  select(-totexp) %>% select(logY, phealth, smoke, everything())

phealth <- meps$phealth
phealth <- case_when(phealth == "Poor" ~ 1, phealth == "Fair" ~ 2,
                     phealth == "Good" ~ 3, phealth == "Very Good" ~ 4,
                     phealth == "Excellent" ~ 5)
meps$phealth <- phealth
rm(phealth)

## Fit ----

formula_y <- logY ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_m <- phealth ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

opts <- Opts(num_burn = 2000, num_thin = 1, num_save = 2000, update_s = FALSE)

set.seed(digest::digest2int("BCMF-Meps-T.R"))

my_fit <- softbart_bcmf(formula_y = formula_y, formula_m = formula_m, 
                        trt = meps$smoke, data = meps, test_data = meps, 
                        opts = opts)

hist(rowMeans(my_fit$delta_train))

saveRDS(my_fit, file = "Data/meps_fit.rds")