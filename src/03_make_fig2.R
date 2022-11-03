library(dplyr)
library(rpart)
library(rpart.plot)

meps <- readRDS('data/meps.rds')

# Individual direct and indirect effects
out_bart <- readRDS('data/out_bart.rds')
indv_zeta <- out_bart$zeta_samples
indv_delta <- out_bart$tau_samples * out_bart$d_samples

meps_post <- meps %>%
  mutate(zeta_hat = colMeans(indv_zeta),
         delta_hat = colMeans(indv_delta)) %>%
  select(age, bmi, edu, income, povlev, region, sex, marital, race, seatbelt, delta_hat, zeta_hat)

rpart.plot(rpart(delta_hat ~ . - zeta_hat, data = meps_post, 
                 control = rpart.control(cp = 0.03)),
           family = 'Times New Roman', nn.family = 'Times New Roman')
