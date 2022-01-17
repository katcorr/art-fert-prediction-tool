# build algorithm on SART CORS data

library(tidyverse)
library(VGAM)           # for fitting right-censored Poisson
library(betareg)        # for beta regression 
library(rpart)          # for regression tree (CART)

# --------------------- load the data ----------------------------------------

sart <- readRDS("Z:/SARTCORS/sart.RDS")

#mosaic::favstats(~num_retrieved, data=sart)
#mosaic::favstats(~num_retrieved, data=filter(sart, num_retrieved>0))


# ----------- get woman-level dataset for step 1 modeling of algorithm --------

# step 1 model: predict number of ET transferred before one live birth 

# cycle order is never missing. good. can order by that
which(is.na(sart$cycle_order)==TRUE)

# cycle id should be unique but there are some duplicated
# for some reason the only thing different between rows with same external_cycle_id
# is the cycle_number
# keep only the row with the min cycle_number
length(unique(sart$external_cycle_id))

test1a <- sart %>%
  group_by(external_cycle_id) %>%
  mutate(num=row_number()) %>%
  filter(num>1)

test1b <- sart %>%
  filter(external_cycle_id %in% test1a$external_cycle_id)

test1c <- test1b %>%
  select(-cycle_order) %>%
  distinct()

length(unique(test1c$external_cycle_id))

test2 <- sart %>%
  select(-cycle_order) %>%
  distinct()

# remove duplicate cycle ids
remove_dups <- sart %>%
  group_by(external_cycle_id) %>%
  mutate(num=row_number()) %>%
  filter(num>1) 

dat_step1_0 <- sart %>%
  # removes the 41 duplicated cycle_id's
  anti_join(remove_dups, by=c("external_patient_id", "external_cycle_id", "cycle_order")) %>%
  select(external_patient_id, external_cycle_id, cycle_order, reporting_year, num_transferred, lb) %>% 
  arrange(external_patient_id, cycle_order, external_cycle_id) %>%
  group_by(external_patient_id) %>%
  mutate(cycnum = row_number()
         , totalET = cumsum(num_transferred)
         , totalLB = cumsum(lb)
         , totalETminus1 = totalET - 1
         , Keep = (lb == 1 & totalLB == 1) |
           (row_number() == max(row_number()) & totalLB==0))

# dat_step1_1 should end up keeping this many rows 
length(unique(sart$external_patient_id))

dat_step1_1 <- dat_step1_0 %>%
  filter(Keep == TRUE) %>%
  mutate(censored = ifelse(totalLB==1,yes=0,no=1)) %>%
  select(external_patient_id, totcyc = cycnum, totalET, totalETminus1, lb, censored) %>%
  ungroup()

# remove women who never had a transfer (for the step1 modeling to determine # ET needed for 1 LB)
dat_step1_2 <- dat_step1_1 %>% 
  filter(totalET > 0)

# mosaic::tally(~censored, data = dat_step1_2)
# mosaic::tally(~censored, data = dat_step1_2, format = "percent")

## next, get first cycle info to use as predictors in the model
first_cyc <- dat_step1_0 %>%
  filter(cycnum==1) %>%
  select(external_patient_id, external_cycle_id, cycle_order, reporting_year)

dat_step1_3 <- sart %>%
  filter(external_cycle_id %in% first_cyc$external_cycle_id) %>%
  distinct() %>%
  select(external_cycle_id, external_patient_id, reporting_year, cycle_order
         , clinic_state, clinic_state2, clinic_region_usa
         , patient_age_at_start, partner_age_at_start_c, partner_age_at_start_missing
         , bmi_c, bmicat4, bmicat6
         , smoker, gravidity_cat4, parity_cat4, max_fsh, fsh_gt10, amh_last_value, amhcat3
         , male_infertility, endometriosis, polycystic_ovaries, diminished_ovarian_reserve
         , tubal_ligation, tubal_hydrosalpinx, tubal_other, dx_tubal
         , uterine, unexplained
         , num_retrieved
         ) %>%
  filter(!is.na(patient_age_at_start) & !is.na(bmicat6) & !is.na(smoker)
          & !is.na(gravidity_cat4) & !is.na(parity_cat4)
          & !is.na(fsh_gt10) & !is.na(dx_tubal)
          & !is.na(uterine) & !is.na(unexplained) & !is.na(num_retrieved))

dat_step1 <- dat_step1_3 %>%
  inner_join(dat_step1_2, by=c("external_patient_id")) %>%
  ungroup() 

#mosaic::favstats(~totalET, data=dat_step1)
#mosaic::favstats(totalET ~ censored, data=dat_step1)
#mosaic::favstats(~num_retrieved, data=dat_step1)
#mosaic::favstats(num_retrieved ~ censored, data=dat_step1)
#count(dat_step1, clinic_state2, clinic_region_usa, clinic_state)

# dat_step1 %>%
#   select(-external_cycle_id, -external_patient_id) %>%
#   gtsummary::tbl_summary(by=censored)

# dat_step1 %>%
#   select(-external_cycle_id, -external_patient_id) %>%
#   gtsummary::tbl_summary(by=reporting_year)

# ------------------- divide into training and testing set ---------------------

n <- nrow(dat_step1)

set.seed(20211119)

dat_step1_r <- dat_step1 %>%
  mutate(randomnum = rnorm(n=n,mean=0,sd=1)) %>%
  arrange(randomnum)

train_step1 <- dat_step1_r %>%
  slice(1:(n/2))

test_step1 <- dat_step1_r %>%
  slice(((n/2)+1):n)

# ------------------------------------------------------------------------------
# ------------------- MODEL 1: Predicting #ET needed for 1 LB -------------------
# ------------------------------------------------------------------------------

# Predict the total number of embryos that will need to be transferred across cycles
# in order to achieve one live birth

pred_step1 <- "patient_age_at_start + I(patient_age_at_start^2)
               + clinic_state2 + bmicat6 + smoker + gravidity_cat4 + parity_cat4 + fsh_gt10 + amhcat3
               + male_infertility + endometriosis + polycystic_ovaries + diminished_ovarian_reserve
               + dx_tubal + uterine + unexplained + num_retrieved"

# right-censored Poisson model (re-characterized)
mod1 <- vglm(as.formula(str_replace_all(paste0("SurvS4(totalETminus1,lb) ~ ", pred_step1), "\n", ""))
                , family = cens.poisson()
                , data = train_step1)

mod1
mod1_coeff <- coefvlm(mod1)

# ----------------------- GET PREDICTIONS ON TRAIN SET ----------------------- 

pred_train <- train_step1 %>%
  mutate(mod1_pred0 = predict(mod1, newdata=train_step1, type="response")
         , mod1_pred = mod1_pred0 + 1
         , mod1_predC = ceiling(mod1_pred)
         ##### differences
         , mod1_diff = mod1_pred - totalET
         , mod1_diffC = mod1_predC - totalET
         ##### directions
         , mod1_dir = case_when(mod1_diff == 0 ~ "Exactly equal"
                                , mod1_diff < 0 ~ "Underestimate"
                                , mod1_diff > 0 ~ "Overestimate")
         , mod1_dirC = case_when(mod1_diffC == 0 ~ "Exactly equal"
                                 , mod1_diffC < 0 ~ "Underestimate"
                                 , mod1_diffC > 0 ~ "Overestimate")
        )

# ----------------------- GET PREDICTIONS ON TEST SET ----------------------- 

pred_test <- test_step1 %>%
  mutate(mod1_pred0 = predict(mod1, newdata=test_step1, type="response")
         , mod1_pred = mod1_pred0 + 1
         , mod1_predC = ceiling(mod1_pred)
         ##### differences
         , mod1_diff = mod1_pred - totalET
         , mod1_diffC = mod1_predC - totalET
         ##### directions
         , mod1_dir = case_when(mod1_diff == 0 ~ "Exactly equal"
                                , mod1_diff < 0 ~ "Underestimate"
                                , mod1_diff > 0 ~ "Overestimate")
         , mod1_dirC = case_when(mod1_diffC == 0 ~ "Exactly equal"
                                 , mod1_diffC < 0 ~ "Underestimate"
                                 , mod1_diffC > 0 ~ "Overestimate")
  )

# ---------------------- SAVE PREDICTIONS FROM STEP 1 -------------------------

pred_ET <- bind_rows(mutate(pred_train, set="Training")
                        , mutate(pred_test, set="Testing"))

# SAVE MODEL 1 PREDICTIONS 
saveRDS(pred_ET, "Z:/SARTCORS/pred_ET.RDS")

# ------------------------------------------------------------------------------
# ------------------- MODEL 2: Predicting fertilization rate  -------------------
# ------------------------------------------------------------------------------

# mosaic::favstats(num_retrieved ~ cycle_type, data=sart)
# mosaic::favstats(~total2pn, data=sart)
# mosaic::favstats(~total2pn, data=filter(sart, total2pn < 3000))
# mosaic::favstats(total2pn ~ cycle_type, data=sart)

#sart %>% count(cycle_type, retrieval_type_autologous_retrieval1)

ggplot(sart, aes(x=num_retrieved)) +
  geom_density()

ggplot(filter(sart, total2pn < 3000), aes(x=total2pn)) +
  geom_density()

train_step2_0 <- sart %>%
  # NEED TO SUBSET ON SAME TRAINING AND TESTING WOMEN AS FOR STEP 1 
  # so can assess final prediction (the ratio) in the end
  filter(external_cycle_id %in% train_step1$external_cycle_id) %>%
  filter(cycle_type == "Fresh") 

train_step2_1 <- train_step2_0 %>%
  # remove 96 cycles (<0.1%) with data-entry errors (total2pn > number retrieved)
  # and cycles missing total2pn (~1,900 / 1%)
  filter(total2pn <= num_retrieved) %>%
  # remove 225 cycles where num retrieved and num 2pn are both 0
  filter(num_retrieved != 0)

# mosaic::favstats(~total2pn, data=train_step2_1)

n_fert <- nrow(train_step2_1)

train_step2 <- train_step2_1 %>%
  # don't have #M2 or #inseminated in SART data so use num_retrieved as proxy
  # expect about 70% of eggs retrieved are mature and/or inseminated (on average)
  # so this proxy will provide conservative estimate
  mutate(partner_age = ifelse(is.na(partner_age_at_start)
                                , yes = -99
                                , no = partner_age_at_start)
         , fert_rate_proxy = total2pn/num_retrieved
         , t1 = (fert_rate_proxy*(n_fert-1) + 0.5)/n_fert
         , t2 = log(t1/(1-t1))) %>%
  select(-pregnancy_outcome_start_date_diff)
  
#length(unique(train_step2$external_patient_id))

# mosaic::favstats(~fert_rate_proxy, data=train_step2)
# mosaic::favstats(~t1, data=train_step2)
# mosaic::favstats(~t2, data=train_step2)

ggplot(train_step2) +
  geom_density(aes(x=fert_rate_proxy))

# train_step2 %>%
#   group_by(clinic_state2) %>%
#   summarize(mean_fert = mean(fert_rate_proxy)
#             , median_fert = median(fert_rate_proxy)) %>%
#   ggplot() +
#   geom_density(aes(x=mean_fert), color="red") +
#   geom_density(aes(x=median_fert), color="blue") 

ggplot(data=train_step2, aes(x=fert_rate_proxy, color=clinic_state2)) +
  geom_boxplot()

# ggplot(sart_fert, aes(x=num_retrieved, y=fert_rate_proxy)) +
#   geom_point()

pred_mod2 <- "patient_age_at_start + I(patient_age_at_start^2) + partner_age + I(partner_age^2) + partner_age_at_start_missing
               + clinic_state2 + bmicat6 + smoker + gravidity_cat4 + parity_cat4 + fsh_gt10 + amhcat3
               + male_infertility + endometriosis + polycystic_ovaries + diminished_ovarian_reserve
               + dx_tubal + uterine + unexplained + num_retrieved"

# https://www.jstatsoft.org/article/view/v034i02
# can't use directly on "fert_rate_proxy" because some obs have value of 1
# get error: Error in betareg(...): invalid dependent variable, all observations must be in (0, 1)
# mod2 <- betareg(fert_rate_proxy ~ patient_age_at_start + I(patient_age_at_start^2)
#                 + state + bmicat6 + smoker + gravidity_cat4 + fsh_gt10 + amhcat3
#                 + male_infertility + endometriosis + polycystic_ovaries + diminished_ovarian_reserve
#                 + dx_tubal + uterine + unexplained
#                 , data=sart_fert)


mod2_bML <- betareg(as.formula(str_replace_all(paste0("t1 ~ ", pred_mod2), "\n", ""))
                     , data=train_step2
                     , type = "ML")
names(mod2_bML)

# bias-corrected 
# mod2_bBC <- betareg(as.formula(str_replace_all(paste0("t1 ~ ", pred_mod2), "\n", ""))
#                     , data=train_step2
#                     , type="BC")
# summary(mod2_bBC)

# bias-reduced 
# mod2_bBR <- betareg(as.formula(str_replace_all(paste0("t1 ~ ", pred_mod2), "\n", ""))
#                   , data=train_step2
#                   , type="BR")
# summary(mod2_bBC)

# linear regression
mod2_lm <- lm(as.formula(str_replace_all(paste0("fert_rate_proxy ~ ", pred_mod2), "\n", ""))
              , data = train_step2)

summary(mod2_lm)



# transformation to whole real line (t2 not constrained to (0,1))
mod2_lmt <- lm(as.formula(str_replace_all(paste0("t2 ~ ", pred_mod2), "\n", ""))
               , data = train_step2)

summary(mod2_lmt)

# CART
# https://uc-r.github.io/regression_trees
# https://www.statmethods.net/advstats/cart.html
# https://cran.r-project.org/web/packages/rpart/vignettes/longintro.pdf

mod2_cart <- rpart(as.formula(str_replace_all(paste0("fert_rate_proxy ~ ", pred_mod2), "\n", ""))
                   , data = train_step2
                   , method = "anova"   # "class" for a classification tree and "anova" for a regression tree
                   , control=rpart.control(minsplit=2, minbucket=1, cp=0.001)
                   )

summary(mod2_cart)

# ----------------------- GET PREDICTIONS ON TRAIN SET ----------------------- 

pred_train_fert <- train_step2 %>%
  mutate(mod2_bML_pred = predict(mod2_bML)
         , mod2_lm_pred = predict(mod2_lm)
         , mod2_lmt_pred00 = predict(mod2_lmt)
         , mod2_lmt_pred0 = exp(mod2_lmt_pred00)/(1+exp(mod2_lmt_pred00))
         , mod2_lmt_pred = (mod2_lmt_pred0*n_fert - 0.5)/(n_fert-1)
         , mod2_cart_pred = predict(mod2_cart, type="vector")
         ##### differences
         , mod2_bML_diff = mod2_bML_pred - fert_rate_proxy
         , mod2_lm_diff = mod2_lm_pred - fert_rate_proxy
         , mod2_lmt_diff = mod2_lmt_pred - fert_rate_proxy
         , mod2_cart_diff = mod2_cart_pred - fert_rate_proxy
         ##### directions
         , across(ends_with("diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                 , . < 0 ~ "Underestimate"
                                                 , . > 0 ~ "Overestimate")
                  , .names = "{.col}_dir"))

# compare predicted dist'ns
# pred_train_fert_long <- pred_train_fert %>%
#   pivot_longer(cols=c("fert_rate_proxy",ends_with("_pred")), names_to="method", values_to="value") %>%
#   select(external_cycle_id, t1, t2, method, value)
# 
# ggplot(pred_train_fert_long, aes(x=value, color=method)) +
#   geom_density()

# ----------------------- GET PREDICTIONS ON TEST SET ----------------------- 

test_step2_0 <- sart %>%
  # NEED TO SUBSET ON SAME TRAINING AND TESTING WOMEN AS FOR STEP 1
  # so can assess final prediction (the ratio) in the end
  filter(external_cycle_id %in% test_step1$external_cycle_id) %>%
  filter(cycle_type == "Fresh") 

test_step2_1 <- test_step2_0 %>%
  filter(total2pn <= num_retrieved) %>%
  filter(num_retrieved != 0)

n_fert_test <- nrow(test_step2_1)

test_step2 <- test_step2_1 %>%
  mutate(partner_age = ifelse(is.na(partner_age_at_start)
                                , yes = -99
                                , no = partner_age_at_start)
         , fert_rate_proxy = total2pn/num_retrieved
         , t1 = (fert_rate_proxy*(n_fert_test-1) + 0.5)/n_fert_test
         , t2 = log(t1/(1-t1))) %>%
  select(-pregnancy_outcome_start_date_diff)

pred_test_fert <- test_step2 %>%
  mutate(mod2_bML_pred = predict(mod2_bML, newdata=test_step2)
         , mod2_lm_pred = predict(mod2_lm, newdata=test_step2)
         , mod2_lmt_pred00 = predict(mod2_lmt, newdata=test_step2)
         , mod2_lmt_pred0 = exp(mod2_lmt_pred00)/(1+exp(mod2_lmt_pred00))
         , mod2_lmt_pred = (mod2_lmt_pred0*n_fert_test - 0.5)/(n_fert_test-1)
         , mod2_cart_pred = predict(mod2_cart, newdata=test_step2, type="vector")
         ##### differences
         , mod2_bML_diff = mod2_bML_pred - fert_rate_proxy
         , mod2_lm_diff = mod2_lm_pred - fert_rate_proxy
         , mod2_lmt_diff = mod2_lmt_pred - fert_rate_proxy
         , mod2_cart_diff = mod2_cart_pred - fert_rate_proxy
         ##### directions
         , across(ends_with("diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                 , . < 0 ~ "Underestimate"
                                                 , . > 0 ~ "Overestimate")
                  , .names = "{.col}_dir"))

# compare predicted dist'ns
# pred_test_fert_long <- pred_test_fert %>%
#   pivot_longer(cols=c("fert_rate_proxy",ends_with("_pred")), names_to="method", values_to="value") %>%
#   select(external_cycle_id, t1, t2, method, value)
# 
# ggplot(pred_test_fert_long, aes(x=value, color=method)) +
#   geom_density()

# ---------------------- SAVE PREDICTIONS FROM MODEL 2 -------------------------

pred_fert <- bind_rows(mutate(pred_train_fert, set="Training")
                        , mutate(pred_test_fert, set="Testing"))

# SAVE PREDICTIONS ON TRAINING SET
#saveRDS(pred_fert, "Z:/SARTCORS/pred_fert.RDS")


# ------------------------------------------------------------------------------
# -------------------       GET FINAL PREDICTIONS         ----------------------
# ------------------------------------------------------------------------------

# ----------------------- GET FINAL PREDICTIONS ON TRAIN SET ------------------- 

names(pred_train)
names(pred_train_fert)
length(unique(pred_train_fert$external_patient_id))

pred_train_final <- pred_train %>%
  select(external_cycle_id, external_patient_id, totalET, anylb=lb, censored, starts_with("mod1_")) %>%
  inner_join(pred_train_fert, by=c("external_cycle_id", "external_patient_id")) %>%
  mutate(pred_pois_cart = ceiling(mod1_pred/mod2_cart_pred)
         , pred_pois_lm = ceiling(mod1_pred/mod2_lm_pred)
         , pred_pois_lmt = ceiling(mod1_pred/mod2_lmt_pred)
         , pred_pois_bML = ceiling(mod1_pred/mod2_bML_pred)
         # differences from total # of embryos transferred across cycles
         , across(starts_with("pred_"), ~ . - totalET, .names = "{.col}_ETdiff")
         , across(ends_with("_ETdiff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                 , . < 0 ~ "Underestimate"
                                                 , . > 0 ~ "Overestimate")
                  , .names = "{.col}_dir")
         # differences from number of oocytes available
         , across(starts_with("pred_") & !ends_with("diff") & !ends_with("_dir")
                  , ~ num_retrieved - ., .names = "{.col}_Oocdiff")
         , across(ends_with("_Oocdiff"), ~ case_when(.== 0 ~ "All oocytes"
                                                    , . > 0 ~ "Surplus oocytes"
                                                    , . < 0 ~ "Not enough oocytes")
                  , .names = "{.col}_dir")
         # differences from 70% of number of oocytes available (since generally only 70% mature / available for insemination)
         , across(starts_with("pred_") & !ends_with("diff") & !ends_with("_dir")
                  , ~ num_retrieved*0.70 - ., .names = "{.col}_Ooc.7diff")
         , across(ends_with("_Ooc.7diff"), ~ case_when(.== 0 ~ "70% of oocytes"
                                                     , . > 0 ~ "Surplus beyond 70% of oocytes"
                                                     , . < 0 ~ "Not enough if need 70% of oocytes")
                  , .names = "{.col}_dir")
         ) %>%
  select(external_patient_id, external_cycle_id, patient_age_at_start
         , censored, anylb, num_retrieved, totalET, fert_rate_proxy
         , starts_with("pred_")
         , everything())


# ----------------------- GET FINAL PREDICTIONS ON TEST SET ------------------- 

pred_test_final <- pred_test %>%
  select(external_cycle_id, external_patient_id, totalET, anylb=lb, censored, starts_with("mod1_")) %>%
  inner_join(pred_test_fert, by=c("external_cycle_id", "external_patient_id")) %>%
  mutate(pred_pois_cart = ceiling(mod1_pred/mod2_cart_pred)
         , pred_pois_lm = ceiling(mod1_pred/mod2_lm_pred)
         , pred_pois_lmt = ceiling(mod1_pred/mod2_lmt_pred)
         , pred_pois_bML = ceiling(mod1_pred/mod2_bML_pred)
         # differences from total # of embryos transferred across cycles
         , across(starts_with("pred_"), ~ . - totalET, .names = "{.col}_ETdiff")
         , across(ends_with("_ETdiff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                    , . < 0 ~ "Underestimate"
                                                    , . > 0 ~ "Overestimate")
                  , .names = "{.col}_dir")
         # differences from number of oocytes available
         , across(starts_with("pred_") & !ends_with("diff") & !ends_with("_dir")
                  , ~ num_retrieved - ., .names = "{.col}_Oocdiff")
         , across(ends_with("_Oocdiff"), ~ case_when(.== 0 ~ "All oocytes"
                                                     , . > 0 ~ "Surplus oocytes"
                                                     , . < 0 ~ "Not enough oocytes")
                  , .names = "{.col}_dir")
         # differences from 70% of number of oocytes available (since generally only 70% mature / available for insemination)
         , across(starts_with("pred_") & !ends_with("diff") & !ends_with("_dir")
                  , ~ num_retrieved*0.70 - ., .names = "{.col}_Ooc.7diff")
         , across(ends_with("_Ooc.7diff"), ~ case_when(.== 0 ~ "70% of oocytes"
                                                       , . > 0 ~ "Surplus beyond 70% of oocytes"
                                                       , . < 0 ~ "Not enough if need 70% of oocytes")
                  , .names = "{.col}_dir")
  ) %>%
  select(external_patient_id, external_cycle_id, patient_age_at_start
         , censored, anylb, num_retrieved, totalET, fert_rate_proxy
         , starts_with("pred_")
         , everything())


# ---------------------- SAVE FINAL PREDICTIONS  -------------------------

pred_final <- bind_rows(mutate(pred_train_final, set="Training")
                       , mutate(pred_test_final, set="Testing"))

# SAVE FINAL PREDICTIONS
#saveRDS(pred_final, "Z:/SARTCORS/pred_final.RDS")
