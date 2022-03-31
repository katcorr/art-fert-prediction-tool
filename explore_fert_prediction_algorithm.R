# fit different modeling approaches to SART CORS data for each step in algorithm 
# fit separate models for each age group to allow different relationships
# run 'assess_oc_explore.R' to compare
# kat correia

library(tidyverse)
library(VGAM)           # for fitting right-censored Poisson
library(survival)       # for survival function
library(betareg)        # for beta regression 
library(rpart)          # for regression tree (CART)
library(quantreg)       # for quantile regression (median regression)

# --------------------------- load data ---------------------------------------

sart_stim <- readRDS("Z:/R03/Data/sart_stim.RDS")

sart_transfers <- readRDS("Z:/R03/Data/sart_transfers.RDS") %>%
  mutate(external_cycle_id = ifelse(is.na(external_cycle_id), yes=frozen_cycle_id
                                    , no=external_cycle_id))

# --------------------- get training/testing assignments -----------------------

set_assign <- readRDS("Z:/R03/Data/set_assign.RDS")

# ------------------------------------------------------------------------------
# ----------- DATASET FOR MODEL 1 (ET) BY TRAIN & TEST SETS --------------------
# ------------------------------------------------------------------------------

# compute total number of embryos transferred acrosss cycles until first live birth 
dat_step1_0 <- sart_transfers %>%
  select(external_patient_id, external_cycle_id, cycle_order, reporting_year, num_transferred, lb) %>% 
  arrange(external_patient_id, cycle_order, external_cycle_id) %>%
  group_by(external_patient_id) %>%
  mutate(cycnum = row_number()
         , totalET = cumsum(num_transferred)
         , totalLB = cumsum(lb)
         , totalETminus1 = totalET - 1
         , Keep = (lb == 1 & totalLB == 1) |
           (row_number() == max(row_number()) & totalLB==0)
         , max_cycle = max(cycle_order)) %>%
  ungroup()

dat_step1_1 <- dat_step1_0 %>%
  ungroup() %>%
  filter(Keep == TRUE) %>%
  mutate(censored = ifelse(totalLB==1,yes=0,no=1)) %>%
  select(external_patient_id, cycle_order, totcyc = cycnum, totalET, totalETminus1, lb, censored)

quantile(dat_step1_1$totalET, p=c(0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99))
mosaic::tally(~censored, data = dat_step1_1)
mosaic::tally(~censored, data = dat_step1_1, format = "percent")

## then take first stimulation cycle info to use as predictors in the model
## need to account for fact that some transfers are frozen so need
## external cycle id needs to be taken from thaw field
first_cyc <- sart_stim %>%
  filter(cycle_order==1) %>%
  select(first_stim_cycle_id = external_cycle_id, external_patient_id
         , stim_reporting_year = reporting_year, stim_cycle_order = cycle_order
         , clinic_state, agegrp, numER
         , patient_age_at_start, partner_age_at_start_c, partner_age_at_start_missing
         , bmi_c, bmicat4, bmicat6
         , smoker, gravidity_cat4, parity_cat4, max_fsh, fsh_gt10, amh_last_value, amhcat3
         , male_infertility, endometriosis, dx_ovulation, diminished_ovarian_reserve
         , tubal_ligation, tubal_hydrosalpinx, tubal_other, dx_tubal
         , uterine, unexplained
         , num_retrieved
  )


# these should all be 1. Good!
count(first_cyc, stim_cycle_order)

# --------- merge w/ first stim cycle parameters & train/test assignment -------

dat_step1 <- dat_step1_1 %>%
  left_join(first_cyc, by="external_patient_id") %>%
  left_join(set_assign, by="external_patient_id")

quantile(dat_step1$totalET, p=c(0.10, 0.25, 0.5, 0.75, 0.90, 0.95, 0.99))

dat_step1_train <- dat_step1 %>%
  filter(set=="Training Set")

dat_step1_test <- dat_step1 %>%
  filter(set=="Testing Set")

# ----------- SPLIT INTO DIFF DATASETS FOR DIFF AGE GROUPS --------------------

# TRAINING DATA
dat_step1_train_a1 <- dat_step1_train %>%
  filter(agegrp=="<35")

dat_step1_train_a2 <- dat_step1_train %>%
  filter(agegrp=="35-37")

dat_step1_train_a3 <- dat_step1_train %>%
  filter(agegrp=="38-40")

dat_step1_train_a4 <- dat_step1_train %>%
  filter(agegrp=="41-42")

dat_step1_train_a5 <- dat_step1_train %>%
  filter(agegrp==">42")

# TESTING DATA
dat_step1_test_a1 <- dat_step1_test %>%
  filter(agegrp=="<35")

dat_step1_test_a2 <- dat_step1_test %>%
  filter(agegrp=="35-37")

dat_step1_test_a3 <- dat_step1_test %>%
  filter(agegrp=="38-40")

dat_step1_test_a4 <- dat_step1_test %>%
  filter(agegrp=="41-42")

dat_step1_test_a5 <- dat_step1_test %>%
  filter(agegrp==">42")

# ------------------------------------------------------------------------------
# ------------------- STEP 1: Predicting #ET needed for 1 LB -------------------
# ------------------------------------------------------------------------------

predictors1 <- "clinic_state + bmicat6 + smoker + gravidity_cat4 + parity_cat4 + fsh_gt10 + amhcat3
               + male_infertility + endometriosis + dx_ovulation + diminished_ovarian_reserve
               + dx_tubal + uterine + unexplained + numER"

predictors2 <- paste("agegrp +", predictors1)

# create function to fit the models to particular age grp
# and get back dataset w predicted values 

fit_mods <- function(dat){
  
  print(mosaic::tally(~censored, data=dat, format="percent"))
  
  if(deparse(substitute(dat))=="dat_step1_train_a45"){
    predictors_a <- predictors2
    
  }
  else{
    predictors_a <- predictors1
  }
  
  print(predictors_a)
  
  # right-censored Poisson model (re-characterized)
  mod1_cp <- vglm(as.formula(str_replace_all(paste0("SurvS4(totalETminus1,lb) ~ ", predictors_a), "\n", ""))
                  , family = cens.poisson()
                  , data = dat)
  
  #print(summary(mod1_cp))
  
  # lm model
  mod1_lm <- lm(as.formula(str_replace_all(paste0("totalET ~ ", predictors_a), "\n", ""))
                , data = dat)
  
  #print(summary(mod1_lm))
  
  # quantile regression model
  mod1_qr <- rq(as.formula(str_replace_all(paste0("totalET ~ ", predictors_a), "\n", ""))
                , data = dat
                #, tau = c(0.5, 0.75, 0.90))
                , tau = 0.5)
  
  #print(summary(mod1_qr))
  
  # ----------------------- GET PREDICTIONS ON TRAIN SET ----------------------- 

    lb_pred_mod <- glm(as.formula(str_replace_all(paste0("lb ~ ", predictors_a, " + totalET"), "\n", ""))
                     , data=dat
                     , family="binomial")
  
  #print(summary(lb_pred_mod))
  
  pred_train_dat <- dat %>%
    mutate(set = "Training" 
           , mod1_cp_pred0 = predict(mod1_cp, newdata=dat, type="response")
           , mod1_cp_pred = mod1_cp_pred0 + 1
           , mod1_cp_predR = round(mod1_cp_pred, digits=0)
           , mod1_cp_predC = ceiling(mod1_cp_pred)
           , mod1_lm_pred = predict(mod1_lm, newdata=dat)
           , mod1_lm_predR = round(mod1_lm_pred, digits=0)
           , mod1_lm_predC = ceiling(mod1_lm_pred)
           , mod1_q50_pred = predict(mod1_qr, newdata=dat)#[,1]
           , mod1_q50_predR = round(mod1_q50_pred, digits=0)
           , mod1_q50_predC = ceiling(mod1_q50_pred)
           # PROBABILITY OF LB
           , lb_pred = predict(lb_pred_mod, newdata=dat, type="response")
           , weight = (1-censored)/lb_pred
           ##### differences
           , across(starts_with("mod1_"), ~ . - totalET, .names = "{.col}_diff")
           ##### directions
           , across(ends_with("diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                   , . < 0 ~ "Underestimate"
                                                   , . > 0 ~ "Overestimate")
                    , .names = "{.col}_dir")
    ) %>%
    select(external_patient_id, totalET, lb, censored, weight, starts_with("mod1_"), everything()) %>%
    select(-mod1_cp_pred0_diff)
  
  # ----------------------- GET PREDICTIONS ON TEST SET ----------------------- 
  
  assign("test_dat", eval(as.name(str_replace(deparse(substitute(dat)), "train", "test"))))
  
  pred_test_dat <- test_dat %>%
    mutate(set = "Testing"
           , mod1_cp_pred0 = predict(mod1_cp, newdata=test_dat, type="response")
           , mod1_cp_pred = mod1_cp_pred0 + 1
           , mod1_cp_predR = round(mod1_cp_pred, digits=0)
           , mod1_cp_predC = ceiling(mod1_cp_pred)
           , mod1_lm_pred = predict(mod1_lm, newdata=test_dat)
           , mod1_lm_predR = round(mod1_lm_pred, digits=0)
           , mod1_lm_predC = ceiling(mod1_lm_pred)
           , mod1_q50_pred = predict(mod1_qr, newdata=test_dat)#[,1]
           , mod1_q50_predR = round(mod1_q50_pred, digits=0)
           , mod1_q50_predC = ceiling(mod1_q50_pred)
           # PROBABILITY OF LB
           , lb_pred = predict(lb_pred_mod, newdata=test_dat, type="response")
           , weight = (1-censored)/lb_pred
           ##### differences
           , across(starts_with("mod1_"), ~ . - totalET, .names = "{.col}_diff")
           ##### directions
           , across(ends_with("diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                   , . < 0 ~ "Underestimate"
                                                   , . > 0 ~ "Overestimate")
                    , .names = "{.col}_dir")
    ) %>%
    select(external_patient_id, totalET, lb, censored, weight, starts_with("mod1_"), everything()) %>%
    select(-mod1_cp_pred0_diff)
  
  # ----------------------- COMBINE PREDICTIONS ----------------------- 
  
  pred_dat <- bind_rows(pred_train_dat, pred_test_dat)
  
  return(pred_dat)
}

# ANALYSES STRATIFIED BY AGE
pred_step1_a1 <- fit_mods(dat=dat_step1_train_a1)    # 29% censored
pred_step1_a2 <- fit_mods(dat=dat_step1_train_a2)    # 40.6% censored
pred_step1_a3 <- fit_mods(dat=dat_step1_train_a3)    # 57% censored
#pred_step1_a4 <- fit_mods(dat=dat_step1_train_a4)   # 76% censored
#pred_step1_a5 <- fit_mods(dat=dat_step1_train_a5)   # 89% censored

# collapse last two groups and include the two age grps as covariate 
dat_step1_train_a45 <- bind_rows(dat_step1_train_a4, dat_step1_train_a5)
dat_step1_test_a45 <- bind_rows(dat_step1_test_a4, dat_step1_test_a5)

pred_step1_a45 <- fit_mods(dat=dat_step1_train_a45)    # 81% censored


# ---------------------- SAVE PREDICTIONS FROM STEP 1 -------------------------

quantile(dat_step1_train$totalET, p=c(0.50,0.75,0.90))

pred_ET_byage <- bind_rows(pred_step1_a1, pred_step1_a2, pred_step1_a3, pred_step1_a45) %>%
  mutate(perc50 = 2
         , perc75 = 3
         , perc90 = 5
         , across(c("perc50", "perc75", "perc90"), ~ . - totalET, .names = "{.col}_diff")
         , across(c("perc50_diff", "perc75_diff", "perc90_diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                                              , . < 0 ~ "Underestimate"
                                                                              , . > 0 ~ "Overestimate")
                  , .names = "{.col}_dir"))

saveRDS(pred_ET_byage, "Z:/R03/Data/pred_ET_byage.RDS")

# ------------------------------------------------------------------------------
# ------------------- MODEL 2: Predicting fertilization rate  ------------------
# ------------------------------------------------------------------------------

dat_step2_0 <- sart_stim %>%
  left_join(set_assign, by="external_patient_id") %>%
  filter(num_retrieved != 0) %>%
  filter(total2pn <= num_retrieved) 

# because need to use sample size in computation of t1, split up now
dat_step2_train0 <- dat_step2_0 %>%
  filter(set=="Training Set")

n_fert_train <- nrow(dat_step2_train0)

dat_step2_test0 <- dat_step2_0 %>%
  filter(set=="Testing Set")

n_fert_test <- nrow(dat_step2_test0)

dat_step2_train <- dat_step2_train0 %>%
   mutate(fert_rate_proxy = total2pn/num_retrieved
         , t1 = (fert_rate_proxy*(n_fert_train-1) + 0.5)/n_fert_train
         , t2 = log(t1/(1-t1))) %>%
  select(-pregnancy_outcome_start_date_diff)


dat_step2_test <- dat_step2_test0 %>%
  mutate(fert_rate_proxy = total2pn/num_retrieved
         , t1 = (fert_rate_proxy*(n_fert_test-1) + 0.5)/n_fert_test
         , t2 = log(t1/(1-t1))) %>%
  select(-pregnancy_outcome_start_date_diff)         


# ----------- SPLIT INTO DIFF DATASETS FOR DIFF AGE GROUPS --------------------

# TRAINING DATA
dat_step2_train_a1 <- dat_step2_train %>%
  filter(agegrp=="<35")

dat_step2_train_a2 <- dat_step2_train %>%
  filter(agegrp=="35-37")

dat_step2_train_a3 <- dat_step2_train %>%
  filter(agegrp=="38-40")

dat_step2_train_a4 <- dat_step2_train %>%
  filter(agegrp=="41-42")

dat_step2_train_a5 <- dat_step2_train %>%
  filter(agegrp==">42")

# TESTING DATA
dat_step2_test_a1 <- dat_step2_test %>%
  filter(agegrp=="<35")

dat_step2_test_a2 <- dat_step2_test %>%
  filter(agegrp=="35-37")

dat_step2_test_a3 <- dat_step2_test %>%
  filter(agegrp=="38-40")

dat_step2_test_a4 <- dat_step2_test %>%
  filter(agegrp=="41-42")

dat_step2_test_a5 <- dat_step2_test %>%
  filter(agegrp==">42")

# ----------- FIT THE MODELS --------------------

fit_fert_mods <- function(dat){
  
  if(deparse(substitute(dat))=="dat_step2_train_a45"){
    predictors_a <- predictors2
    
  }
  else{
    predictors_a <- predictors1
  }
  
  print(predictors_a)
  
  # can't use directly on "fert_rate_proxy" because some obs have value of 1
  # and all observations must be in (0, 1) for beta reg
  mod2_bML <- betareg(as.formula(str_replace_all(paste0("t1 ~ ", predictors_a), "\n", ""))
                      , data=dat
                      , type = "ML")
  
  # linear regression
  mod2_lm <- lm(as.formula(str_replace_all(paste0("fert_rate_proxy ~ ", predictors_a), "\n", ""))
                , data = dat)
  
  # transformation to whole real line (t2 not constrained to (0,1))
  mod2_lmt <- lm(as.formula(str_replace_all(paste0("t2 ~ ", predictors_a), "\n", ""))
                 , data = dat)
  
  # quantile regression model
  mod2_qr <- rq(as.formula(str_replace_all(paste0("fert_rate_proxy ~ ", predictors_a), "\n", ""))
                , data = dat
                , tau = c(0.10, 0.25, 0.50))
  
  # CART
  mod2_cart <- rpart(as.formula(str_replace_all(paste0("fert_rate_proxy ~ ", predictors_a), "\n", ""))
                     , data = dat
                     , method = "anova"   # "class" for a classification tree and "anova" for a regression tree
                     , control=rpart.control(minsplit=2, minbucket=1, cp=0.001)
  )
  
  pred_step2_train <- dat %>%
    mutate(set="Training"
           , mod2_bML_pred = predict(mod2_bML, newdata=dat)
           , mod2_lm_pred = predict(mod2_lm, newdata=dat)
           , mod2_lmt_pred00 = predict(mod2_lmt, newdata=dat)
           , mod2_lmt_pred0 = exp(mod2_lmt_pred00)/(1+exp(mod2_lmt_pred00))
           , mod2_lmt_pred = (mod2_lmt_pred0*n_fert_train - 0.5)/(n_fert_train-1)
           , mod2_qr10_pred = predict(mod2_qr, newdata=dat)[,1]
           , mod2_qr25_pred = predict(mod2_qr, newdata=dat)[,2]
           , mod2_qr50_pred = predict(mod2_qr, newdata=dat)[,3]
           , mod2_cart_pred = predict(mod2_cart, newdata=dat, type="vector")
           ##### differences
           , across(starts_with("mod2_"), ~ . - fert_rate_proxy, .names = "{.col}_diff")
           ##### directions
           , across(ends_with("_diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                    , . < 0 ~ "Underestimate"
                                                    , . > 0 ~ "Overestimate")
                    , .names = "{.col}_dir"))
  
  # ----------------------- GET PREDICTIONS ON TEST SET ----------------------- 
  
  assign("test_dat", eval(as.name(str_replace(deparse(substitute(dat)), "train", "test"))))
  
  pred_step2_test <- test_dat %>%
    mutate(set="Testing"
           , mod2_bML_pred = predict(mod2_bML, newdata=test_dat)
           , mod2_lm_pred = predict(mod2_lm, newdata=test_dat)
           , mod2_lmt_pred00 = predict(mod2_lmt, newdata=test_dat)
           , mod2_lmt_pred0 = exp(mod2_lmt_pred00)/(1+exp(mod2_lmt_pred00))
           , mod2_lmt_pred = (mod2_lmt_pred0*n_fert_train - 0.5)/(n_fert_train-1)
           , mod2_qr10_pred = predict(mod2_qr, newdata=test_dat)[,1]
           , mod2_qr25_pred = predict(mod2_qr, newdata=test_dat)[,2]
           , mod2_qr50_pred = predict(mod2_qr, newdata=test_dat)[,3]
           , mod2_cart_pred = predict(mod2_cart, newdata=test_dat, type="vector")
           ##### differences
           , across(starts_with("mod2_"), ~ . - fert_rate_proxy, .names = "{.col}_diff")
           ##### directions
           , across(ends_with("_diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                    , . < 0 ~ "Underestimate"
                                                    , . > 0 ~ "Overestimate")
                    , .names = "{.col}_dir"))
  
  # ----------------------- COMBINE PREDICTIONS ----------------------- 
  
  pred_dat <- bind_rows(pred_step2_train, pred_step2_test)
  
  return(pred_dat)
}

# ANALYSES STRATIFIED BY AGE
pred_step2_a1 <- fit_fert_mods(dat=dat_step2_train_a1)    
pred_step2_a2 <- fit_fert_mods(dat=dat_step2_train_a2)   
pred_step2_a3 <- fit_fert_mods(dat=dat_step2_train_a3)   
#pred_step2_a4 <- fit_fert_mods(dat=dat_step2_train_a4)    
#pred_step2_a5 <- fit_fert_mods(dat=dat_step2_train_a5)   

# collapse last two groups
dat_step2_train_a45 <- bind_rows(dat_step2_train_a4, dat_step2_train_a5)
dat_step2_test_a45 <- bind_rows(dat_step2_test_a4, dat_step2_test_a5)

pred_step2_a45 <- fit_fert_mods(dat=dat_step2_train_a45)


# ---------------------- SAVE PREDICTIONS FROM STEP 2 -------------------------

quantile(dat_step2_train$fert_rate_proxy, p=c(0.10, 0.25, 0.50))

pred_fert_byage <- bind_rows(pred_step2_a1, pred_step2_a2, pred_step2_a3
                             , pred_step2_a45) %>%
  mutate(perc10 = 0.20
         , perc25 = 0.40
         , perc50 = 0.58
         , across(c("perc10", "perc25", "perc50"), ~ . - fert_rate_proxy, .names = "{.col}_diff")
         , across(c("perc10_diff", "perc25_diff", "perc50_diff"), ~ case_when(.== 0 ~ "Exactly equal"
                                                                              , . < 0 ~ "Underestimate"
                                                                              , . > 0 ~ "Overestimate")
                  , .names = "{.col}_dir"))

saveRDS(pred_fert_byage, "Z:/R03/Data/pred_fert_byage.RDS")
