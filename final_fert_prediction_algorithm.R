# FINAL MODELS
# March 14, 2022

sart_stim <- readRDS("Z:/R03/Data/sart_stim.RDS")

sart_transfers <- readRDS("Z:/R03/Data/sart_transfers.RDS") %>%
  mutate(external_cycle_id = ifelse(is.na(external_cycle_id), yes=frozen_cycle_id
                                    , no=external_cycle_id))

set_assign <- readRDS("Z:/R03/Data/set_assign.RDS")

predictors1 <- "clinic_state + bmicat6 + smoker + gravidity_cat4 + parity_cat4 + fsh_gt10 + amhcat3
               + male_infertility + endometriosis + dx_ovulation + diminished_ovarian_reserve
               + dx_tubal + uterine + unexplained + numER"

predictors2 <- paste("agegrp +", predictors1)

# ------------------------------------------------------------------------------
# -----------                   ET MODEL                    --------------------
# ------------------------------------------------------------------------------

# ------------------- get transfers --------------------------------------------

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

# --------- merge w/ first stim cycle parameters & train/test assignment -------

dat_step1 <- dat_step1_1 %>%
  left_join(first_cyc, by="external_patient_id") %>%
  left_join(set_assign, by="external_patient_id")

dat_step1_train <- dat_step1 %>%
  filter(set=="Training Set")

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

# ------------------------- fit models ---------------------------------------

fit_ET_mod <- function(dat){
  
  print(mosaic::tally(~censored, data=dat, format="percent"))
  
  if(deparse(substitute(dat))=="dat_step1_train_a45"){
    predictors_a <- predictors2
    
  }
  else{
    predictors_a <- predictors1
  }
  
  print(predictors_a)
  
  # quantile regression model
  mod1_qr <- rq(as.formula(str_replace_all(paste0("totalET ~ ", predictors_a), "\n", ""))
                , data = dat
                , tau = 0.5)
  
  # remove data from the model object 
  mod1_qr$model <- NULL
  mod1_qr$x <- NULL
  mod1_qr$y <- NULL 

  print(names(mod1_qr))
  
  return(mod1_qr)
}

# ANALYSES STRATIFIED BY AGE
mod_ET_a1 <- fit_ET_mod(dat=dat_step1_train_a1)    # 29% censored
mod_ET_a2 <- fit_ET_mod(dat=dat_step1_train_a2)    # 40.6% censored
mod_ET_a3 <- fit_ET_mod(dat=dat_step1_train_a3)    # 57% censored

# try collapsing last two groups
dat_step1_train_a45 <- bind_rows(dat_step1_train_a4, dat_step1_train_a5)

mod_ET_a45 <- fit_ET_mod(dat=dat_step1_train_a45)    # 81% censored


# ------------------------------------------------------------------------------
# -----------                   FERT MODEL                  --------------------
# ------------------------------------------------------------------------------

dat_step2_0 <- sart_stim %>%
  left_join(set_assign, by="external_patient_id") %>%
  filter(num_retrieved != 0) %>%
  filter(total2pn <= num_retrieved) 

# because need to use sample size in computation of t1, split up now
dat_step2_train0 <- dat_step2_0 %>%
  filter(set=="Training Set")

n_fert_train <- nrow(dat_step2_train0)

dat_step2_train <- dat_step2_train0 %>%
  mutate(fert_rate_proxy = total2pn/num_retrieved) %>%
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


# ----------- FIT THE MODELS --------------------

fit_fert_mod <- function(dat){
  
  if(deparse(substitute(dat))=="dat_step2_train_a45"){
    predictors_a <- predictors2
    
  }
  else{
    predictors_a <- predictors1
  }
  
  print(predictors_a)
  
  # CART
  mod2_cart <- rpart(as.formula(str_replace_all(paste0("fert_rate_proxy ~ ", predictors_a), "\n", ""))
                     , data = dat
                     , method = "anova"   # "class" for a classification tree and "anova" for a regression tree
                     , control=rpart.control(minsplit=2, minbucket=1, cp=0.001)
  )
  
  # remove data from the model object 
  mod2_cart$y <- NULL 
  
  print(names(mod2_cart))
  
  return(mod2_cart)
}

# ANALYSES STRATIFIED BY AGE
mod_fert_a1 <- fit_fert_mod(dat=dat_step2_train_a1)    
mod_fert_a2 <- fit_fert_mod(dat=dat_step2_train_a2)   
mod_fert_a3 <- fit_fert_mod(dat=dat_step2_train_a3)   
#mod_fert_a4 <- fit_fert_mods(dat=dat_step2_train_a4)    
#mod_fert_a5 <- fit_fert_mods(dat=dat_step2_train_a5)   

# try collapsing last two groups
dat_step2_train_a45 <- bind_rows(dat_step2_train_a4, dat_step2_train_a5)

mod_fert_a45 <- fit_fert_mod(dat=dat_step2_train_a45)


# ------------------------------------------------------------------------------
# -----------                SAVE MODELS                    --------------------
# ------------------------------------------------------------------------------

save(mod_ET_a1, mod_ET_a2, mod_ET_a3, mod_ET_a45
    , mod_fert_a1, mod_fert_a2, mod_fert_a3, mod_fert_a45
    , file="Z:/R03/Data/final_models.RData")
