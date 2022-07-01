# FINAL MODELS
# June 2022

library(tidyverse)
library(survival)
library(ranger)     

sart_stim <- readRDS("Z:/R03/Data/sart_stim.RDS")

sart_transfers <- readRDS("Z:/R03/Data/sart_transfers.RDS") 

set_assign <- readRDS("Z:/R03/Data/set_assign.RDS")

predictors_slim <- "agegrp6 + clinic_state"

predictors_mid <- "agegrp6 + clinic_state + amhcat3_new + diminished_ovarian_reserve + num_retrieved" 

predictors_full <- "agegrp6 + clinic_state + amhcat3_new + num_retrieved + bmicat6 + 
                    gravidity_cat4 + parity_cat4 + fsh_gt10 +
                    male_infertility + dx_tubal + endometriosis + uterine + dx_ovulation + 
                    diminished_ovarian_reserve + unexplained"


# ------------------------------------------------------------------------------
# -----------         STAGE 1: DAY 3 OR DAY 5 ET             -------------------
# ------------------------------------------------------------------------------

sart_transfers_set <- sart_transfers %>%
  left_join(set_assign, by="external_patient_id") %>%
  mutate(d5 = ifelse(d5_include==1, yes=1, no=0)) %>%
  filter(!is.na(d5))

sart_transfers_set %>% count(set)

dat_dayET_train <- sart_transfers_set %>%
  filter(set=="Training Set") 

# ------------------------------- FIT & SAVE FINAL MODEL --------------------- #

# prevent the model.frame, y, and model.matrix from being returned in model object
mod_dayET <- glm(as.formula(str_replace_all(paste0("d5 ~  ", predictors_mid), "\n", ""))
                 , data=dat_dayET_train, family="binomial"
                 , model = FALSE, y = FALSE, x = FALSE) 

# remove data from model object
mod_dayET$data <- NULL

save(mod_dayET, file="Z:/R03/Data/mod_dayET.Rdata")



# ------------------------------------------------------------------------------
# -----------         STAGE 2 - MODEL 1: PROP BLAST          -------------------
# ------------------------------------------------------------------------------

dat_blast <- sart_stim %>%
  filter(d5_include==1) %>%
  left_join(set_assign, by="external_patient_id") %>%
  filter(num_retrieved != 0 & reason_for_no_transfer!="Inability to obtain sperm" & reason_for_no_transfer != "PGT") %>%
  mutate(num_TC = case_when(transfer_attempted=="N" & (num_transferred==0 | is.na(num_transferred)) ~ num_fresh_embryos_cryoed
                            , TRUE ~ num_TC)
         , num_TC0 = case_when(is.na(num_fresh_embryos_cryoed) & transfer_attempted=="N" ~ NA_real_
                               , is.na(num_fresh_embryos_cryoed) ~ num_transferred
                               , is.na(num_transferred) ~ num_fresh_embryos_cryoed
                               , TRUE ~ num_TC)
         , blast_rate = num_TC/num_retrieved
         , blast_rate0 = num_TC0/num_retrieved) %>%
  # remove 3 cycles where num_TC0 > num_retrieved (one or the other has a data entry error...)
  filter(num_TC0 <= num_retrieved) %>%
  select(blast_rate, num_TC, blast_rate0, num_TC0, num_blasts, num_transferred, num_fresh_embryos_cryoed
         , d5_include, set, reporting_year, external_patient_id, external_cycle_id
         , agegrp6, clinic_state, amhcat3, num_retrieved, bmicat6, gravidity_cat4, parity_cat4 
         , fsh_gt10, male_infertility, dx_tubal, endometriosis, uterine, dx_ovulation
         , diminished_ovarian_reserve, unexplained
         , transfer_attempted, reason_for_no_transfer)

dat_blast_train <- dat_blast %>%
  filter(set=="Training Set") 

# ------------------------------- FIT & SAVE FINAL MODEL --------------------- #

mod_d5blast <- ranger(as.formula(str_replace_all(paste0("blast_rate0 ~ ", predictors_mid), "\n", ""))
                      , data=dat_blast_train)

names(mod_d5blast)

save(mod_d5blast, file="Z:/R03/Data/mod_d5blast.Rdata")


# ------------------------------------------------------------------------------
# -----------         STAGE 2 - MODEL 2: NUM ET NEEDED       -------------------
# ------------------------------------------------------------------------------

dat_numET0 <- sart_transfers %>%
  filter(!is.na(d5_include)) %>%
  select(external_patient_id, external_cycle_id, cycle_order, reporting_year, d5_include, num_transferred, lb) %>% 
  arrange(external_patient_id, cycle_order, external_cycle_id) %>%
  group_by(external_patient_id) %>%
  mutate(cycnum = row_number()
         , minD5 = min(d5_include)
         , maxD5 = max(d5_include)
         , totalET = cumsum(num_transferred)
         , totalLB = cumsum(lb)
         , totalETminus1 = totalET - 1
         , Keep = (lb == 1 & totalLB == 1) |
           (row_number() == max(row_number()) & totalLB==0)
         , max_cycle = max(cycle_order)) %>%
  ungroup()

dat_numET1 <- dat_numET0 %>%
  ungroup() %>%
  filter(Keep == TRUE) %>%
  mutate(censored = ifelse(totalLB==1,yes=0,no=1)) %>%
  select(external_patient_id, minD5, maxD5, cycle_order, totcyc = cycnum, totalET, totalETminus1, lb, censored)

# 6.9% are a mix of day 3 and day 5 transfers; 19% are day 3 only; and 74% are day 5 only
dat_numET1 %>% count(minD5, maxD5) %>% mutate(perc=n/sum(n))

first_cyc <- sart_stim %>%
  filter(cycle_order==1) %>%
  select(first_stim_cycle_id = external_cycle_id, external_patient_id
         , stim_reporting_year = reporting_year, stim_cycle_order = cycle_order
         , clinic_state, agegrp6, numER
         , patient_age_at_start, partner_age_at_start_c, partner_age_at_start_missing
         , bmi_c, bmicat4, bmicat6
         , smoker, gravidity_cat4, parity_cat4, max_fsh, fsh_gt10, amh_last_value, amhcat3
         , male_infertility, endometriosis, dx_ovulation, diminished_ovarian_reserve
         , tubal_ligation, tubal_hydrosalpinx, tubal_other, dx_tubal
         , uterine, unexplained
         , num_retrieved
  )

# --------- merge w/ first stim cycle parameters & train/test assignment -------

dat_numET <- dat_numET1 %>%
  left_join(first_cyc, by="external_patient_id") %>%
  left_join(set_assign, by="external_patient_id")

dat_numETd5 <- dat_numET %>%
  filter(maxD5==1) 

dat_numETd5_train <- dat_numETd5 %>%
  filter(set=="Training Set")

# ------------------------------- FIT & SAVE FINAL MODELS --------------------- #

mod_ETd5 <- survreg(as.formula(paste0("Surv(totalET, lb) ~ ", predictors_mid))
                    , data=dat_numETd5_train, dist='loglogistic')

save(mod_ETd5, file="Z:/R03/Data/mod_numETd5.Rdata")

