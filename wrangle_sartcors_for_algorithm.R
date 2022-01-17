# prepare sart cors database for analysis (derive vars, wrangle etc)
# kat correia

library(tidyverse)

# -------------------- load data ----------------------------------------

# SART CORS sends data as csv file; imported that file into R, used
# janitor::clean_names() to clean the variable names, and saved as RDS for
# future faster loading
sart_full_raw <- readRDS("Z:/SARTCORS/sart_full_raw.RDS")

# -------------------- wrangle data ----------------------------------------

# "cannot allocate vector of size ..." errors + R keeps crashing; 
# try breaking up the steps and re-arranging

# derive variables and make sure no GCs or eggs donors are included
sart_full_a <- sart_full_raw %>%
  # exclude non-ivf cycles and cycles where both fresh and frozen were transferred
  filter(ivf=="Y" & gestational_carrier=="N" & donor_oocyte=="N" &
           !(fresh_embryo=="Y" & thawed_embryo=="Y")) %>%
  mutate(partner_age_at_start = parse_number(partner_age_at_start)
         , partner_age_at_start_c = ifelse(partner_age_at_start < 18 | partner_age_at_start > 50
                                          , yes=NA_real_
                                          , no=partner_age_at_start)
         , partner_age_at_start_missing = ifelse(is.na(partner_age_at_start), yes=1, no=0)
         , height_inches = parse_number(height_inches)
         , weight_pounds = parse_number(weight_pounds)
         , bmi = parse_number(bmi)
         , bmi_computed = (weight_pounds/(height_inches^2))*703
         , bmi_c = case_when(bmi > 13 & bmi < 70 ~ bmi
                             , bmi_computed > 13 & bmi_computed < 70 ~ bmi_computed)
         , bmicat4 = factor(case_when(is.na(bmi_c) ~ "Missing"
                               , bmi_c < 18.5 ~ "Underweight"
                               , bmi_c >= 18.5 & bmi_c < 25 ~ "Healthy weight"
                               , bmi_c >= 25 & bmi_c < 30 ~ "Overweight"
                               , bmi_c >= 30 ~ "Obese")
                           , levels = c("Underweight", "Healthy weight"
                                        , "Overweight", "Obese", "Missing"))
         , bmicat6 = factor(case_when(is.na(bmi_c) ~ "Missing"
                               , bmi_c < 18.5 ~ "Underweight"
                               , bmi_c >= 18.5 & bmi_c < 25 ~ "Healthy weight"
                               , bmi_c >= 25 & bmi_c < 30 ~ "Overweight"
                               , bmi_c >= 30 & bmi_c < 35 ~ "Obese, Class I"
                               , bmi_c >= 35 & bmi_c < 40 ~ "Obese, Class II"
                               , bmi_c >= 40 ~ "Obese, Class III")
                            , levels = c("Underweight", "Healthy weight"
                                         , "Overweight", "Obese, Class I"
                                         , "Obese, Class II", "Obese, Class III"
                                         , "Missing"))
         , gravidity_cat4 = factor(case_when(gravidity %in% c("Unknown", "0", "1", "2") ~ gravidity
                                      , TRUE ~ "3+")
                                   , levels=c("0", "1", "2", "3+", "Unknown"))
         , full_term_birthsn = case_when(full_term_births %in% c("Unknown", "Not Entered") ~ NA_real_
                                          , full_term_births == ">10" ~ 11
                                          , TRUE ~ as.numeric(full_term_births))
         , pre_term_birthsn = case_when(pre_term_births %in% c("Unknown", "Not Entered") ~ NA_real_
                                        , pre_term_births == ">10" ~ 11
                                        , TRUE ~ as.numeric(pre_term_births))
         , parity = full_term_birthsn + pre_term_birthsn
         , parity_cat4 = factor(case_when(is.na(parity) ~ "Unknown"
                                  , 0 <= parity & parity <= 2 ~ as.character(parity)
                                  , parity >= 3 ~ "3+")
                                , levels =c("0", "1", "2", "3+", "Unknown"))
         # smoking fields changed in 2016 (see data dictionary)
         , smoke3mo = case_when(reporting_year < 2016 & !(smoke_prior3months %in% c("None", "Unknown")) ~ "Yes"
                                , reporting_year < 2016 & smoke_prior3months=="None" ~ "No"
                                , reporting_year < 2016 & smoke_prior3months=="Unknown" ~ "Unknown"
                                , reporting_year >= 2016 ~ smoke3months_prior_to_start)
         , months_attempting_pregnancy = parse_number(months_attempting_pregnancy)
         , max_fsh = parse_number(max_fsh)
         , fsh_gt10 = factor(case_when(is.na(max_fsh) ~ "Missing"
                                , max_fsh <= 10 ~ "<= 10"
                                , max_fsh > 10 ~ "> 10")
                             , levels=c("<= 10", "> 10", "Missing"))
         , amhcat3 = factor(case_when(is.na(amh_last_value) | amh_last_value == 0 ~ "Missing"
                               , amh_last_value < 1 ~ "< 1"
                               , amh_last_value >= 1 & amh_last_value < 4 ~ "1 - <4"
                               , amh_last_value >= 4 ~ ">= 4")
                            , levels=c("< 1", "1 - <4", ">= 4", "Missing"))
         , clomiphene_dosage = parse_number(clomiphene_dosage)
         , fsh_dosage = parse_number(fsh_dosage)
         , aromatase_inhibitors_dosage = parse_number(aromatase_inhibitors_dosage)
         , letrozole_dosage = parse_number(letrozole_dosage)
         , dx_tubal = ifelse(tubal_ligation=="Y" | tubal_hydrosalpinx=="Y"
                           | tubal_other=="Y", yes="Y", no="N")
         , any_icsi = case_when(icsi_none=="Y" ~ "No"
                                , icsi_some_oocytes=="Y" | icsi_all_mature_oocytes=="Y" ~ "Yes"
                                , icsi_not_entered=="Y" | icsi_unknown=="Y" ~ "Missing")
         , any_ah = case_when(assisted_hatching=="None" ~ "No"
                              , assisted_hatching %in% c("All Transferred Embryos", "Some Embryos") ~ "Yes"
                              , assisted_hatching %in% c("Not Entered", "Unknown") ~ "Missing")
         , total2pn = parse_number(total2pn)
         , cycle_type = case_when(retrieval_type_autologous_retrieval1=="Fresh" ~ "Fresh"
                                  , retrieval_type_autologous_retrieval1=="Embryo Banking" ~ "Embryo Banking"
                                  , retrieval_type_autologous_retrieval1=="NULL" ~ "Cancelled"
                                  , TRUE ~ "Frozen")
         , preg = case_when(treatment_outcome_clinical_intrauterine_gestation=="Y" ~ 1
                              , treatment_outcome_clinical_intrauterine_gestation=="N" ~ 0)
         , lb = case_when(pregnancy_outcome_live_birth=="Y" ~ 1
                         , pregnancy_outcome_live_birth=="N" ~ 0)
         , clinic_state = ifelse(i_clinic_state_restricted=="NULL", yes="Other", no=i_clinic_state_restricted)
         , clinic_state2 = ifelse(clinic_state=="Other"
                                  , yes = paste0("Other - ", clinic_region_usa)
                                  , no = clinic_state)) %>%
  # only keep relevant variables 
  select(reporting_year, external_patient_id, cycle_order, external_cycle_id
         , starts_with("linked_"), patient_age_at_start
         , partner_age_at_start, partner_age_at_start_c, partner_age_at_start_missing
         , i_clinic_state_restricted, clinic_state, clinic_state2, clinic_region_usa
         , height_inches, weight_pounds, bmi, bmi_c, bmicat4, bmicat6
         , smoker, smoke_prior3months, smoke3months_prior_to_start, smoke3mo
         , gravidity, gravidity_cat4, full_term_births, pre_term_births, parity, parity_cat4
         , prior_gonadotropin_cycles, prior_fresh_cycles, prior_frozen_cycles
         , max_fsh, fsh_gt10, amh_last_value, amhcat3, male_infertility, prior_vasectomy
         , endometriosis, polycystic_ovaries, diminished_ovarian_reserve
         , tubal_ligation, tubal_hydrosalpinx, tubal_other, dx_tubal
         , uterine, unexplained
         , other_rfa, other_pgd 
         , clomiphene, clomiphene_dosage, fsh, fsh_dosage, aromatase_inhibitors
         , aromatase_inhibitors_dosage, minimal_stimulation, unstimulated
         , medication_with_lhhcg_activity, letrozole, letrozole_dosage
         , agonist_suppression, agonist_flare, antagonist_suppression
         , cycle_cancelled, day_of_cancellation, starts_with("cancel_reason")
         , starts_with("comp_")
         , sperm_source_not_entered, sperm_source_partner, sperm_source_donor
         , sperm_source_mixed, sperm_source_unknown, sperm_status
         , starts_with("retrieval_type_autologous_retrieval"), cycle_type
         , starts_with("number_retrieved_autologous_retrieval")
         , starts_with("number_thawed_autologous_retrieval")
         , total2pn, assisted_hatching, any_ah, pgd, starts_with("pgd_")
         , starts_with("icsi_"), any_icsi
         , laboratory_routine, endometrial_thickness
         , starts_with("embryos_to_uterus"), num_suitable_for_transfer
         , transfer_attempted, reason_for_no_transfer, no_transfer_other_text
         , elective_single_embryo_transfer
         , starts_with("treatment_outcome"), starts_with("pregnancy_outcome")
         , pregnancy_loss_abortion
         , starts_with("grade_embryo"), starts_with("stage_embryo")
         , number_live_born
         , preg, lb)

# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_full_a")) { rm(sart_full_raw) }

nrow(sart_full_a)

count(sart_full_a, clinic_region_usa, clinic_state)
mosaic::favstats(~patient_age_at_start, data=sart_full_a)

# limit to women ages 18-45 
sart_full_b <- sart_full_a %>%
  filter(patient_age_at_start >= 18 & patient_age_at_start <= 45)

# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_full_b")) { rm(sart_full_a) }

# nrow(sart_full_b)
# mosaic::favstats(~patient_age_at_start, data=sart_full_b)
# mosaic::favstats(~partner_age_at_start, data=sart_full_b)
# mosaic::favstats(~partner_age_at_start_c, data=sart_full_b)
# count(sart_full_b, partner_age_at_start_missing)

# now compute number of oocytes retrieved and number of embryos transferred in two data steps: 
# first, make all the fields numeric
# then, sum them up
sart_full_c <- sart_full_b %>%
  mutate(across(starts_with(c("embryos_to_uterus","number_retrieved_autologous_retrieval")), ~ parse_number(.x)))
 
# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_full_c")) { rm(sart_full_b) }

sart <- sart_full_c %>%
  mutate(num_retrieved = rowSums(.[grep("number_retrieved_autologous_retrieval", names(.))], na.rm = TRUE)
         , num_transferred = rowSums(.[grep("embryos_to_uterus", names(.))], na.rm = TRUE))


# TO SAVE SPACE, remove each dataset once have next
if (exists("sart")) { rm(sart_full_c) }

# check derivations
# mosaic::favstats(~num_retrieved, data=sart)
# mosaic::favstats(~num_transferred, data=sart)
# sart %>% count(bmicat4)
# sart %>% count(bmicat6)
# sart %>% count(fsh_gt10)
# sart %>% count(amhcat3)
# mosaic::favstats(amh_last_value ~ amhcat3, data=sart)
# mosaic::favstats(max_fsh ~ fsh_gt10, data=sart)
# sart %>% count(gravidity, gravidity_cat4)
# sart %>% count(parity, parity_cat4)
# sart %>% count(retrieval_type_autologous_retrieval1, cycle_type)
# sart %>% count(cycle_type, num_transferred)
# sart %>% count(elective_single_embryo_transfer, num_transferred)
# sart %>% count(dx_tubal, tubal_ligation, tubal_hydrosalpinx, tubal_other)
# sart %>% count(pgd)
# sart %>% count(transfer_attempted, pgd)

# sart %>%
#   select(cycle_type, num_transferred, starts_with("embryos_to_uterus")) %>%
#   slice(1:20)

# -------------------- save cleaned data -------------------------------------

saveRDS(sart, file = "Z:/SARTCORS/sart.RDS")

