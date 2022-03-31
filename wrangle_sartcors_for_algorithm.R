# kat correia - feb 2022
# prepare sart cors database for analysis (derive vars, wrangle etc)
# and assign women to training and testing sets
# only include women with first stim cycle between 1/1/14 and 12/31/19
# that is, exclude transfers from women that had stim cycle prior  to 1/1/13

library(tidyverse)

# -------------------- convert csv file to RDS --------------------------------
# -------------------- ONLY DO THIS ONCE --------------------------------------

# sart_raw0 <- read_csv("Z:/R03/Data/FB94219_20220114.csv")
# 
# sart_raw <- sart_raw0 %>%
#   janitor::clean_names()
# 
# saveRDS(sart_raw, "Z:/R03/Data/sart_raw.RDS")
# 
# sart_raw_metainfo <- data.frame(orig_names = names(sart_raw0)
#                          , type = sapply(sart_raw0, class)
#                          , new_names = names(sart_raw))
# 
# write_csv(sart_raw_metainfo, "Z:/R03/Data/sart_raw_metainfo.csv")


# -------------------- load data ----------------------------------------

sart_full_raw <- readRDS("Z:/SARTCORS/sart_full_raw.RDS")

# -------------------- wrangle data ----------------------------------------

# "cannot allocate vector of size ..." errors + R keeps crashing; 
# try breaking up the steps and re-arranging

# derive variables and make sure no GCs or eggs donors are included
sart_a <- sart_raw %>%
  # exclude non-ivf cycles and cycles where both fresh and frozen were transferred
  filter(pgd=="N" & gift=="N" & zift=="N" & gestational_carrier=="N" & donor_oocyte=="N" &
           donor_embryo=="N") %>%
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
         , any_prior_sab = case_when(spontaneous_abortions == "Not Entered" & gravidity=="0" ~ 0
                                     , spontaneous_abortions=="Unknown" ~ NA_real_
                                     , gravidity=="Unknown" & spontaneous_abortions=="Not Entered" ~ NA_real_
                                     , spontaneous_abortions=="0" ~ 0
                                     , TRUE ~ 1)
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
         , fsh_dosage = parse_number(fsh_dosage)
         , dx_tubal = ifelse(tubal_ligation=="Y" | tubal_hydrosalpinx=="Y"
                             | tubal_other=="Y", yes="Y", no="N")
         , dx_ovulation = ifelse(polycystic_ovaries=="Y" | ovulation_disorders == "Y" | other_ovulation_disorders == "Y"
                                 | hypothalamic_amenorrhea=="Y", yes="Y", no="N")
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
         #, endometrial_thickness = parse_number(endometrial_thickness)
         , preg = case_when(reason_for_no_transfer %in% c("PGT","Planned multiple cycles (Egg/Embryo Accumulation) with PGT"
                                                          , "Planned multiple cycles (Egg/Embryo Accumulation) without PGT"
                                                          , "Intent for Freeze-All, Unspecific"
                                                          , "Thaw for Refreeze, No Intent to Transfer") ~ NA_real_
                            , treatment_outcome_clinical_intrauterine_gestation=="Y" ~ 1
                            , treatment_outcome_clinical_intrauterine_gestation=="N" ~ 0)
         , lb = case_when(reason_for_no_transfer %in% c("PGT","Planned multiple cycles (Egg/Embryo Accumulation) with PGT"
                                                        , "Planned multiple cycles (Egg/Embryo Accumulation) without PGT"
                                                        , "Intent for Freeze-All, Unspecific"
                                                        , "Thaw for Refreeze, No Intent to Transfer") ~ NA_real_
                          , pregnancy_outcome_live_birth=="Y" ~ 1
                          , pregnancy_outcome_live_birth=="N" ~ 0)
         , clinic_state = ifelse(clinic_state_restricted=="NULL", yes="Other", no=clinic_state_restricted)) %>%
  # only keep relevant variables 
  select(reporting_year, external_patient_id, cycle_order, external_cycle_id
         , starts_with("linked_"), patient_age_at_start
         , partner_age_at_start, partner_age_at_start_c, partner_age_at_start_missing
         , clinic_state_restricted, clinic_state
         , starts_with("pt_"), starts_with("part_")
         , height_inches, weight_pounds, bmi, bmi_c, bmicat4, bmicat6
         , smoker, smoke_prior3months, smoke3months_prior_to_start, smoke3mo
         , gravidity, gravidity_cat4, full_term_births, pre_term_births, parity, parity_cat4
         , spontaneous_abortions, any_prior_sab, prior_fresh_cycles, prior_frozen_cycles
         , max_fsh, fsh_gt10, amh_last_value, amhcat3, male_infertility, prior_vasectomy
         , endometriosis, diminished_ovarian_reserve
         , tubal_ligation, tubal_hydrosalpinx, tubal_other, dx_tubal
         , polycystic_ovaries, hypothalamic_amenorrhea, ovulation_disorders, other_ovulation_disorders, dx_ovulation
         , uterine, unexplained, recurrent_pregnancy_loss
         , other_rfa, other_pgd 
         , abnormal_semen_parameters, abnormal_sperm_parameters, azoospermia_obstructive
         , azoospermia_nonobstructive, oligospermia_moderate, oligospermia_severe
         , low_motility, low_morphology, very_severe_male_factor
         , clomiphene, fsh, fsh_dosage, aromatase_inhibitors, minimal_stimulation, unstimulated
         , medication_with_lhhcg_activity, letrozole, ovulation_trigger
         , agonist_suppression, agonist_flare, antagonist_suppression
         , embryo_banking, oocyte_banking, patient_freeze_all
         , cycle_cancelled, starts_with("cancel_reason")
         #, starts_with("comp_")
         , sperm_source_not_entered, sperm_source_partner, sperm_source_donor
         , sperm_source_mixed, sperm_source_unknown, sperm_status
         , starts_with("retrieval_type_autologous_retrieval"), cycle_type
         , starts_with("number_retrieved_autologous_retrieval")
         , starts_with("number_thawed_autologous_retrieval")
         , total2pn, assisted_hatching, any_ah, pgd, starts_with("pgd_")
         , starts_with("icsi_"), any_icsi
         , laboratory_routine
         , starts_with("embryos_to_uterus")
         , transfer_attempted, reason_for_no_transfer
         , elective_single_embryo_transfer
         , starts_with("fresh_embryos_cryoed"), starts_with("blastocyst_transfer")
         # embryo morphology
         # can't use these since they are about the quality of the embryos transferred to the uterus
         # , starts_with("grade_embryo_morphology"), starts_with("stage_embryo_morphology")
         # , starts_with("fragmentation_embryo_morphology")
         # , starts_with("morphology_embryo_morphology")
         # , starts_with("symmetry_embryo_morphology")
         # , starts_with("inner_cell_embryo_morphology")
         # , starts_with("trophoblast_embryo_morphology")
         , starts_with("treatment_outcome"), starts_with("pregnancy_outcome")
         , pregnancy_loss_abortion
         , starts_with("grade_embryo"), starts_with("stage_embryo")
         , number_live_born
         , preg, lb, cycle_classification)

# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_full_a")) { rm(sart_full_raw) }

nrow(sart_full_a)

count(sart_full_a, clinic_region_usa, clinic_state)
mosaic::favstats(~patient_age_at_start, data=sart_full_a)

# now compute number of oocytes retrieved and 
# number of embryos transferred in two data steps: 
# first, make all the fields numeric
# then, sum them up
sart_b <- sart_a %>%
  mutate(across(starts_with(c("embryos_to_uterus","number_retrieved_autologous_retrieval")), ~ parse_number(.x)))

# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_b")) { rm(sart_a) }

sart_for_R03 <- sart_b %>%
  mutate(num_retrieved = rowSums(.[grep("number_retrieved_autologous_retrieval", names(.))], na.rm = TRUE)
         , num_transferred = rowSums(.[grep("embryos_to_uterus", names(.))], na.rm = TRUE)
         , agegrp =  case_when(patient_age_at_start < 35 ~ "<35"
                               , patient_age_at_start >= 35 & patient_age_at_start < 38 ~ "35-37"
                               , patient_age_at_start >= 38 & patient_age_at_start < 41 ~ "38-40"
                               , patient_age_at_start >= 41 & patient_age_at_start < 43 ~ "41-42"
                               , patient_age_at_start >= 43 ~ ">42")
         , numER = case_when(0 <= num_retrieved & num_retrieved <= 5 ~ "<6"
                             , 5 < num_retrieved & num_retrieved <= 10 ~ "6-10"
                             , 10 < num_retrieved & num_retrieved <= 15 ~ "11-15"
                             , 15 < num_retrieved & num_retrieved <= 20 ~ "16-20"
                             , 20 < num_retrieved & num_retrieved <= 30 ~ "21-30"
                             , 30 < num_retrieved & num_retrieved <= 40 ~ "31-40"
                             , 40 < num_retrieved & num_retrieved <= 50 ~ "41-50"
                             , num_retrieved > 50 ~ ">50"))


# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_for_R03")) { rm(sart_b) }

# check derivations
# mosaic::favstats(~num_retrieved, data=sart_for_R03)
# mosaic::favstats(~num_transferred, data=sart_for_R03)
# sart_for_R03 %>% count(bmicat4)
# sart_for_R03 %>% count(bmicat6)
# sart_for_R03 %>% count(fsh_gt10)
# sart_for_R03 %>% count(amhcat3)
# mosaic::favstats(amh_last_value ~ amhcat3, data=sart_for_R03)
# mosaic::favstats(max_fsh ~ fsh_gt10, data=sart_for_R03)
# sart_for_R03 %>% count(gravidity_cat4, gravidity)
# sart_for_R03 %>% count(parity_cat4, parity)
# sart_for_R03 %>% count(retrieval_type_autologous_retrieval1, cycle_type)
# sart_for_R03 %>% count(cycle_type, num_transferred)
# sart_for_R03 %>% count(num_transferred, elective_single_embryo_transfer)
# sart_for_R03 %>% count(dx_tubal, tubal_ligation, tubal_hydrosalpinx, tubal_other)
# sart_for_R03 %>% count(dx_ovulation, polycystic_ovaries, hypothalamic_amenorrhea, ovulation_disorders, other_ovulation_disorders)
# sart_for_R03 %>% count(pgd)
# sart_for_R03 %>% count(transfer_attempted, pgd)

# sart %>%
#   select(cycle_type, num_transferred, starts_with("embryos_to_uterus")) %>%
#   slice(1:20)

# -------------------- save cleaned data -------------------------------------

saveRDS(sart_for_R03, file = "Z:/R03/Data/sart_for_R03.RDS")


# ----------------------------------------------------------------------------
# ------------------ create stim dataset + transfer dataset ------------------
# ----------------------------------------------------------------------------

# -------------------- stimulation cycles -------------------------------------

# if a cycle is a thaw cycle, it will have a non-NULL link_source_external_cycle_id
# if fresh, that field will be NULL
sart_stim0 <- sart_for_R03 %>%
  filter(linked_source_external_cycle_i_ds=="NULL" & reporting_year >= 2014)

sart_stim0 %>% count(reporting_year)
sart_stim0 %>% count(cycle_order)
length(unique(sart_stim0$external_patient_id))

# identify women who had a first cycle 2014 or later 
# and were between 18-45 yo at that first cycle
sart_stim_first <- sart_stim0 %>%
  group_by(external_patient_id) %>%
  mutate(first_cycle = min(cycle_order)) %>%
  filter(first_cycle==1 & cycle_order==first_cycle & 
           patient_age_at_start >= 18 & patient_age_at_start <= 45)

mosaic::favstats(~patient_age_at_start, data=sart_stim_first)

# only keep women who had a first cycle 2014 or later 
# and were between 18-45 yo at that first cycle
sart_stim <- sart_stim0 %>%
  filter(external_patient_id %in% sart_stim_first$external_patient_id)

length(unique(sart_stim$external_patient_id))

# -------------------- transfer cycles (fresh + frozen) -----------------------

# ONLY KEEP TRANSFERS THAT RESULTED FROM STIMS IN 2014+

sart_frozen <- sart_for_R03 %>%
  filter(linked_source_external_cycle_i_ds != "NULL") %>%
  rename(frozen_cycle_id = external_cycle_id) %>%
  mutate(frozen = 1
         # hand-correct one super long weird linked ID that repeats
         , linked_source_external_cycle_i_ds=ifelse(linked_source_external_cycle_i_ds ==
                                                      "76162769, 76162770, 76162771, 76162771, 76162772, 76162774, 76162774, 76162769, 76162769, 76162770, 76162770, 76162771, 76162772, 76162772, 76162774"
                                                    , yes = "76162769, 76162770, 76162771, 76162772, 76162774"
                                                    , no = linked_source_external_cycle_i_ds)) %>%
  # the ones with 5 and 6 end up repeating the linked stim cycles...
  separate(linked_source_external_cycle_i_ds, into=paste0("link",1:6)
           , remove=FALSE, convert=TRUE) %>%
  left_join(select(sart_for_R03, external_cycle_id, stim_year1 = reporting_year)
            , by=c("link1"="external_cycle_id")) %>%
  left_join(select(sart_for_R03, external_cycle_id, stim_year2 = reporting_year)
            , by=c("link2"="external_cycle_id")) %>%
  left_join(select(sart_for_R03, external_cycle_id, stim_year3 = reporting_year)
            , by=c("link3"="external_cycle_id")) %>%
  left_join(select(sart_for_R03, external_cycle_id, stim_year4 = reporting_year)
            , by=c("link4"="external_cycle_id")) %>%
  left_join(select(sart_for_R03, external_cycle_id, stim_year5 = reporting_year)
            , by=c("link5"="external_cycle_id")) %>%
  left_join(select(sart_for_R03, external_cycle_id, stim_year6 = reporting_year)
            , by=c("link6"="external_cycle_id")) %>%
  # only keep frozen transfers w/ embryos produced from stim cycle in 2014 or later
  filter((!is.na(link1) & stim_year1 >= 2014) 
         & (stim_year2 >= 2014 | (is.na(link2) & is.na(stim_year2))) 
         & (stim_year3 >= 2014 | (is.na(link3) & is.na(stim_year3)))
         & (stim_year4 >= 2014 | (is.na(link4) & is.na(stim_year4)))
         & (stim_year5 >= 2014 | (is.na(link5) & is.na(stim_year5)))
         & (stim_year6 >= 2014 | (is.na(link6) & is.na(stim_year6)))) %>%
  filter(external_patient_id %in% sart_stim$external_patient_id) %>%
  filter(num_transferred != 0)


# (of the 50,975 excluded, >45,000 had the corresponding stim cyc < 2014 or age 
# at corresponding stim cyc outside of range; the remaining ~5,600 had the corresponding
# stim cyc >=2014 but belonged to patients whose *first* stim cyc was < 2014 and are
# thus being excluded...)
# 
# sart_frozen %>% slice(19199)
# sart_frozen %>% slice(136794)
# sart_frozen %>% count(num_transferred)
# sart_frozen %>% count(reporting_year)
# sart_frozen %>% count(stim_year1)
# sart_frozen %>% count(stim_year2)
# sart_frozen %>% count(stim_year3)
# sart_frozen %>% count(stim_year4)
# sart_frozen %>% count(stim_year5)
# sart_frozen %>% count(stim_year6)
# sart_frozen %>% count(cycle_order)

#sart_stim %>% count(num_transferred)
#sart_stim %>% count(num_transferred, cycle_type)

# some fresh transfers correspond to id's with transfers before 2014 
# , i.e. initiated first stim before 2014 so need to exclude stims from those women
sart_fresh <- sart_stim %>%
  filter(num_transferred > 0) %>%
  mutate(frozen = 0)

sart_transfers <- bind_rows(sart_fresh, sart_frozen) 


# -------------------- SAVE DATASETS -------------------------------------

# stimulation cycles
saveRDS(sart_stim, file = "Z:/R03/Data/sart_stim.RDS")

# transfer cycles
saveRDS(sart_transfers, file = "Z:/R03/Data/sart_transfers.RDS")


# -----------------------------------------------------------------------------
# ---------------- SPLIT WOMEN into TRAINING and TESTING SETS -----------------
# -----------------------------------------------------------------------------

# since all women in transfers have a stim cycle, but not the other way around,
# use stim cycle dataset to split women into training and testing sets

all_ids <- unique(sart_stim$external_patient_id)

set.seed(20221102)
training_ids <- sample(all_ids, size=ceiling(length(all_ids)/2), replace=FALSE)
testing_ids <- all_ids[which(!(all_ids %in% training_ids))]

set_assign <- data.frame(external_patient_id=training_ids
                         , set="Training Set") %>%
  bind_rows(data.frame(external_patient_id=testing_ids
                       , set="Testing Set"))

saveRDS(set_assign, "Z:/R03/Data/set_assign.RDS")
