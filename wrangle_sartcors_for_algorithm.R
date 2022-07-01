# kat correia - feb 2022
# prepare sart cors database for analysis (derive vars, wrangle etc)
# and assign women to training and testing sets
# only include women with first stim cycle between 1/1/14 and 12/31/19

library(tidyverse)

path <- "~/cluster-scratch/R03/Data/"

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

sart_raw <- readRDS(paste0(path, "sart_raw.RDS"))

#sart_raw %>% count(donor_oocyte)
#sart_raw %>% count(donor_embryo)
#sart_raw %>% count(gestational_carrier)

# -------------------- wrangle data ----------------------------------------

# memory errors + R keeps crashing; try breaking up the steps and re-arranging

# 800,738 to 783,391
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
         , total2pn, pgd, starts_with("pgd_")
         , laboratory_routine
         , freeze_all
         , starts_with("embryos_to_uterus")
         , transfer_attempted, reason_for_no_transfer
         , elective_single_embryo_transfer
         , starts_with(c("retrieval_day_of_transfer", "fresh_embryos_cryoed"
                         , "oocytes_cryoed", "blastocyst_transfer"
                         , "start_date_freeze_date_diff_autologous_retrieval"))
         , starts_with("treatment_outcome"), starts_with("pregnancy_outcome")
         , pregnancy_loss_abortion
         , starts_with("grade_embryo"), starts_with("stage_embryo")
         , number_live_born
         , preg, lb, cycle_classification)


# sart_a %>%
#   select(-external_patient_id, -external_cycle_id) %>%
#   gtsummary::tbl_summary(missing_text="Unknown"
#                          #, type =  all_continuous() ~ "continuous2"
#                          , statistic = list(all_continuous2() ~ c("{mean} ({sd})"
#                                                                   , "{median} ({p25}, {p75})"
#                                                                   , "{min} - {max}")
#                                             , all_categorical() ~ "{n} ({p}%)"))


# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_a")) { rm(sart_raw) }

# now compute number of oocytes retrieved and 
# number of embryos transferred in two data steps: 
# first, make all the fields numeric
# then, sum them up
# n= 783,343
sart_b <- sart_a %>%
  mutate(across(starts_with(c("embryos_to_uterus","number_retrieved_autologous_retrieval"
                              , "retrieval_day_of_transfer"
                              , "fresh_embryos_cryoed", "oocytes_cryoed")), ~ parse_number(.x))) %>%
  # there are some repeats of same cycle ID for 16 cycles?? some as many as 8 times.
  # keep first one
  group_by(external_patient_id, external_cycle_id) %>%
  mutate(cyctemp=row_number()) %>%
  ungroup() %>%
  filter(cyctemp==1)

# TO SAVE SPACE, remove each dataset once have next
if (exists("sart_b")) { rm(sart_a) }

sart_c1 <- sart_b %>%
  select(external_patient_id, external_cycle_id, starts_with(c("number_retrieved_autologous_retrieval", "embryos_to_uterus", "fresh_embryos_cryoed"
                                                               , "oocytes_cryoed"))) %>%
  pivot_longer(cols=starts_with(c("number_retrieved_autologous_retrieval", "embryos_to_uterus", "fresh_embryos_cryoed"
                                  , "oocytes_cryoed")), names_to = "variable", values_to = "value") %>%
  mutate(var2 = case_when(str_detect(variable, "number_retrieved_autologous_retrieval") ~ "num_retrieved"
                          , str_detect(variable, "embryos_to_uterus") ~ "num_transferred"
                          , str_detect(variable, "fresh_embryos_cryoed") ~ "num_fresh_embryos_cryoed"
                          , str_detect(variable, "oocytes_cryoed") ~ "num_oocytes_cryoed")) %>%
  group_by(external_patient_id, external_cycle_id, var2) %>%
  summarize(sumvalue0 = sum(value, na.rm=TRUE)
            , nmiss = sum(is.na(value))) %>%
  ## there are 6 retrieved variables, 24 ET variables, 6 cryod variables 
  mutate(sumvalue = case_when(var2=="num_retrieved" & nmiss==6 ~ NA_real_
                              , var2=="num_transferred" & nmiss==24 ~ NA_real_
                              , var2=="num_fresh_embryos_cryoed" & nmiss==6 ~ NA_real_
                              , var2=="num_oocytes_cryoed" & nmiss==6 ~ NA_real_
                              , TRUE ~ sumvalue0)) %>%
  pivot_wider(id_cols = c(external_patient_id, external_cycle_id), names_from=var2, values_from = sumvalue) %>%
  ungroup()

sart_c2 <- sart_b %>%
  select(external_patient_id, external_cycle_id, starts_with(c("retrieval_day_of_transfer"))) %>%
  pivot_longer(cols=starts_with(c("retrieval_day_of_transfer")), names_to = "variable", values_to = "value") %>%
  group_by(external_patient_id, external_cycle_id) %>%
  summarize(min_dayET = min(value, na.rm=TRUE)
            , max_dayET = max(value, na.rm=TRUE)) %>%
  ungroup()

sart_for_R03 <- sart_b %>%
  left_join(sart_c1, by=c("external_patient_id", "external_cycle_id")) %>%
  left_join(sart_c2, by=c("external_patient_id", "external_cycle_id")) %>%
  mutate(min_dayET = ifelse(is.finite(min_dayET), yes=min_dayET, no=NA_real_)
         , max_dayET = ifelse(is.finite(max_dayET), yes=max_dayET, no=NA_real_)
         , num_TC = num_transferred + num_fresh_embryos_cryoed
         , agegrp =  factor(case_when(patient_age_at_start < 35 ~ "<35"
                                      , patient_age_at_start >= 35 & patient_age_at_start < 38 ~ "35-37"
                                      , patient_age_at_start >= 38 & patient_age_at_start < 41 ~ "38-40"
                                      , patient_age_at_start >= 41 & patient_age_at_start < 43 ~ "41-42"
                                      , patient_age_at_start >= 43 ~ ">42")
                            , levels=c("<35", "35-37", "38-40", "41-42", ">42"))
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


# -------------------- save cleaned data --------------------------------------

saveRDS(sart_for_R03, file = paste0(path, "sart_for_R03.RDS"))


# ----------------------------------------------------------------------------
# ------------------ create stim dataset + transfer dataset ------------------
# ----------------------------------------------------------------------------

# -------------------- stimulation cycles -------------------------------------

# if a cycle is a thaw cycle, it will have a non-NULL link_source_external_cycle_id
# if fresh, that field will be NULL
# n=527,091 stimulation cycles
# n=495,704 stimulation cycles that started in 2014 or later
sart_stim0 <- sart_for_R03 %>%
  filter(linked_source_external_cycle_i_ds=="NULL" & reporting_year >= 2014)

sart_stim0 %>% count(reporting_year)
sart_stim0 %>% count(cycle_order)
length(unique(sart_stim0$external_patient_id))

# remove stim cycles from women whose first stim was before 2014
sart_stim_first <- sart_stim0 %>%
  filter(cycle_order==1 & 
           patient_age_at_start >= 18 & patient_age_at_start <= 45)

mosaic::favstats(~patient_age_at_start, data=sart_stim_first)

sart_stim1 <- sart_stim0 %>%
  filter(external_patient_id %in% sart_stim_first$external_patient_id)


# -------------------- transfer cycles (fresh + frozen) -----------------------

# ONLY KEEP TRANSFERS THAT RESULTED FROM STIMS IN 2014+

# frozen cycles with > 0 ET & w corresponding stim cycles in 2014+
sart_frozen0 <- sart_for_R03 %>%
  filter(linked_source_external_cycle_i_ds != "NULL") %>%
  mutate(frozen = 1
         # hand-correct one super long weird linked ID that repeats
         , linked_source_external_cycle_i_ds=ifelse(linked_source_external_cycle_i_ds ==
                                                      "76162769, 76162770, 76162771, 76162771, 76162772, 76162774, 76162774, 76162769, 76162769, 76162770, 76162770, 76162771, 76162772, 76162772, 76162774"
                                                    , yes = "76162769, 76162770, 76162771, 76162772, 76162774"
                                                    , no = linked_source_external_cycle_i_ds))


sart_frozen1 <- sart_frozen0 %>%
  # the ones with 5 and 6 end up repeating the linked stim cycles...
  separate(linked_source_external_cycle_i_ds, into=paste0("link",1:6)
           , remove=FALSE, convert=TRUE) %>%
  select(external_cycle_id, linked_source_external_cycle_i_ds, link1, link2, link3, link4, link5, link6) %>%
  pivot_longer(cols=-c(external_cycle_id, linked_source_external_cycle_i_ds)
               , names_to="name", values_to="source_id") %>%
  filter(!is.na(source_id)) %>%
  select(-name, -linked_source_external_cycle_i_ds) %>%
  left_join(sart_frozen0, by=c("external_cycle_id")) %>%
  # filter on women who meet inclusion criteria (first stim cycle 2014 or later)
  filter(external_patient_id %in% sart_stim_first$external_patient_id) %>%
  # also filter on cycles where SOURCE cycle is in sart_for_R03 so no PGD cycles included
  # (some women had multiple stim cycles, some without PGD and some with...)
  filter(source_id %in% sart_for_R03$external_cycle_id) %>%
  group_by(external_cycle_id) %>%
  mutate(checknum=row_number()) %>%
  select(checknum, source_id, external_cycle_id, external_patient_id, frozen, min_dayET, max_dayET, num_transferred, everything()) %>%
  ungroup()

sart_frozen <- sart_frozen1 %>%
  # the info on multiple rows will be the same (since coming from external_cycle_id, NOT coming from source row)
  # but need sart_frozen1 by source row so can get info on day ET for stim cycles below . . . 
  filter(checknum==1) %>%
  # 2,452 frozen cycles with no ET
  filter(num_transferred != 0) %>%
  # for those missing latest amh in frozen transfer cycle, grab latest amh from fresh cycle
  left_join(select(sart_stim1, external_cycle_id, amhcat3_fresh=amhcat3)
            , by=c("source_id"="external_cycle_id"))


# (of the 50,975 excluded, ~48,000 had the corresponding stim cyc < 2014 or age 
# at corresponding stim cyc outside of range; the remaining ~5,600 had the corresponding
# stim cyc >=2014 but belonged to patients whose *first* stim cyc was < 2014 and are
# thus being excluded...)

# 262,605 fresh transfers but some correspond to id's with transfers
# before 2014 / initiated first stim before 2014 so need to exclude stims from those women
sart_fresh <- sart_stim1 %>%
  filter(num_transferred > 0) %>%
  mutate(frozen = 0)

sart_transfers <- bind_rows(sart_fresh, sart_frozen) %>%
  mutate(amhcat3_new = factor(case_when(amhcat3 != "Missing" ~ amhcat3
                                        , amhcat3=="Missing" & !is.na(amhcat3_fresh) ~ amhcat3_fresh
                                        , amhcat3=="Missing" & is.na(amhcat3_fresh) ~ as.factor("Missing"))
                              , labels=c("< 1", "1 - <4", ">= 4", "Missing"))
         , blast_trans=case_when(blastocyst_transfer_autologous_retrieval1=="Y" |
                                   blastocyst_transfer_autologous_retrieval2=="Y" |
                                   blastocyst_transfer_autologous_retrieval3=="Y" |
                                   blastocyst_transfer_autologous_retrieval4=="Y" |
                                   blastocyst_transfer_autologous_retrieval5=="Y" |
                                   blastocyst_transfer_autologous_retrieval6=="Y" 
                                 ~ 1
                                 , TRUE ~ 0)
         , d5_include = case_when(blast_trans==1 ~ 1
                                  , is.na(min_dayET) ~ NA_real_
                                  , min_dayET >= 5 ~ 1
                                  , TRUE ~ 0)
         , num_blasts = ifelse(d5_include==1, yes=num_TC, no=NA_real_)
         , agegrp6 =  factor(case_when(patient_age_at_start < 32 ~ "<32"
                                       , patient_age_at_start >= 32 & patient_age_at_start < 35 ~ "32-34"
                                       , patient_age_at_start >= 35 & patient_age_at_start < 38 ~ "35-37"
                                       , patient_age_at_start >= 38 & patient_age_at_start < 41 ~ "38-40"
                                       , patient_age_at_start >= 41 & patient_age_at_start < 43 ~ "41-42"
                                       , patient_age_at_start >= 43 ~ ">42")
                             , levels=c("<32", "32-34", "35-37", "38-40", "41-42", ">42"))) %>%
  select(reporting_year, frozen, d5_include, min_dayET, blast_trans, everything())


# -------------- back to stimulation cycles - get day frozen for frozen -------

get_day <- sart_frozen %>% 
  select(source_id, min_dayET, max_dayET) %>%
  # same source can contribute to multiple thaw cycles
  group_by(source_id) %>%
  summarize(min_dayfreeze0=min(min_dayET)
            , max_dayfreeze0=max(max_dayET))

sart_stim <- sart_stim1 %>%
  left_join(get_day, by=c("external_cycle_id"="source_id")) %>%
  mutate(day_freeze = case_when(!is.na(min_dayET) ~ min_dayET
                                , TRUE ~ min_dayfreeze0)
         , blast_trans=case_when(blastocyst_transfer_autologous_retrieval1=="Y" |
                                   blastocyst_transfer_autologous_retrieval2=="Y" |
                                   blastocyst_transfer_autologous_retrieval3=="Y" |
                                   blastocyst_transfer_autologous_retrieval4=="Y" |
                                   blastocyst_transfer_autologous_retrieval5=="Y" |
                                   blastocyst_transfer_autologous_retrieval6=="Y" 
                                 ~ 1
                                 , TRUE ~ 0)
         , d5_include = case_when(blast_trans==1 ~ 1
                                  , is.na(day_freeze) ~ NA_real_
                                  , day_freeze >= 5 ~ 1
                                  , TRUE ~ 0)
         , any_embryos_frozen = ifelse(num_fresh_embryos_cryoed >= 1, yes=1, no=0)
         , fert_rate = ifelse(num_retrieved != 0, yes = total2pn/num_retrieved, no=NA_real_)
         , num_blasts = ifelse(d5_include==1, yes=num_TC, no=NA_real_)
         , blast_rate = num_blasts/num_retrieved
         , fert_to_blast = num_blasts/total2pn
         , agegrp6 =  factor(case_when(patient_age_at_start < 32 ~ "<32"
                                       , patient_age_at_start >= 32 & patient_age_at_start < 35 ~ "32-34"
                                       , patient_age_at_start >= 35 & patient_age_at_start < 38 ~ "35-37"
                                       , patient_age_at_start >= 38 & patient_age_at_start < 41 ~ "38-40"
                                       , patient_age_at_start >= 41 & patient_age_at_start < 43 ~ "41-42"
                                       , patient_age_at_start >= 43 ~ ">42")
                             , levels=c("<32", "32-34", "35-37", "38-40", "41-42", ">42"))) %>%
  select(reporting_year, external_patient_id, external_cycle_id, day_freeze
         , min_dayET, max_dayET, min_dayfreeze0, max_dayfreeze0, num_transferred
         , num_retrieved, num_TC, total2pn, num_blasts, fert_rate, blast_rate, fert_to_blast
         , everything()) 


# -------------------- SAVE DATASETS -------------------------------------

saveRDS(sart_stim, file = paste0(path, "sart_stim.RDS"))

saveRDS(sart_transfers, file = paste0(path, "sart_transfers.RDS"))


# -----------------------------------------------------------------------------
# ---------------- SPLIT women up into training and testing sets        -------
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

saveRDS(set_assign, paste0(path, "set_assign.RDS"))
