# Assess operating characteristics for each model and for final prediction
# kat correia

# -----------------------------------------------------------------------------
# ------------------- ASSESS OC FOR MODEL 1: NUMBER ET FOR 1 LB ---------------
# -----------------------------------------------------------------------------

pred_ET <- readRDS("Z:/SARTCORS/pred_ET.RDS")

mosaic::tally(censored~set, data=pred_ET)
mosaic::tally(censored~set, data=pred_ET, format="percent")

# ----------------------- compare oc's by age group ---------------------------

pred_ET_age_dat <- pred_ET %>%
  mutate(agegrp = case_when(patient_age_at_start < 32 ~ "<32"
                            , patient_age_at_start >= 32 & patient_age_at_start < 38 ~ "32-37"
                            , patient_age_at_start >= 38 & patient_age_at_start < 40 ~ "38-39"
                            , patient_age_at_start >= 40 & patient_age_at_start < 43 ~ "40-42"
                            , patient_age_at_start >= 43 ~ "43+")) %>%
  bind_rows(pred_ET %>%
              mutate(agegrp="Overall"))


pred_ET_age1 <- pred_ET_age_dat %>%
  pivot_longer(cols=ends_with(c("_dir","_dirC")), names_to="model", values_to="measure") %>%
  group_by(set, censored, agegrp, model) %>%
  count(measure) %>%
  mutate(percent=n/sum(n)) %>%
  arrange(measure) 

pred_ET_age2 <- pred_ET_age_dat %>%
  pivot_longer(cols=ends_with(c("_diff","_diffC")), names_to="model", values_to="value") %>%
  group_by(set, censored, agegrp, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value")

pred_ET_age <- pred_ET_age1 %>%
  rename(value=percent) %>%
  bind_rows(pred_ET_age2) %>%
  mutate(model2=str_replace_all(model,"_dir|_diff","")) %>%
  pivot_wider(id_cols=c(set, censored,agegrp,measure), names_from=model2, values_from=value) %>%
  arrange(set, censored, agegrp)

pred_ET_age %>%
  filter(measure=="RMSE") %>%
  pivot_longer(cols=starts_with("mod1"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="RMSE") +
  facet_grid(set~censored)

pred_ET_age %>%
  filter(measure=="MAE") %>%
  pivot_longer(cols=starts_with("mod1"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="MAE") +
  facet_grid(set~censored)

pred_ET_age %>%
  filter(measure=="Underestimate") %>%
  pivot_longer(cols=starts_with("mod1"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion of Underestimates") +
  facet_grid(set~censored)

pred_ET_age %>%
  filter(measure=="Overestimate") %>%
  pivot_longer(cols=starts_with("mod1"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion of Overestimates") +
  facet_grid(set~censored)

pred_ET_age %>%
  filter(measure=="Exactly equal") %>%
  pivot_longer(cols=starts_with("mod1"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion of Overestimates") +
  facet_grid(set~censored)

pred_ET_age %>%
  filter(measure %in% c("Exactly equal", "Underestimate", "Overestimate")) %>%
  pivot_longer(cols="mod1", names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, fill=measure)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age Group", y="Proportion of Under/Over estimates") +
  facet_grid(~set)

# New facet label names for censor variable
# censor_lbl_rmse <- c("Live Birth Observed", "Censored \n (RMSE is a minimum RMSE)")
# names(censor_lbl_rmse) <- c("0", "1")
censor_lbl_mae <- c("Live Birth Observed", "Censored \n (MAE is a minimum MAE)")
names(censor_lbl_mae) <- c("0", "1")
# censor_lbl <- c("Live Birth Observed", "Censored")
# names(censor_lbl) <- c("0", "1")

pred_ET_age %>%
  filter(measure=="MAE") %>%
  pivot_longer(cols="mod1", names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=set)) +
  geom_point() +
  labs(x = "Age group", y="Mean Absolute Error") +
  facet_grid(~censored) +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  facet_grid(~censored, labeller = labeller(censored = censor_lbl_mae)) 


# ----------------------- compare oc's by diagnosis ---------------------------

# allow cycles to contribute to all diagnoses selected
pred_ET_dx_dat <- bind_rows(
  pred_ET %>% filter(male_infertility=="Y") %>% mutate(infertdx="Male factor")
  , pred_ET %>% filter(endometriosis=="Y") %>% mutate(infertdx="Endometriosis")
  , pred_ET %>% filter(polycystic_ovaries=="Y") %>% mutate(infertdx="PCOS")
  , pred_ET %>% filter(diminished_ovarian_reserve =="Y") %>% mutate(infertdx="DOR")
  , pred_ET %>% filter(dx_tubal=="Y") %>% mutate(infertdx="Tubal factor")
  , pred_ET %>% filter(uterine=="Y") %>% mutate(infertdx="Uterine factor")
  , pred_ET %>% filter(unexplained=="Y") %>% mutate(infertdx="Unexplained")
  , pred_ET %>% mutate(infertdx="Overall")) %>%
  mutate(infertdx=factor(infertdx, levels=c("Male factor", "Tubal factor", "Uterine factor"
                                            , "Endometriosis"
                                            , "PCOS", "DOR", "Unexplained", "Overall")))

pred_ET_dx_dat %>% count(infertdx)

pred_ET_dx <- pred_ET_dx_dat %>%
  pivot_longer(cols=ends_with(c("_diff","_diffC")), names_to="model", values_to="value") %>%
  group_by(set, censored, infertdx, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value") 

pred_ET_dx %>%
  filter(measure=="MAE" & model=="mod1_diff") %>%
  ggplot(aes(x=infertdx, y=value, color=set)) +
  geom_point() +
  labs(x = "Infertility Diagnosis", y="Mean Absolute Error") +
  facet_grid(~censored) +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  facet_grid(~censored, labeller = labeller(censored = censor_lbl_mae)) +
  theme(axis.text.x = element_text(angle = 50)) +
  ylim(c(1.5,3))



# ---------------- compare oc's by number of oocytes retrieved ----------------

pred_ET_er_dat <- pred_ET %>%
  filter(num_retrieved > 0) %>%
  mutate(numER = case_when(1 <= num_retrieved & num_retrieved <= 5 ~ "1-5"
                           , 5 < num_retrieved & num_retrieved <= 10 ~ "6-10"
                           , 10 < num_retrieved & num_retrieved <= 15 ~ "11-15"
                           , 15 < num_retrieved & num_retrieved <= 20 ~ "16-20"
                           , 20 < num_retrieved & num_retrieved <= 30 ~ "21-30"
                           , 30 < num_retrieved & num_retrieved <= 40 ~ "31-40"
                           , 40 < num_retrieved & num_retrieved <= 50 ~ "41-50"
                           , num_retrieved > 50 ~ ">50")) %>%
  bind_rows(pred_ET %>%
              mutate(numER="Overall")) %>%
  mutate(numER = factor(numER, levels=c("1-5", "6-10", "11-15", "16-20", "21-30", "31-40"
                                        , "41-50", ">50", "Overall")))

pred_ET_er_dat %>% count(numER)

pred_ET_er <- pred_ET_er_dat %>%
  pivot_longer(cols=ends_with(c("_diff","_diffC")), names_to="model", values_to="value") %>%
  group_by(set, censored, numER, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value") 

pred_ET_er %>%
  filter(measure=="MAE" & model=="mod1_diff") %>%
  ggplot(aes(x=numER, y=value, color=set)) +
  geom_point() +
  labs(x = "Number of oocytes retrieved", y="Mean Absolute Error") +
  facet_grid(~censored) +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  facet_grid(~censored, labeller = labeller(censored = censor_lbl_mae)) +
  theme(axis.text.x = element_text(angle = 50)) +
  ylim(c(1.5,3))

# -----------------------------------------------------------------------------
# ------------------------ FERTILIZATION RATE ---------------------------------
# -----------------------------------------------------------------------------

pred_fert <- readRDS("Z:/SARTCORS/pred_fert.RDS")

# ----------------------- compare oc's by age group ---------------------------

pred_fert_age_dat <- pred_fert %>%
  mutate(agegrp = case_when(patient_age_at_start < 32 ~ "<32"
                            , patient_age_at_start >= 32 & patient_age_at_start < 38 ~ "32-37"
                            , patient_age_at_start >= 38 & patient_age_at_start < 40 ~ "38-39"
                            , patient_age_at_start >= 40 & patient_age_at_start < 43 ~ "40-42"
                            , patient_age_at_start >= 43 ~ "43+")) %>%
  bind_rows(pred_fert %>%
              mutate(agegrp="Overall"))

pred_fert_age1 <- pred_fert_age_dat %>%
  pivot_longer(cols=ends_with("_dir"), names_to="model", values_to="measure") %>%
  group_by(set, agegrp, model) %>%
  count(measure) %>%
  mutate(percent=n/sum(n)) %>%
  arrange(measure) 

pred_fert_age2 <- pred_fert_age_dat %>%
  pivot_longer(cols=ends_with("_diff"), names_to="model", values_to="value") %>%
  group_by(set, agegrp, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value")

pred_fert_age <- pred_fert_age1 %>%
  rename(value=percent) %>%
  bind_rows(pred_fert_age2) %>%
  mutate(model2=str_replace_all(model,"_dir|_diff","")) %>%
  pivot_wider(id_cols=c(set,agegrp,measure), names_from=model2, values_from=value) %>%
  arrange(agegrp)

pred_fert_age %>%
  filter(measure=="RMSE") %>%
  pivot_longer(cols=starts_with("mod2"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="RMSE") +
  facet_grid(~set)

pred_fert_age %>%
  filter(measure=="MAE") %>%
  pivot_longer(cols=starts_with("mod2"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="MAE") +
  facet_grid(~set)

pred_fert_age %>%
  filter(measure=="MAE") %>%
  pivot_longer(cols="mod2_cart", names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=set)) +
  geom_point() +
  labs(x = "Age group", y="Mean Absolute Error") +
  guides(color=guide_legend(title=NULL)) # remove label from legend

pred_fert_age %>%
  filter(measure=="Underestimate") %>%
  pivot_longer(cols=starts_with("mod2"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion of Underestimates") +
  facet_grid(~set)

pred_fert_age %>%
  filter(measure=="Overestimate") %>%
  pivot_longer(cols=starts_with("mod2"), names_to="model", values_to="value") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion of Overestimates") +
  facet_grid(~set)

pred_fert_age %>%
  filter(measure %in% c("Exactly equal", "Underestimate", "Overestimate")) %>%
  pivot_longer(cols="mod2_cart", names_to="model", values_to="value") %>%
  mutate(model2=str_replace(model, "mod2_", "")) %>%
  ggplot(aes(x=agegrp, y=value, fill=measure)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age group", y="Proportion of Under/Over estimates") +
  facet_grid(~set)


# ----------------------- compare oc's by reason for ART -----------------------

pred_fert_dx_dat <- bind_rows(
  pred_fert %>% filter(male_infertility=="Y") %>% mutate(infertdx="Male factor")
  , pred_fert %>% filter(endometriosis=="Y") %>% mutate(infertdx="Endometriosis")
  , pred_fert %>% filter(polycystic_ovaries=="Y") %>% mutate(infertdx="PCOS")
  , pred_fert %>% filter(diminished_ovarian_reserve =="Y") %>% mutate(infertdx="DOR")
  , pred_fert %>% filter(dx_tubal=="Y") %>% mutate(infertdx="Tubal factor")
  , pred_fert %>% filter(uterine=="Y") %>% mutate(infertdx="Uterine factor")
  , pred_fert %>% filter(unexplained=="Y") %>% mutate(infertdx="Unexplained")
  , pred_fert %>% mutate(infertdx="Overall")) %>%
  mutate(infertdx=factor(infertdx, levels=c("Male factor", "Tubal factor", "Uterine factor"
                                            , "Endometriosis"
                                            , "PCOS", "DOR", "Unexplained", "Overall")))

pred_fert_dx_dat %>% count(infertdx)

pred_fert_dx <- pred_fert_dx_dat %>%
  pivot_longer(cols=ends_with(c("_diff","_diffC")), names_to="model", values_to="value") %>%
  group_by(set, infertdx, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value") 

pred_fert_dx %>%
  filter(measure=="MAE" & model=="mod2_cart_diff") %>%
  ggplot(aes(x=infertdx, y=value, color=set)) +
  geom_point() +
  labs(x = "Infertility Diagnosis", y="Mean Absolute Error") +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  theme(axis.text.x = element_text(angle = 50))


# --------------- compare oc's by number of oocytes retrieved ------------------

pred_fert_er_dat <- pred_fert %>%
  filter(num_retrieved > 0) %>%
  mutate(numER = case_when(1 <= num_retrieved & num_retrieved <= 5 ~ "1-5"
                           , 5 < num_retrieved & num_retrieved <= 10 ~ "6-10"
                           , 10 < num_retrieved & num_retrieved <= 15 ~ "11-15"
                           , 15 < num_retrieved & num_retrieved <= 20 ~ "16-20"
                           , 20 < num_retrieved & num_retrieved <= 30 ~ "21-30"
                           , 30 < num_retrieved & num_retrieved <= 40 ~ "31-40"
                           , 40 < num_retrieved & num_retrieved <= 50 ~ "41-50"
                           , num_retrieved > 50 ~ ">50")) %>%
  bind_rows(pred_fert %>%
              mutate(numER="Overall")) %>%
  mutate(numER = factor(numER, levels=c("1-5", "6-10", "11-15", "16-20", "21-30", "31-40"
                                        , "41-50", ">50", "Overall")))

pred_fert_er_dat %>% count(numER)

pred_fert_er <- pred_fert_er_dat %>%
  pivot_longer(cols=ends_with(c("_diff","_diffC")), names_to="model", values_to="value") %>%
  group_by(set, numER, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value") 

pred_fert_er %>%
  filter(measure=="MAE" & model=="mod2_cart_diff") %>%
  ggplot(aes(x=numER, y=value, color=set)) +
  geom_point() +
  labs(x = "Number of oocytes retrieved", y="Mean Absolute Error") +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  theme(axis.text.x = element_text(angle = 50)) 


# ------------------------------------------------------------------------------
# -------------------       Final Prediction    --------------------------------
# ------------------------------------------------------------------------------

pred_final <- readRDS("Z:/SARTCORS/pred_final.RDS")

# see dist'n of final predictions
pred_final %>%
  select(set, censored, starts_with("pred_")) %>%
  select(-ends_with(c("diff", "_dir"))) %>%
  pivot_longer(cols=-c(set, censored), names_to="prediction", values_to="value") %>%
  filter(prediction=="pred_pois_cart") %>%
  ggplot(aes(x=value)) +
  geom_density() +
  facet_grid(censored ~ set) +
  labs(x="Predicted number of oocytes to inseminate", title="Final prediction")


# compare predictions
pred_final %>%
  select(set, censored, starts_with("pred_")) %>%
  select(-ends_with(c("diff", "_dir"))) %>%
  pivot_longer(cols=-c(set, censored), names_to="prediction", values_to="value") %>%
  ggplot(aes(x=value, color=prediction)) +
  geom_density() +
  facet_grid(censored ~ set) +
  labs(x="Predicted number of oocytes to inseminate")

pred_final %>%
  select(set, starts_with("pred_")) %>%
  select(-ends_with(c("diff", "_dir"))) %>%
  pivot_longer(cols=-set, names_to="prediction", values_to="value") %>%
  ggplot(aes(x=value, color=set)) +
  geom_density() +
  facet_wrap(~prediction) +
  guides(color="none") +
  labs(x="Predicted number of oocytes to inseminate")

# zoom in on reasonable values
pred_final %>%
  select(set, censored, starts_with("pred_")) %>%
  select(-ends_with(c("diff", "_dir"))) %>%
  pivot_longer(cols=-c(set,censored), names_to="prediction", values_to="value") %>%
  filter(value < 25) %>%
  ggplot(aes(x=value, color=prediction)) +
  geom_density() +
  facet_wrap(censored~set) +
  labs(x="Predicted number of oocytes to inseminate")

pred_final %>%
  select(set, censored, starts_with("pred_")) %>%
  select(-ends_with(c("diff", "_dir"))) %>%
  pivot_longer(cols=-c(set,censored), names_to="prediction", values_to="value") %>%
  filter(value < 25) %>%
  ggplot(aes(x=value, color=set)) +
  geom_density() +
  facet_wrap(~prediction) +
  guides(color="none") +
  labs(x="Predicted number of oocytes to inseminate")

mosaic::favstats(num_retrieved ~ set, data=pred_final)

mosaic::favstats(pred_pois_cart ~ set, data=pred_final)
mosaic::favstats(~pred_pois_cart_ETdiff, data=pred_final)
mosaic::favstats(pred_pois_cart_ETdiff ~ pred_pois_cart_ETdiff_dir, data=pred_final)

mosaic::favstats(~pred_pois_cart_Oocdiff, data=pred_final)
mosaic::favstats(pred_pois_cart_Oocdiff ~ pred_pois_cart_Oocdiff_dir, data=pred_final)

mosaic::favstats(~pred_pois_cart_Ooc.7diff, data=pred_final)
mosaic::favstats(pred_pois_cart_Ooc.7diff ~ pred_final_Ooc.7diff_dir, data=pred_final)

# find 90th and 95th percentiles for number needed to transfer before one lb
# extremely uncommon for a woman to need more than 6 embryos to transfer
quantile(pred_final$totalET, p=0.90)
quantile(pred_final$totalET, p=0.95)

pred_final_lbs <- pred_final %>% filter(censored==0)
quantile(pred_final_lbs$totalET, p=0.90)
quantile(pred_final_lbs$totalET, p=0.95)

# total surplus embryos avoided
surplus <- pred_final %>%
  filter(pred_pois_cart_Oocdiff_dir=="Surplus oocytes" & censored != 1)

sum(surplus$pred_pois_cart_Oocdiff)

# from how many diff women? these two should be equal. yes, good.
nrow(surplus)
length(unique(surplus$external_patient_id))

ggplot(data=surplus) +
  geom_density(aes(x=pred_pois_cart_Oocdiff), color="red") +
  geom_density(aes(x=num_retrieved), color="blue")

ggplot(data=surplus) +
  geom_density(aes(x=pred_pois_cart_Oocdiff))

ggplot(data=surplus) +
  geom_density(aes(x=pred_pois_cart_ETdiff))

# total surplus embryos avoided for test set ("new" patients)
surplus_test <- pred_final %>%
  filter(set=="Testing" & pred_pois_cart_Oocdiff_dir=="Surplus oocytes" & censored != 1)

sum(surplus_test$pred_pois_cart_Oocdiff)

# from how many diff women? these two should be equal. yes, good.
nrow(surplus_test)
length(unique(surplus_test$external_patient_id))

ggplot(data=surplus_test) +
  geom_density(aes(x=pred_pois_cart_Oocdiff), color="red") +
  geom_density(aes(x=num_retrieved), color="blue")

ggplot(data=surplus_test) +
  geom_density(aes(x=pred_pois_cart_ETdiff))
