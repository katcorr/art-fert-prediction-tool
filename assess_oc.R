# compare operating characteristics for each model within each step
# kat correia

library(tidyverse)

# -----------------------------------------------------------------------------
# ------------------- COMPARE OC FOR MODEL 1: NUMBER ET FOR 1 LB ---------------
# -----------------------------------------------------------------------------

pred_ET_byage <- readRDS("Z:/R03/Data/pred_ET_byage.RDS") %>%
  mutate(agegrp=factor(agegrp, levels=c("<35", "35-37", "38-40", "41-42", ">42")))

mosaic::tally(censored~set, data=pred_ET_byage)
mosaic::tally(censored~set, data=pred_ET_byage, format="percent")

var_pred_vec <- c("mod1_cp_pred", "mod1_cp_predR", "mod1_cp_predC"
                  , "mod1_lm_pred", "mod1_lm_predR", "mod1_lm_predC"
                  , "mod1_q50_pred", "mod1_q50_predR", "mod1_q50_predC"
                  , "perc50", "perc75", "perc90")

# get truncated prediction error curve for total ET between 5 and 30
pe_curve0 <- pred_ET_byage %>%
  rowwise() %>%
  mutate(across(all_of(var_pred_vec), ~ min(totalET, 5) - min(., 5), .names="l5_{.col}")
        , across(all_of(var_pred_vec), ~ min(totalET, 10) - min(., 10), .names="l10_{.col}")
        , across(all_of(var_pred_vec), ~ min(totalET, 20) - min(., 20), .names="l20_{.col}")
        , across(all_of(var_pred_vec), ~ min(totalET, 30) - min(., 30), .names="l30_{.col}")) %>%
  ungroup()

pe_curve0 %>%
  slice(1:10) %>%
  select(weight, censored, totalET, perc50, l30_perc50, mod1_q50_pred, l30_mod1_q50_pred)

pe_curve <- pe_curve0 %>%
  select(set, agegrp, weight, censored, totalET, all_of(var_pred_vec), starts_with("l30_")) %>%
  pivot_longer(cols=c("l30_perc50", "l30_mod1_q50_pred") #, "l30_mod1_lm_pred", "l30_mod1_cp_pred)
               , names_to="model", values_to="value") %>%
  mutate(weighted_abs_value = abs(value)*weight) %>%
  #select(model, value, weight, weighted_value, everything())
  group_by(set, agegrp, model) %>%
  summarize(n=n()
            , PE = sum(abs(value)/n)
            , PE_IPW = sum(weighted_abs_value)/n)

ggplot(data=pe_curve, aes(x=agegrp, y=PE_IPW, color=model, group=model))+ 
  geom_point() +
  geom_line() +
  facet_wrap(~set)

ggplot(data=pe_curve, aes(x=agegrp, y=PE, color=model, group=model))+ 
  geom_point() +
  geom_line() +
  facet_wrap(~set)

# can easily change grouping variable to look at performance among specific groups
pred_ET_age_dat <- pred_ET_byage %>%
  bind_rows(pred_ET_byage %>%
              mutate(agegrp="Overall")) %>%
  mutate(agegrp=factor(agegrp, levels=c("<35", "35-37", "38-40", "41-42", ">42")))

pred_ET_age1 <- pred_ET_age_dat %>%
  pivot_longer(cols=ends_with(c("_dir","_dirC")), names_to="model", values_to="measure") %>%
  group_by(set, censored, agegrp, model) %>%
  count(measure) %>%
  mutate(percent=n/sum(n)) %>%
  arrange(measure) %>%
  ungroup()

pred_ET_age2 <- pred_ET_age_dat %>%
  pivot_longer(cols=c("perc50_diff", "mod1_q50_pred_diff", "mod1_lm_pred_diff")# ,"mod1_cp_pred_diff")
               , names_to="model", values_to="value") %>%
  group_by(set, censored, agegrp, model) %>%
  summarize(n=n()
            , RMSE = sqrt(sum(value^2)/n)
            , MAE = sum(abs(value))/n) %>%
  pivot_longer(cols=c("RMSE", "MAE"), names_to="measure", values_to="value") %>%
  ungroup()

pred_ET_age2 %>%
  filter(measure=="MAE") %>%
  ggplot(aes(x=agegrp, y=value, color=model)) +
  geom_point() +
  labs(x = "Age group", y="MAE") +
  facet_grid(set~censored)

pred_ET_age1 %>%
  filter(measure=="Underestimate" & !str_detect(model, "C|R")) %>%
  ggplot(aes(x=agegrp, y=percent, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion of Underestimates") +
  facet_grid(set~censored)

pred_ET_age1 %>%
  filter(measure=="Exactly equal" & str_detect(model, "R")) %>%
  ggplot(aes(x=agegrp, y=percent, color=model)) +
  geom_point() +
  labs(x = "Age group", y="Proportion Exactly Equal") +
  facet_grid(set~censored)

pred_ET_age1 %>%
  filter(model=="mod1_lm_predR_diff_dir") %>%
  filter(measure %in% c("Exactly equal", "Underestimate", "Overestimate")) %>%
  ggplot(aes(x=agegrp, y=percent, fill=measure)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age Group", y="Proportion of Under/Over estimates", title="LM Rounded") +
  facet_grid(set ~ censored)

pred_ET_age1 %>%
  filter(model=="mod1_lm_predC_diff_dir") %>%
  filter(measure %in% c("Exactly equal", "Underestimate", "Overestimate")) %>%
  ggplot(aes(x=agegrp, y=percent, fill=measure)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age Group", y="Proportion of Under/Over estimates", title="LM Ceiling") +
  facet_grid(set ~ censored)

pred_ET_age1 %>%
  filter(model=="mod1_q50_pred_diff_dir") %>%
  filter(measure %in% c("Exactly equal", "Underestimate", "Overestimate")) %>%
  ggplot(aes(x=agegrp, y=percent, fill=measure)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age Group", y="Proportion of Under/Over estimates"
       ,  title="Median Regression") +
  facet_grid(set ~ censored)

pred_ET_byage %>%
  pivot_longer(cols=c("mod1_cp_predR", "mod1_lm_predR", "mod1_q50_pred", "perc50")
               , names_to="model", values_to="value") %>%
  ggplot(aes(x=value, fill=model)) +
  geom_histogram() +
  facet_wrap(~model)

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

pred_ET_er_dat <- pred_ET_byage %>%
  filter(num_retrieved > 0) %>%
  mutate(numER = case_when(1 <= num_retrieved & num_retrieved <= 5 ~ "1-5"
                           , 5 < num_retrieved & num_retrieved <= 10 ~ "6-10"
                           , 10 < num_retrieved & num_retrieved <= 15 ~ "11-15"
                           , 15 < num_retrieved & num_retrieved <= 20 ~ "16-20"
                           , 20 < num_retrieved & num_retrieved <= 30 ~ "21-30"
                           , 30 < num_retrieved & num_retrieved <= 40 ~ "31-40"
                           , 40 < num_retrieved & num_retrieved <= 50 ~ "41-50"
                           , num_retrieved > 50 ~ ">50")) %>%
  bind_rows(pred_ET_byage %>%
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
  filter(measure=="MAE") %>%
  ggplot(aes(x=numER, y=value, color=model)) +
  geom_point() +
  labs(x = "Number of oocytes retrieved", y="Mean Absolute Error") +
  facet_grid(set~censored) +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  theme(axis.text.x = element_text(angle = 50)) 

# -----------------------------------------------------------------------------
# ------------------- COMPARE OC FOR MODEL 2: FERTILIZATION RATE --------------
# -----------------------------------------------------------------------------

pred_fert_byage <- readRDS("Z:/R03/Data/pred_fert_byage.RDS") %>%
  select(set, agegrp, num_retrieved, male_infertility, endometriosis, polycystic_ovaries
         , diminished_ovarian_reserve, dx_tubal, uterine, unexplained
         , starts_with(c("mod2", "perc")))

# ----------------------- compare oc's by age group ---------------------------

pred_fert_age_dat <- pred_fert_byage %>%
  bind_rows(pred_fert_byage %>%
              mutate(agegrp="Overall")) %>%
  mutate(agegrp=factor(agegrp, levels=c("<35", "35-37", "38-40", "41-42", ">42")))

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
  pivot_longer(cols="mod2_cart_pred", names_to="model", values_to="value") %>%
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
  pivot_longer(cols="mod2_cart_pred", names_to="model", values_to="value") %>%
  mutate(model2=str_replace(model, "mod2_", "")) %>%
  ggplot(aes(x=agegrp, y=value, fill=measure)) +
  geom_bar(position="fill", stat="identity") +
  labs(x = "Age group", y="Proportion of Under/Over estimates") +
  facet_grid(~set)


# ----------------------- compare oc's by reason for ART -----------------------

pred_fert_dx_dat <- bind_rows(
  pred_fert_byage %>% filter(male_infertility=="Y") %>% mutate(infertdx="Male factor")
  , pred_fert_byage %>% filter(endometriosis=="Y") %>% mutate(infertdx="Endometriosis")
  , pred_fert_byage %>% filter(polycystic_ovaries=="Y") %>% mutate(infertdx="PCOS")
  , pred_fert_byage %>% filter(diminished_ovarian_reserve =="Y") %>% mutate(infertdx="DOR")
  , pred_fert_byage %>% filter(dx_tubal=="Y") %>% mutate(infertdx="Tubal factor")
  , pred_fert_byage %>% filter(uterine=="Y") %>% mutate(infertdx="Uterine factor")
  , pred_fert_byage %>% filter(unexplained=="Y") %>% mutate(infertdx="Unexplained")
  , pred_fert_byage %>% mutate(infertdx="Overall")) %>%
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
  filter(measure=="MAE" & model=="mod2_cart_pred_diff") %>%
  ggplot(aes(x=infertdx, y=value, color=set)) +
  geom_point() +
  labs(x = "Infertility Diagnosis", y="Mean Absolute Error") +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  theme(axis.text.x = element_text(angle = 50))


# --------------- compare oc's by number of oocytes retrieved ------------------

pred_fert_er_dat <- pred_fert_byage %>%
  filter(num_retrieved > 0) %>%
  mutate(numER = case_when(1 <= num_retrieved & num_retrieved <= 5 ~ "1-5"
                           , 5 < num_retrieved & num_retrieved <= 10 ~ "6-10"
                           , 10 < num_retrieved & num_retrieved <= 15 ~ "11-15"
                           , 15 < num_retrieved & num_retrieved <= 20 ~ "16-20"
                           , 20 < num_retrieved & num_retrieved <= 30 ~ "21-30"
                           , 30 < num_retrieved & num_retrieved <= 40 ~ "31-40"
                           , 40 < num_retrieved & num_retrieved <= 50 ~ "41-50"
                           , num_retrieved > 50 ~ ">50")) %>%
  bind_rows(pred_fert_byage %>%
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
  filter(measure=="MAE" & model=="mod2_cart_pred_diff") %>%
  ggplot(aes(x=numER, y=value, color=set)) +
  geom_point() +
  labs(x = "Number of oocytes retrieved", y="Mean Absolute Error") +
  guides(color=guide_legend(title=NULL)) + # remove label from legend
  theme(axis.text.x = element_text(angle = 50)) 
