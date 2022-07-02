# compare operating characteristics for each model within each stage
# kat correia

library(tidyverse)
library(survival)
library(ranger)     
library(pROC)       

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
# -----------         STAGE 1: DAY 3 OR DAY 5 ET            --------------------
# ------------------------------------------------------------------------------

sart_transfers_set <- sart_transfers %>%
  left_join(set_assign, by="external_patient_id") %>%
  mutate(d5 = ifelse(d5_include==1, yes=1, no=0)) %>%
  filter(!is.na(d5))

dat_dayET_train <- sart_transfers_set %>%
  filter(set=="Training Set") 


# --------------------- LOGISTIC REGRESSION MODELS --------------------------- #

dayET_mod1_slim <- glm(as.formula(str_replace_all(paste0("d5 ~  ", predictors_slim), "\n", ""))
                       , data=dat_dayET_train, family="binomial") 

summary(dayET_mod1_slim)

dayET_mod1_mid <- glm(as.formula(str_replace_all(paste0("d5 ~  ", predictors_mid), "\n", ""))
                      , data=dat_dayET_train, family="binomial") 


summary(dayET_mod1_mid)

dayET_mod1_full <- glm(as.formula(str_replace_all(paste0("d5 ~  ", predictors_full), "\n", ""))
                       , data=dat_dayET_train, family="binomial") 


summary(dayET_mod1_full)


# ---------------------------    RANDOM FOREST     --------------------------- #

dayET_mod2_slim <- ranger(as.formula(str_replace_all(paste0("d5 ~  ", predictors_slim), "\n", ""))
                          , data=dat_dayET_train, classification=TRUE)
dayET_mod2_slim

dayET_mod2_mid <- ranger(as.formula(str_replace_all(paste0("d5 ~  ", predictors_mid), "\n", ""))
                         , data=dat_dayET_train, classification=TRUE)
dayET_mod2_mid

dayET_mod2_full <- ranger(as.formula(str_replace_all(paste0("d5 ~  ", predictors_full), "\n", ""))
                          , data=dat_dayET_train, classification=TRUE)
dayET_mod2_full

# compute predictions at diff thresholds for logistic model
# to get confusion matrix, accuracy, sensitivity, specificity, etc.
dayET_measures <- sart_transfers_set %>%
  select(d5, set, reporting_year, external_cycle_id, external_patient_id
         , agegrp6, clinic_state, amhcat3_new, num_retrieved, bmicat6 
         , gravidity_cat4, parity_cat4, fsh_gt10 
         , male_infertility, dx_tubal, endometriosis, uterine, dx_ovulation 
         , diminished_ovarian_reserve, unexplained) %>%
  mutate(mod1_slim_prob = predict(dayET_mod1_slim, newdata=sart_transfers_set, type = "response")
         , mod1_mid_prob = predict(dayET_mod1_mid, newdata=sart_transfers_set, type = "response")
         , mod1_full_prob = predict(dayET_mod1_full, newdata=sart_transfers_set, type = "response")
         , mod2_slim_pred = predict(dayET_mod2_slim, data=sart_transfers_set)$predictions
         , mod2_mid_pred = predict(dayET_mod2_mid, data=sart_transfers_set)$predictions
         , mod2_full_pred = predict(dayET_mod2_full, data=sart_transfers_set)$predictions
         # compare different thresholds for logistic model
         , mod1_slim_50pred = ifelse(mod1_slim_prob >= .50, yes=1, no=0) 
         , mod1_mid_50pred = ifelse(mod1_mid_prob >= .50, yes=1, no=0) 
         , mod1_full_50pred = ifelse(mod1_full_prob >= .50, yes=1, no=0) 
         , mod1_slim_40pred = ifelse(mod1_slim_prob >= .40, yes=1, no=0) 
         , mod1_mid_40pred = ifelse(mod1_mid_prob >= .40, yes=1, no=0) 
         , mod1_full_40pred = ifelse(mod1_full_prob >= .40, yes=1, no=0) 
         , mod1_slim_30pred = ifelse(mod1_slim_prob >= .30, yes=1, no=0) 
         , mod1_mid_30pred = ifelse(mod1_mid_prob >= .30, yes=1, no=0) 
         , mod1_full_30pred = ifelse(mod1_full_prob >= .30, yes=1, no=0) 
         , mod1_slim_20pred = ifelse(mod1_slim_prob >= .20, yes=1, no=0) 
         , mod1_mid_20pred = ifelse(mod1_mid_prob >= .20, yes=1, no=0) 
         , mod1_full_20pred = ifelse(mod1_full_prob >= .20, yes=1, no=0) 
         , mod1_slim_60pred = ifelse(mod1_slim_prob >= .60, yes=1, no=0) 
         , mod1_mid_60pred = ifelse(mod1_mid_prob >= .60, yes=1, no=0) 
         , mod1_full_60pred = ifelse(mod1_full_prob >= .60, yes=1, no=0) 
         , mod1_slim_70pred = ifelse(mod1_slim_prob >= .70, yes=1, no=0) 
         , mod1_mid_70pred = ifelse(mod1_mid_prob >= .70, yes=1, no=0) 
         , mod1_full_70pred = ifelse(mod1_full_prob >= .70, yes=1, no=0) 
         , mod1_slim_75pred = ifelse(mod1_slim_prob >= .75, yes=1, no=0) 
         , mod1_mid_75pred = ifelse(mod1_mid_prob >= .75, yes=1, no=0) 
         , mod1_full_75pred = ifelse(mod1_full_prob >= .75, yes=1, no=0) 
         , mod1_slim_80pred = ifelse(mod1_slim_prob >= .80, yes=1, no=0) 
         , mod1_mid_80pred = ifelse(mod1_mid_prob >= .80, yes=1, no=0) 
         , mod1_full_80pred = ifelse(mod1_full_prob >= .80, yes=1, no=0)
  )

dayET_measures_test <- dayET_measures %>%
  filter(set=="Testing Set")

auc(dayET_measures_test$d5, predict(dayET_mod1_slim, dayET_measures_test, type="response")) # AUC is 0.720
auc(dayET_measures_test$d5, predict(dayET_mod1_mid, dayET_measures_test, type="response")) # AUC is 0.813
auc(dayET_measures_test$d5, predict(dayET_mod1_full, dayET_measures_test, type="response")) # AUC is 0.816

# create roc curve on test set
roc_slim <- roc(dayET_measures_test$d5, predict(dayET_mod1_slim, dayET_measures_test, type="response"))
plot(roc_slim)
auc(roc_slim)

roc_mid <- roc(dayET_measures_test$d5, predict(dayET_mod1_mid, dayET_measures_test, type="response"))
plot(roc_mid)
auc(roc_mid)

roc_full <- roc(dayET_measures_test$d5, predict(dayET_mod1_full, dayET_measures_test, type="response"))
plot(roc_full)
auc(roc_full)

plot(roc_slim, col="red", print.auc = TRUE, print.auc.y = .4)
plot(roc_mid, add = TRUE, col="green", print.auc = TRUE, print.auc.y = .5)
plot(roc_full, add = TRUE, col="blue", lty="dashed", print.auc = TRUE, print.auc.y = .6)


measures_summary <- dayET_measures %>%
  pivot_longer(cols=c(ends_with("pred")), names_to="model", values_to="pred") %>%
  mutate(temp = case_when(d5==0 & pred==0 ~ "TN"
                          , d5==1 & pred==0 ~ "FN"
                          , d5==0 & pred==1 ~ "FP"
                          , d5==1 & pred==1 ~ "TP")) %>%
  group_by(model, set) %>%
  summarize(TN=sum(temp=="TN")
            , FN=sum(temp=="FN")
            , FP=sum(temp=="FP")
            , TP=sum(temp=="TP")) %>%
  mutate(total=TN+FN+FP+TP
         , accuracy=(TP+TN)/total
         , sens=TP/(TP+FN)
         , spec=TN/(TN+FP)
         , ppv=TP/(TP+FP)
         , npv=TN/(TN+FN)
         , f1=2*(ppv*sens)/(ppv+sens)) %>%
  separate(col=model, into=c("mod", "predictors", "logistic_threshold")
           , remove=FALSE) %>%
  mutate(logistic_threshold=parse_number(logistic_threshold))

ggplot(data=measures_summary, aes(x=logistic_threshold,y=accuracy, color=predictors)) +
  geom_point() +
  geom_hline(data=filter(measures_summary, mod=="mod2"), aes(yintercept=accuracy)
             , lty="dashed", color="black") +
  facet_grid(set ~ predictors) +
  labs(title="Accuracy"
       , x="Logistic regression threshold value") 

ggplot(data=measures_summary, aes(x=logistic_threshold,y=sens, color=predictors)) +
  geom_point() +
  geom_hline(data=filter(measures_summary, mod=="mod2"), aes(yintercept=sens)
             , lty="dashed", color="black") +
  facet_grid(set ~ predictors) +
  labs(title="Sensitivity"
       , x="Logistic regression threshold value") 

ggplot(data=measures_summary, aes(x=logistic_threshold,y=spec, color=predictors)) +
  geom_point() +
  geom_hline(data=filter(measures_summary, mod=="mod2"), aes(yintercept=spec)
             , lty="dashed", color="black") +
  facet_grid(set ~ predictors) +
  labs(title="Specificity"
       , x="Logistic regression threshold value") 

ggplot(data=measures_summary, aes(x=logistic_threshold,y=ppv, color=predictors)) +
  geom_point() +
  geom_hline(data=filter(measures_summary, mod=="mod2"), aes(yintercept=ppv)
             , lty="dashed", color="black") +
  facet_grid(set ~ predictors) +
  labs(title="PPV"
       , x="Logistic regression threshold value") 

ggplot(data=measures_summary, aes(x=logistic_threshold,y=npv, color=predictors)) +
  geom_point() +
  geom_hline(data=filter(measures_summary, mod=="mod2"), aes(yintercept=npv)
             , lty="dashed", color="black") +
  facet_grid(set ~ predictors) +
  labs(title="NPV"
       , x="Logistic regression threshold value")

ggplot(data=measures_summary, aes(x=logistic_threshold,y=f1, color=predictors)) +
  geom_point() +
  geom_hline(data=filter(measures_summary, mod=="mod2"), aes(yintercept=f1)
             , lty="dashed", color="black") +
  facet_grid(set ~ predictors) +
  labs(title="F1"
       , x="Logistic regression threshold value") 

measures_summary %>%
  filter((logistic_threshold %in% c(60,70,75,80) | mod=="mod2") & predictors %in% c("full", "mid")) %>%
  mutate(threshold=ifelse(is.na(logistic_threshold),yes="RF",no=paste0("L (0.",round(logistic_threshold,1),")"))) %>%
  pivot_longer(cols=c("accuracy", "sens", "spec", "ppv", "npv", "f1")
               , names_to="measure", values_to="value") %>%
  ggplot(aes(x=threshold, y=value, color=predictors)) +
  geom_point() +
  facet_grid(set ~ measure) +
  geom_hline(yintercept=0.5,lty="dashed", color="grey")

measures_summary %>%
  ungroup() %>%
  filter((logistic_threshold %in% c(60,70,75,80) | mod=="mod2") & predictors=="mid" & set=="Testing Set") %>%
  mutate(threshold=ifelse(is.na(logistic_threshold),yes="RF",no=paste0("L (0.",round(logistic_threshold,1),")"))) %>%
  select(threshold, accuracy, sens, spec, ppv, npv, f1)

# how are the measures doing over time?
# good, accuracy is *improving* over time . . .
measures_summary_byyear <- dayET_measures %>%
  pivot_longer(cols=c("mod1_mid_60pred", "mod1_mid_70pred", "mod1_mid_75pred", "mod1_mid_80pred", "mod2_mid_pred")
               , names_to="model", values_to="pred") %>%
  mutate(temp = case_when(d5==0 & pred==0 ~ "TN"
                          , d5==1 & pred==0 ~ "FN"
                          , d5==0 & pred==1 ~ "FP"
                          , d5==1 & pred==1 ~ "TP")) %>%
  group_by(model, set, reporting_year) %>%
  summarize(TN=sum(temp=="TN")
            , FN=sum(temp=="FN")
            , FP=sum(temp=="FP")
            , TP=sum(temp=="TP")) %>%
  mutate(total=TN+FN+FP+TP
         , accuracy=(TP+TN)/total
         , sens=TP/(TP+FN)
         , spec=TN/(TN+FP)
         , ppv=TP/(TP+FP)
         , npv=TN/(TN+FN)
         , f1=2*(ppv*sens)/(ppv+sens)
         , x_roc=1-spec) %>%
  separate(col=model, into=c("mod", "predictors", "logistic_threshold")
           , remove=FALSE) %>%
  mutate(logistic_threshold=parse_number(logistic_threshold))

# SEE PREDICTIONS BY AGE GROUP
mosaic::tally(mod1_mid_75pred ~ agegrp6, data=dayET_measures, format="percent")
mosaic::tally(mod1_mid_75pred ~ amhcat3_new, data=dayET_measures, format="percent")
mosaic::tally(mod1_mid_75pred ~ diminished_ovarian_reserve, data=dayET_measures, format="percent")
mosaic::favstats(num_retrieved~mod1_mid_75pred, data=dayET_measures, format="percent")


# ------------------------------------------------------------------------------
# -----------         STAGE 2 - MODEL 1: BLAST RATE         --------------------
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

# ------------------------------ RANDOM FOREST --------------------------------

mod_blast <- rpart(as.formula(str_replace_all(paste0("blast_rate0 ~ ", predictors_mid), "\n", ""))
                   , data = dat_blast_train
                   , method = "anova"   # "class" for a classification tree and "anova" for a regression tree
                   , control=rpart.control(minsplit=2, minbucket=1, cp=0.001))

mod_blast$y <- NULL 
environment(mod_blast$terms) <- NULL
mod_blast$where <- NULL

summary(mod_blast)
predict(mod_blast, newdata=dat_blast_train)

mod_blast2 <- ranger(as.formula(str_replace_all(paste0("blast_rate0 ~ ", predictors_mid), "\n", ""))
                     , data=dat_blast_train)

predict(mod_blast2, data=dat_blast_train)$prediction

# pretty similar predictions
cbind(predict(mod_blast, newdata=dat_blast_train), predict(mod_blast2, data=dat_blast_train)$prediction)

mod_blast2_slim <- ranger(as.formula(str_replace_all(paste0("blast_rate0 ~ ", predictors_slim), "\n", ""))
                          , data=dat_blast_train)

mod_blast2_full <- ranger(as.formula(str_replace_all(paste0("blast_rate0 ~ ", predictors_full), "\n", ""))
                          , data=dat_blast_train)

# see which one minimizes MAE and MSPE on TEST set
pred_d5blast <- dat_blast %>%
  select(set, reporting_year, external_patient_id, external_cycle_id
         , blast_rate, blast_rate0, num_TC, num_TC0
         , agegrp6, clinic_state, amhcat3, diminished_ovarian_reserve, num_retrieved
         , bmicat6, gravidity_cat4, parity_cat4 
         , fsh_gt10, male_infertility, dx_tubal, endometriosis, uterine, dx_ovulation
         , unexplained, transfer_attempted, reason_for_no_transfer) %>%
  mutate(pred_blast_mod1 = predict(mod_blast, newdata=dat_blast)
         , pred_blast_mod2_mid = predict(mod_blast2, data=dat_blast)$prediction
         , pred_blast_mod2_slim = predict(mod_blast2_slim, data=dat_blast)$prediction
         , pred_blast_mod2_full = predict(mod_blast2_full, data=dat_blast)$prediction
         , diff_mod1 = pred_blast_mod1 - blast_rate0
         , diff_mod2_mid = pred_blast_mod2_mid - blast_rate0
         , diff_mod2_slim = pred_blast_mod2_slim - blast_rate0
         , diff_mod2_full = pred_blast_mod2_full - blast_rate0)

pred_d5blast %>%
  group_by(set) %>%
  summarize(n=n()
            , mae1 = sum(abs(diff_mod1))/n
            , mae2_mid = sum(abs(diff_mod2_mid))/n
            , mae2_slim = sum(abs(diff_mod2_slim))/n
            , mae2_full = sum(abs(diff_mod2_full))/n
            , mse1 = sqrt(sum(diff_mod1^2)/n)
            , mse2_mid = sqrt(sum(diff_mod2_mid^2)/n)
            , mse2_slim = sqrt(sum(diff_mod2_slim^2)/n)
            , mse2_full = sqrt(sum(diff_mod2_full^2)/n))


# ------------------------------------------------------------------------------
# -----------         STAGE 2 - MODEL 2: NUMBER ET NEEDED   --------------------
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
  # include the mixed d3/d5 in d5 (good prognosis) group (see reasoning above)
  filter(maxD5==1) 

dat_numETd5_train <- dat_numETd5 %>%
  filter(set=="Training Set")

# ---------------------- fit KM model ----------------------------------------

mod_ETd5_km <- survfit(Surv(totalET, lb) ~  agegrp6 + clinic_state, data = dat_numETd5)

median_km_dat <- summary(mod_ETd5_km)$table[, "median"] %>% 
  as.data.frame() %>%
  rownames_to_column() %>%
  separate(rowname, into=c("blah", "agegrp6", "blah2", "clinic_state")
           , sep="=|,", remove=FALSE) %>%
  mutate(clinic_state=str_trim(clinic_state)) %>%
  select(-rowname, -blah, -blah2)

names(median_km_dat)[3] <- "median_KM"

# ----------------------  AFT models -----------------------------------

mod_ETd5_aft_slim <- survreg(as.formula(paste0("Surv(totalET, lb) ~ ", predictors_slim))
                             , data=dat_numETd5, dist='loglogistic')

mod_ETd5_aft_mid <- survreg(as.formula(paste0("Surv(totalET, lb) ~ ", predictors_mid))
                            , data=dat_numETd5, dist='loglogistic')

mod_ETd5_aft_full <- survreg(as.formula(str_replace_all(paste0("Surv(totalET, lb) ~ ", predictors_full), "\n", ""))
                             , data=dat_numETd5, dist='loglogistic')

concordance(mod_ETd5_aft_slim, mod_ETd5_aft_mid, mod_ETd5_aft_full)
extractAIC(mod_ETd5_aft_slim)
extractAIC(mod_ETd5_aft_mid)
extractAIC(mod_ETd5_aft_full)
anova(mod_ETd5_aft_slim, mod_ETd5_aft_mid)
anova(mod_ETd5_aft_mid, mod_ETd5_aft_full)

# check loglogistic dist'n assumption
fitted_values <- mod_ETd5_aft_mid$linear.predictors
resids <- (log(mod_ETd5_aft_mid$y[, 1]) - fitted_values) / mod_ETd5_aft_mid$scale

resKM <- survfit(Surv(resids, lb) ~ 1, data = dat_numETd5)
plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "No Live Birth Probability")
xx <- seq(min(resids), max(resids), length.out = 35)
yy <- plogis(xx, lower.tail = FALSE)
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")

pred_ETd5 <- dat_numETd5 %>%
  select(set, external_patient_id, agegrp6, clinic_state, amhcat3, diminished_ovarian_reserve, num_retrieved
         , totalET, censored
         , bmicat6, gravidity_cat4, parity_cat4, fsh_gt10, male_infertility, dx_tubal, endometriosis
         , uterine, dx_ovulation, unexplained) %>%
  left_join(median_km_dat, by=c("agegrp6", "clinic_state")) %>%
  mutate(pred_ETd5_mod_slim = predict(mod_ETd5_aft_slim, newdata=dat_numETd5)
         , pred_ETd5_mod_mid = predict(mod_ETd5_aft_mid, newdata=dat_numETd5)
         , pred_ETd5_mod_full = predict(mod_ETd5_aft_full, newdata=dat_numETd5)
         , diff_mod_slim = pred_ETd5_mod_slim - totalET
         , diff_mod_mid = pred_ETd5_mod_mid - totalET
         , diff_mod_full = pred_ETd5_mod_full - totalET
         , diff_mod_km = median_KM - totalET) 

ggplot(data=pred_ETd5) +
  geom_density(aes(x=pred_ETd5_mod_slim), color="red") +
  geom_density(aes(x=pred_ETd5_mod_mid), color="green") +
  geom_density(aes(x=pred_ETd5_mod_full), color="blue") +
  geom_density(aes(x=median_KM), color="purple") 

pred_ETd5 %>%
  group_by(set, censored) %>%
  summarize(n=n()
            , mae_slim = sum(abs(diff_mod_slim))/n
            , mae_mid = sum(abs(diff_mod_mid))/n
            , mae_full = sum(abs(diff_mod_full))/n
            , mae_km = sum(abs(diff_mod_km), na.rm=TRUE)/n
            , rmse_slim = sqrt(sum(diff_mod_slim^2)/n)
            , rmse_mid = sqrt(sum(diff_mod_mid^2)/n)
            , rmse_full = sqrt(sum(diff_mod_full^2)/n)
            , rmse_km = sqrt(sum(diff_mod_km^2, na.rm=TRUE)/n))



