# rm(list = ls())
library(data.table)
library(tidyverse)
library(survey)
library(binom)
library( magrittr)
library(feather)

setwd("/home/laurette/Desktop/Github/MSF_data")

source('~/Desktop/Github/utility/SourceLibs.R', echo = FALSE)
source('~/Desktop/Github/SurveySimulations/estimated_prevalences2.R')
source('~/Desktop/Github/SurveySimulations/grebe_fittedprevalencemodels.R')

# msfdata2013old <- read.csv("/home/laurette/Downloads/pone.0203638s003.csv")

msfdata2018 <- read.csv("/home/laurette/Desktop/Github/MSF_data/07032020_hivdata_toshare_sacema.csv")


# msfdata2018 <- read.csv("/home/laurette/Desktop/Github/MSF_data/07032020_hivdata_toshare_sacema.csv") %>% 
#   select(unique_id, cluster_id, gender, finalage, hivtest_results_final,
#          Lag_Classification_Final, VL_result_new, Odn_Final,vl1000_90, Final_Classification) %>% 
#   mutate(recent = ifelse(VL_result_new > 75 & Odn_Final< 2,1,0 ))

mortalitydata <- read_feather("/home/laurette/Desktop/Github/MSF_data/incidence-combined-method-master/mortality_table.feather") 

# hist(msfdata2018$)

# goes from 2010 upto 2015
# convert data online from Grebe et al to a data.frame 

msfdata2013 <- read_feather("/home/laurette/Desktop/Github/incidence-combined-method-master/kzn.feather") %>% 
  select(ID, cluster, ward, sex, age_years, hiv, naat_final,
         LAgFinal, biorad_fin, viral_load_fin, recent)

msfdata2013$intdate <- "2013-06-06"


############################################################################

# msfdata2018f <- msfdata2018 %>%
#   summarise(total = sum(ifelse(gender == 2, 1, 0)),
#             specimens = n(),
#             difference = specimens - total )
  
msfdata2018_mod <- msfdata2018 %>%
  transmute(household_id = household_id,
            id = as.factor(X), cluster = cluster_id,
            ward = NA, sex = ifelse(gender == 2, "Female", "Male"),
            age_years = as.numeric(as.character(finalage)),
            hiv_status = ifelse(hivtest_results_final == 1, 1,
                                          ifelse(hivtest_results_final == 2,  0, NA)),
            LAgODn = Odn_Final,
            viral_load = VL_result_new,
            intdate = as.character(intdate),
            # recent = ifelse(Final_Classification == "RECENT", 1,
            #                 ifelse(Final_Classification == "LT", 0, NA)),
            recent = ifelse(viral_load > 75 & LAgODn < 2, 1,0 ),
            dates.x = 2018)


# hist(msfdata2018_mod$LAgODn[msfdata2018_mod$LAgODn<=3], stats.bins = 45)

msfdata2013_mod <- msfdata2013 %>%
  transmute(id = ID, cluster = cluster,
            ward = ward, sex = sex,
            age_years = as.integer( age_years),
            hiv_status = hiv,
            naat_final = naat_final,
            LAgODn = LAgFinal,
            biorad_fin = biorad_fin,
            viral_load = viral_load_fin,
            intdate = intdate,
            recent = ifelse(viral_load > 75 & LAgODn < 2, 1,0 ), #recent, # for combining the data from 2013 to 2018 the recency case 
            # defination is redefined leaving out the NAAT and BioRAD.
            dates.x = 2013) 


msfdata2013_original <- msfdata2013 %>%
  transmute(id = ID, cluster = cluster,
            ward = ward, sex = sex,
            age_years = as.integer( age_years),
            hiv_status = hiv,
            naat_final = naat_final,
            LAgODn = LAgFinal,
            biorad_fin = biorad_fin,
            viral_load = viral_load_fin,
            intdate = intdate,
            recent = recent, # for combing the data the recency case 
            # defination is redefined to NAAT and BioRAD.
            dates.x = 2013)


# combined data 
 msfdata <- bind_rows(msfdata2018_mod, msfdata2013_mod) %>% 
   filter(age_years<=59)
 write_feather(msfdata2013_mod, "kzn.feather")

############################################################################################################################################
############################################################################################################################################

# Summarise and structure msf survey data into required data structure!!

############################################################################################################################################
############################################################################################################################################

msfdata_summary <- msfdata %>% 
   # filter(dates.x == 2013) %>% 
  group_by(dates.x, age_years) %>% 
  summarise(n_specimens = n(),
            positive = sum(hiv_status, na.rm = T), 
            negative  = n_specimens - positive,
            recent  = sum(ifelse(recent  == 1, 1, 0), na.rm = T),
            nonrecent = positive - recent) %>% 
  mutate(age = age_years, overallweight = 1,
         prevalence = positive/ (positive + negative))

msfdata_summary$age = as.numeric(as.character(msfdata_summary$age_years))



dim(msfdata_summary)


populationest <- msfdata_summary %>% 
  group_by(age_years) %>%
  # filter(age_years >= 15 & age_years <= 35) %>% 
  summarise(total = sum( n_specimens)) %>% 
  mutate(props = total/sum(total)) 
  

names(populationest) <- c("age", "total", "props")



# msfdata %>% 
#   ggplot()
  


msfdata_summary_date <- msfdata %>% 
  group_by(dates.x) %>% 
  summarise(n_specimens = n()) 


newdata <- bind_rows(msfdata2018_mod, msfdata2013_mod) %>%
  filter(!is.na(hiv_status))

  ggplot(newdata, aes(x = age_years, fill = as.factor(hiv_status)))+
    geom_histogram()+
    facet_grid(rows = vars(dates.x), cols = vars(sex))+
    scale_x_continuous(limits=c(15,105), breaks = seq(20,100, 10))+
    theme_bw(base_size = 32, base_family = "") +
    labs(x = "age", y = "frequency ", fill = "")+
    scale_fill_manual(labels = c( "negative", "positive"),
                      values = c( "turquoise3", "red"))+ 
    theme(legend.position = "bottom")
  
  ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/histogram.png", 
         w=15,h= 10)

  
  
  
  
  
 msfdata <- bind_rows(msfdata2013_mod) %>% 
   group_by(age_years) %>% 
    summarise(Hivptve = sum(hiv_status, na.rm = T ),
              Recent = sum(recent, na.rm = T ))
 
 # %>% 
    mutate(props = total/sum(total)) 

############################################################################################################################################
# ############################################################################################################################################
# 
# # Analysis parameters
# 
# ############################################################################################################################################
# ############################################################################################################################################
# 
# predicttimes <- sum(unique(msfdata_summary$dates.x), na.rm = T)/2
# anchorage = age_predict = unique(msfdata_summary$age)[unique(msfdata_summary$age) <= 50]
# 
# MDRI = 207 # comment on 180 
#   
# FRR = 0.0018/100
# RSE_MDRI = 0.06909014; RSE_FRR = 0.25; epsi <- (MDRI/365.25)/2; BigT = 365.25
# 
# ############################################################################################################################################
# ############################################################################################################################################
# 
# #Females 
# 
# msfdata_summary2018 <- msfdata %>% 
#    filter(dates.x == 2013 ) %>% 
#   group_by(dates.x) %>%
#   # filter( age_years <= 30) %>% 
#   summarise(n_specimens = n(),
#             positive = sum(hiv_status, na.rm = T), 
#             negative  = n_specimens - positive,
#             recent  = sum( recent, na.rm = T),
#             nonrecent = positive - recent) 
# 
# 
# prev = binom.exact(x = msfdata_summary2018$positive, msfdata_summary2018$positive +msfdata_summary2018$negative)
# prev$sd <- (prev$lower - prev$upper)/(2 *1.96)
# 
# 
# rec = binom.exact(x = msfdata_summary2018$recent, msfdata_summary2018$positive )
# rec$sd <- (rec$lower - rec$upper)/(2 * 1.96)
# 
# incidenceK_females <- tryCatch(incprops(PrevH = prev$mean, PrevR = rec$mean,
#                                 RSE_PrevH = prev$sd/ prev$mean,
#                                 RSE_PrevR = rec$sd / rec$mean,
#                                 MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR =  FRR,
#                                 RSE_FRR = RSE_FRR, BigT = BigT), error = function(e) NULL)
# 
# incidenceK_females_est <- as.numeric(levels(incidenceK_females$Incidence.Statistics$Incidence))[incidenceK_females$Incidence.Statistics$Incidence]
# incidenceK_females_se <- as.numeric(levels(incidenceK_females$Incidence.Statistics$RSE))[incidenceK_females$Incidence.Statistics$RSE]* incidenceK_females_est
# 
# 
# ################################################################################################################
# 
# males_msfdata_summary2018 <- msfdata %>% 
#   filter(dates.x == 2018 & sex == "Male") %>% 
#   group_by(dates.x) %>%
#   filter( age_years <= 30) %>% 
#   summarise(n_specimens = n(),
#             positive = sum(hiv_status, na.rm = T), 
#             negative  = n_specimens - positive,
#             recent  = sum( recent, na.rm = T),
#             nonrecent = positive - recent) 
# 
# 
# prev_males = binom.exact(x = males_msfdata_summary2018$positive, males_msfdata_summary2018$positive + males_msfdata_summary2018$negative)
# prev_males$sd <- (prev_males$lower - prev_males$upper)/(2*1.96)
# 
# 
# rec_males = binom.exact(x = males_msfdata_summary2018$recent, males_msfdata_summary2018$positive )
# rec_males$sd <- (rec_males$lower - rec_males$upper)/(2*1.96)
# 
# incidenceK_males <- tryCatch(incprops(PrevH = prev_males$mean, PrevR = rec_males$mean,
#                                 RSE_PrevH = prev_males$sd/ prev_males$mean,
#                                 RSE_PrevR = rec_males$sd / rec_males$mean,
#                                 MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR =  FRR,
#                                 RSE_FRR = RSE_FRR, BigT = BigT), error = function(e) NULL)
# 
# incidenceK_males_est <- as.numeric(levels(incidenceK_males$Incidence.Statistics$Incidence))[incidenceK_males$Incidence.Statistics$Incidence]
# incidenceK_males_se <- as.numeric(levels(incidenceK_males$Incidence.Statistics$RSE))[incidenceK_males$Incidence.Statistics$RSE]* incidenceK_males_est
# 
# 
# ################################################################################################################
# 
# combined_msfdata_summary2018 <- msfdata %>% 
#   filter(dates.x == 2018 ) %>% 
#   group_by(dates.x) %>%
#   filter( age_years <= 30) %>% 
#   summarise(n_specimens = n(),
#             positive = sum(hiv_status, na.rm = T), 
#             negative  = n_specimens - positive,
#             recent  = sum( recent, na.rm = T),
#             nonrecent = positive - recent) 
# 
# 
# prev_combined = binom.exact(x = combined_msfdata_summary2018$positive, combined_msfdata_summary2018$positive + combined_msfdata_summary2018$negative)
# prev_combined$sd <- (prev_combined$lower - prev_combined$upper)/(2*1.96)
# 
# 
# rec_combined = binom.exact(x = combined_msfdata_summary2018$recent, combined_msfdata_summary2018$positive )
# rec_combined$sd <- (rec_combined$lower - rec_combined$upper)/(2*1.96)
# 
# Combined_incidenceK <- tryCatch(incprops(PrevH = prev_combined$mean, PrevR = rec_combined$mean,
#                                 RSE_PrevH = prev_combined$sd/ prev_combined$mean,
#                                 RSE_PrevR = rec_combined$sd / rec_combined$mean,
#                                 MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR =  FRR,
#                                 RSE_FRR = RSE_FRR, BigT = BigT), error = function(e) NULL)
# 
# 
# Combined_incidenceKest <- as.numeric(levels(Combined_incidenceK$Incidence.Statistics$Incidence))[Combined_incidenceK$Incidence.Statistics$Incidence]
# Combined_incidenceKse <- as.numeric(levels(Combined_incidenceK$Incidence.Statistics$RSE))[Combined_incidenceK$Incidence.Statistics$RSE]* Combined_incidenceKest
# 
# 
# ######################################################################################################
# incidenceestimates <- data.frame(pop_group = c("females", "males", "combined"),
#                                  incidence = c(incidenceK_females_est, incidenceK_males_est, Combined_incidenceKest),
#                                  std_error = c( incidenceK_males_se, incidenceK_males_se,Combined_incidenceKse) )
