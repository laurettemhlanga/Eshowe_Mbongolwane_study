library(tidyverse)
library(survey)
library(binom)
library( magrittr)
library(feather)
library(foreach)
library(doParallel)

setwd("/home/laurette/Desktop/Github/MSF_data")

source('~/Desktop/Github/utility/SourceLibs.R', echo = FALSE)
source('~/Desktop/Github/SurveySimulations/estimated_prevalences2.R')
source('~/Desktop/Github/SurveySimulations/grebe_fittedprevalencemodels.R')
source('~/Desktop/Github/inctools/inctools/R/incidence.R', echo=F)
source('~/Desktop/Github/MSF_data/MSF_datacleaning.R')
source('~/Desktop/Github/inctools/inctools/R/incidence.R', echo=F)
source('~/Desktop/Github/SurveySimulations/log_fittedprevalencemodels.R')
source('~/Desktop/Github/msfanalysis20211110/bootstrap_surveydata.R')
source('~/Desktop/Github/msfanalysis20211110/mortality_sensitivity.R')

################################################################################################################
# uses the data provided by Eduard - kznfeather
# REQUIRED uses the new sampling strategy to estimate the 
# mid point and survey times incidence estimates.
# fits one model on all the data from 15:45 
# the link function uses the logit for prevalence and 
# clog log for the prevalence of recency function.
# anchorage is set at 25 to ensure use of all data points between 15:35
################################################################################################################

mortalitydata <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/mortality2018") 

excess_mortality <-  mortalitydata %>% 
  filter(Sex == "Male") %>%
  transmute(age = Age, sex = Sex,
            X2013 = X2013,
            X2018 = X2018) %>%
  melt(id.vars = c("age","sex")) %>%
  filter(age <= 40)

names(excess_mortality) <- c("age", "sex", "dates", "mortality")

excess_mortality$dates = ifelse(excess_mortality$dates == "X2013", 2013,
                                ifelse(excess_mortality$dates == "X2015", 2015.5, 2018))

excessmortality <-  excess_mortality %>% 
  transmute(age = age,
            dates = as.numeric(as.character(dates)),
            mortality = mortality)

######################################
# input parameters
######################################

timecutoff = 6; agetimetrans = F; points_weights = 1; 
inclusionmethod = "distance"; overall_weights = FALSE
exponentterm = F; recency = T; prevlink = "logit"; 
prevtaylororder = 4; epsilon = 0.01; rectaylororder = 1
previnclusion_dist = 20; recinclusion_dist = 0;  
anchorage <- 20; age_predict <- 20:30; oneparam = F; 
agetimetrans = F; len_anchorage <- length(anchorage)
midpoint <- 2015.5; predicttimes <- c(2013, 2018)

######################################
mortality_factor <- seq(0.5,2, 0.25)
######################################
# spliting the survey data
######################################
surveyiterations <- 10000

clustersin2013 <- unique(msfdata[msfdata$dates.x == 2013,]$cluster);
clustersin2018 <- unique(msfdata[msfdata$dates.x == 2018,]$cluster)

ncluster13 = length(clustersin2013)
ncluster18 = length(clustersin2018)
nclusters = ncluster18 + ncluster13
nward = 14

bootstrap_strategy = "both"

msfdata2013mod <- msfdata2013_original %>% 
  filter(sex == "Male")


msfdata2018mod <- msfdata2018_mod %>% 
  filter(sex == "Male")

msfdata <- bind_rows(msfdata2013mod, msfdata2018mod)

# preprocess 2013
msfdata2013_split <- split(msfdata2013mod, list(msfdata2013mod$ward))

clusters_in_ward <- lapply(seq_along(msfdata2013_split), 
                           function(x) unique(msfdata2013_split[[x]]$cluster))

#preprocess 2018
msfdata2018_split <- split(msfdata2018mod, list(msfdata2018mod$cluster))


registerDoParallel(detectCores())
getDoParWorkers()
######################################
#  Bootstraps 
######################################


start_time <- Sys.time()

bootstrap_incidence <- foreach(seqindex  =  seq_along(mortality_factor)) %:%
  foreach(index = seq_along(1:surveyiterations), 
          .combine = "rbind") %dopar% {
                                 
             excessmortality <-  excess_mortality %>% 
               transmute(age = age,
                         dates = as.numeric(as.character(dates)),
                         mortality = mortality * mortality_factor[seqindex])
             
             analysis_data <- bootstrap_data(data_2013  = msfdata2013_split, 
                                             data_2018 = msfdata2018_split,
                                             bootstrap_strategy = bootstrap_strategy)
             
             incidence <- sensitivity_mortality(msfdata_summary  = analysis_data)
             
             incidence$mortality_factor = mortality_factor[seqindex] 
             
             incidence$index = index 
             
             incidence
          }

end_time <- Sys.time()
duration <- end_time - start_time


bootstrapincidence <- bind_rows(bootstrap_incidence)


write.csv(bootstrapincidence, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013mort_sensitivity_boot")



bootstrap_pivoted <-  bootstrapincidence  %>% 
  group_by(age, dates,  index, mortality_factor) %>%
  dplyr::select( incidencem) %>%
  pivot_wider(names_from = c(dates), 
              values_from = c(incidencem)) 

Incidence_diff <- bootstrap_pivoted %>%
  dplyr::select(age, index, mortality_factor, `2013`, `2018`) %>%
  transmute(age, index,mortality_factor,
            Mahiane_diff = `2018` - `2013`) %>% 
  group_by(age, mortality_factor) %>% 
  mutate(Mahiane_diffse = sd(Mahiane_diff))




write.csv(Incidence_diff, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/220182013mort_sensitivity_sd")


######################################
#  incidence  difference estimates 
######################################

inc_estimate <- list()

for (ii in seq_along(mortality_factor)){

  excessmortality <-  excess_mortality %>% 
    transmute(age = age,
              dates = as.numeric(as.character(dates)),
              mortality = mortality * mortality_factor[ii])
  
  
  incestimate <- sensitivity_mortality(msfdata_summary = msfdata)
  
  incestimate$mortality_factor = mortality_factor[ii] 
  
  
  inc_estimate[[ii]] = incestimate
}

incestimate <- bind_rows(inc_estimate)

incestimate_pivoted <-  incestimate  %>% 
  group_by(age, dates,mortality_factor) %>%
  dplyr::select( incidencem) %>%
  pivot_wider(names_from = c(dates), 
              values_from = c(incidencem)) 

Incidence_diff_estimate <- incestimate_pivoted %>%
  dplyr::select(age, mortality_factor, `2013`, `2018`) %>%
  transmute(age, mortality_factor,
            Mahiane_diff_est = `2018` - `2013`) 


############################################
#difference
############################################


populationest <- read.csv("/home/laurette/Desktop/Github/MSF_data/poulation_est.csv")%>%
  mutate(age  = as.numeric(age)) 


delt = 2; centralage <- 27

    population_est <- populationest %>% 
      dplyr::select(age, year,  males) %>% 
      filter(year == 2013) %>% 
      transmute(age = as.numeric(as.character(age)),
                population = males)%>% 
      filter(age >= centralage - delt & age <= centralage + delt)
    
    total_pop <- sum(population_est$population)
    
    bootstrap_pivoted <-  Incidence_diff  %>% 
      filter(age >= centralage - delt & age <= centralage + delt) %>% 
      group_by(age,  index, mortality_factor) %>%
      dplyr::select(Mahiane_diff) %>% 
      inner_join(population_est) 
    
    inc_diff <- bootstrap_pivoted %>%
      group_by(index, age,mortality_factor) %>%
      inner_join(Incidence_diff_estimate) %>% 
      transmute(age, index, mortality_factor,
                Mahianeinc_diff_est = (Mahiane_diff_est*population)/total_pop,
                Mahianeinc_diff = (Mahiane_diff*population)/total_pop)


    Mahianeeinc_diff_estimates <-  inc_diff %>%
      group_by(index, mortality_factor) %>%
      summarise(Mave_est = sum(Mahianeinc_diff, na.rm = T),
                Mave_estse = sum(Mahianeinc_diff, na.rm = T)) %>% 
      group_by(mortality_factor) %>% 
      summarise(M_ave = mean(Mave_est),
                M_avese = sd(Mave_estse),
                M_pvalue = 2*pnorm(-abs(M_ave/M_avese))) %>% 
      mutate(group = "25-29")
    


    
    newdata <- bootstrap_pivoted %>%
      filter(age == 27) %>% 
      group_by(mortality_factor) %>% 
  summarise(M_ave = mean(Mahiane_diff),
            M_avese = sd(Mahiane_diff),
            M_pvalue = 2*pnorm(-abs(M_ave/M_avese)),
            group = "27") 


bind_rows(Mahianeeinc_diff_estimates, newdata) %>% 
ggplot()+
  geom_line(aes(x = mortality_factor, y = M_pvalue, color = group)) +
  labs(x = "mortality factor", y = "p value ", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, 0.1))+
  theme(legend.position = "bottom")

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/plots20211119/mvarying_excess_mort.png", w = 15, h = 10)

write.csv(bind_rows(Mahianeeinc_diff_estimates, newdata), "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/pvaluemort_sensitivity_sd")
