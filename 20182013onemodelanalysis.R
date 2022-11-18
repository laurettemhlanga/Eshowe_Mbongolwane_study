rm(list = ls())
library(data.table)
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
  filter(Sex == "Female") %>%
  transmute(age = Age, sex = Sex,
            X2013 = X2013,
            X2015 = X2015,
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
  prevtaylororder = 3; epsilon = 0.01; rectaylororder = 3
  previnclusion_dist = 20; recinclusion_dist = 20;  
  anchorage <- 25; age_predict <- 15:40; oneparam = F; 
  agetimetrans = F; len_anchorage <- length(anchorage);
  reclink = "cloglog"; recinclusion_dist = 10; MDRI = 207; FRR = 0.0018
  RSE_MDRI = 0.01 ; RSE_FRR = 0.25; 
  epsi <- (MDRI/365.25)/2; BigT = 730
  midpoint <- 2015.5; predicttimes <- c(2013, 2015.5, 2018)
  
  ######################################
  
  ######################################
  # spliting the survey data
  ######################################
  surveyiterations <- 1000
  
  clustersin2013 <- unique(msfdata[msfdata$dates.x == 2013,]$cluster);
  clustersin2018 <- unique(msfdata[msfdata$dates.x == 2018,]$cluster)
  
  ncluster13 = length(clustersin2013)
  ncluster18 = length(clustersin2018)
  nclusters = ncluster18 + ncluster13
  nward = 14
  
  
  msfdata2013mod <- msfdata2013_mod %>% 
    filter(sex == "Female")
    
  
  msfdata2018mod <- msfdata2018_mod %>% 
    filter(sex == "Female")
  
  # preprocess 2013
  msfdata2013_split <- split(msfdata2013mod, list(msfdata2013mod$ward))
  
  clusters_in_ward <- lapply(seq_along(msfdata2013_split), 
                             function(x) unique(msfdata2013_split[[x]]$cluster))
  
  #preprocess 2018
  msfdata2018_split <- split(msfdata2018mod, list(msfdata2018mod$cluster))


  ######################################
  #  Bootstraps 
  ######################################
  
  
  registerDoParallel(detectCores() - 3)
  getDoParWorkers()
  
  start_time <- Sys.time()

  bootstrap_incidence <- foreach(index = seq_along(1:surveyiterations), 
                                 .combine = "rbind") %dopar% {
  
  ###################################### SAMPLE 2013 ################################################################
  
    sampledclusters <- lapply(1:nward, function(x) sample(x = clusters_in_ward[[x]], 
                                                          size = length(clusters_in_ward[[x]]),
                                                          replace = T, prob = NULL))
    
    msfdata2013_splitcluster <- lapply(1:nward, function(x) 
      lapply(1:length(sampledclusters[[x]]), function(i) msfdata2013_split[[x]] %>%
               filter(cluster == sampledclusters[[x]][i])))
    
    msf2013 <-  bind_rows(msfdata2013_splitcluster)
    
    ############################################ SAMPLE 2018 ##########################################################
  
    sampledclusters <-  sample(x = 1:ncluster18, size = ncluster18, replace = T, prob = NULL)
  
    unique_household <- lapply(1:ncluster18, function(x) 
      unique(msfdata2018_split[[sampledclusters[x]]]$household_id))
    
    # sample hhlds 
    sampledhouseholds <- lapply(1:ncluster18, function(x)
      sample(x = unique_household[[x]], size = length(unique_household[[x]]), 
             replace = T, prob = NULL))
    
    sampled_households <- lapply(1: ncluster18, function(x)
      lapply(1:length(sampledhouseholds[[x]]), function(i) msfdata2018_split[[sampledclusters[x]]] %>% 
               filter(household_id == sampledhouseholds[[x]][i])))
    
    msf2018 <- bind_rows(sampled_households)
    
    ############################################ COMBINE DATA ##########################################################
    
    msfdata_summary <- bind_rows(msf2018, msf2013) %>%
      mutate(age_years = as.integer(age_years)) %>% 
      group_by(dates.x, age_years) %>% 
      summarise(n_specimens = n(),
                positive = sum(hiv_status, na.rm = T), 
                negative  = n_specimens - positive,
                recent  = sum(recent , na.rm = T),
                nonrecent = positive - recent) %>%
      mutate(age = age_years, overallweight = 1)
    
    
    msfdata_onesvy <- msfdata_summary
    
    ###################################################### Dynamic ######################################################
    
    filter_data <-   inclusionexclusion(ageradius = previnclusion_dist, age_0 = anchorage,
                                        agetimetrans =  agetimetrans,
                                        exponentterm = exponentterm, points_weights = points_weights,
                                        timeradius = timecutoff, recency = recency,
                                        time_0 = midpoint,
                                        inclusion_method = inclusionmethod,
                                        glmformatteddata = msfdata_onesvy)
    
    # fit models prevalence and recency 
    prevalencefittedmodel <-  tryCatch(fitprevmultitaylor(glmformatteddata = filter_data,
                                                          link = prevlink, 
                                                          prevtaylororder = prevtaylororder), 
                                       error = function(e) NA)
    
    
    
    recentfittedmodel <-  tryCatch(fitrecmultitaylor(glmformatteddata = filter_data,
                                                     link = reclink, 
                                                     rectaylororder = rectaylororder), 
                                   error = function(e) NA)
    
    
    
    # extract the prevalences for Mahiane 
    prev_estimatedprevalencedata <- bind_rows(lapply(seq_along(predicttimes), function(x) tryCatch(
      extractprevmultitaylor2(prevalencefittedmodel = prevalencefittedmodel, 
                              recencyfittedmodel = recentfittedmodel,
                              agetimetrans = agetimetrans, epsilon = epsilon,
                              anchorage = anchorage, anchortime = midpoint,
                              prevtaylororder = prevtaylororder,
                              age_predict = age_predict,  #oneparam = oneparam,
                              time_predict = predicttimes[x]),
      error = function(e) NA)))
    
    
    prev_estimated_prevalencedata <- prev_estimatedprevalencedata %>% 
      inner_join(excessmortality, by = c("age", "dates"))
    
    
    prev_estimated_prevalencedata$incidencem <- mincidence_pe(prev = prev_estimated_prevalencedata$prevalence, 
                                                              slope = prev_estimated_prevalencedata$prevslope,
                                                              excess_mortality = prev_estimated_prevalencedata$mortality)
    
    
    incidenceK <- tryCatch(incprops(PrevH = prev_estimated_prevalencedata$prevalence, PrevR = prev_estimated_prevalencedata$prevrecency,
                                    RSE_PrevH = prev_estimated_prevalencedata$prevalence_stderror/ prev_estimated_prevalencedata$prevalence,
                                    RSE_PrevR = prev_estimated_prevalencedata$prevrecency_stderror / prev_estimated_prevalencedata$prevrecency,
                                    MDRI = MDRI, RSE_MDRI = RSE_MDRI, FRR =  FRR,
                                    RSE_FRR = RSE_FRR, BigT = BigT), error = function(e) NULL)
    
    prev_estimated_prevalencedata$incidencek <- as.numeric(levels(incidenceK$Incidence.Statistics$Incidence))[incidenceK$Incidence.Statistics$Incidence]
    
    prev_estimated_prevalencedata$index = index
    ###################################################### Estimates ######################################################
    
    incestimate <- prev_estimated_prevalencedata
    
    }
  
  end_time <- Sys.time()
  overalltime =  end_time - start_time 
  
  
write.csv(bootstrap_incidence, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/20182013bootsonemodel_sd")
  
  
  bootstrap_sd <- bootstrap_incidence %>%
    group_by(age, dates) %>% 
    summarise(covariance = cov(matrix(incidencem, incidencek, nrow = n(), ncol = 2))[2,1],
              prevse = sd(prevalence),
              mse = sd(incidencem),
              kse = sd(incidencek),
              correlation = covariance/ (mse *kse),
              weight = (kse^2 - covariance)/(mse^2 + kse^2 - 2 * covariance),
              optse = sqrt((weight * mse)^2 + ((1 - weight) * kse)^2 + 
                             (2 * covariance*(weight *(1 - weight))))) 
    
  
write.csv(bootstrap_sd, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/20182013kznonemodel_sd")


##############################################################
# difference 
##############################################################
  
  
bootstrap_estimates<- bootstrap_incidence %>%
  group_by(age, dates) %>% 
  mutate(covariance = cov(matrix(incidencem, incidencek, nrow = n(), ncol = 2))[2,1],
            prevse = sd(prevalence),
            mse = sd(incidencem),
            kse = sd(incidencek),
            weight = (kse^2 - covariance)/(mse^2 + kse^2 - 2 * covariance),
            incidenceo = (weight *incidencem) + (1 - weight) *incidencek,
            optse = sqrt((weight * mse)^2 + ((1 - weight) * kse)^2 + 
                           (2 * covariance*(weight *(1 - weight)))))%>% 
  select(age, dates, index, incidencek, incidencem, incidenceo) %>% 
  filter(dates != 2015.5)%>%
  group_by(age, dates,  index) %>%
  select( incidencem, incidencek, incidenceo) %>%
  pivot_wider(names_from = c(dates), values_from = c(incidencem, incidencek, incidenceo)) 


incidence_diff <- bootstrap_estimates %>%
  mutate(Mahiane = incidencem_2018 - incidencem_2013, 
         Kassanjee = incidencek_2018 - incidencek_2013,
         Weighted = incidenceo_2018 - incidenceo_2013) %>% 
  group_by(age) %>% 
  summarise(mse = sd(Mahiane, na.rm = T),
            kse = sd(Kassanjee, na.rm = T),
            wse = sd(Weighted, na.rm = T))

