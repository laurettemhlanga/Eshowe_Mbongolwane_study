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

################################################################################################################
# uses the data provided by Eduard - kznfeather
# REQUIRED uses the new sampling strategy to estimate the 
# mid point and survey times incidence estimates.
# fits one model on all the data from 15:45 
# the link function uses the logit for prevalence and 
# clog log for the prevalence of recency function.
# anchorage is set at 25 to ensure use of all data points between 15:35


######################################
# input parameters
######################################


previnclusion_dist = 5; recinclusion_dist = 5;  
anchorage <- 25; 
######################################

######################################
# spliting the survey data
######################################
surveyiterations <- 10000

msfdata2013mod <- msfdata2013_original %>% 
  filter(sex == "Male" ) 



msfdata2018mod <- msfdata2018_mod %>% 
  filter(sex == "Male")

clustersin2013 <- unique(msfdata2013mod[msfdata2013mod$dates.x == 2013,]$cluster);
clustersin2018 <- unique(msfdata2018mod[msfdata2018mod$dates.x == 2018,]$cluster)

ncluster13 = length(clustersin2013)
ncluster18 = length(clustersin2018)
nclusters = ncluster18 + ncluster13
nward = 14

bootstrap_strategy = "both"



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


bootstrap_incidence <- foreach(index = seq_along(1:surveyiterations)) %dopar% {
                                 
                                 analysis_data <- bootstrap_data(data_2013  = msfdata2013_split, 
                                                                 data_2018 = msfdata2018_split,
                                                                 bootstrap_strategy = bootstrap_strategy)
                                 
                                 incidence <- tryCatch(naive_estimates(msfdata_summary  = analysis_data),
                                                       error = function(e) NA)
                                 
                                 # incidence$index = ifelse(is.na(incidence), NA, index) 
                                 
                                 incidence
                               }

end_time <- Sys.time()
duration <- end_time - start_time


incidence_pnt <- bind_rows(bootstrap_incidence[(which(!is.na(bootstrap_incidence)))]) %>% 
  group_by(dates) %>% 
  summarise(estimate = mean(incidence),
            se = sd(incidence),
            lower  = estimate - (1.96 *se),
            upper  = estimate + (1.96 *se))

write.csv(incidence_pnt, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/naivepooling_estimates")


# sde <- bootstrap_incidence %>%
#   group_by(age, dates) %>% 
#   summarise(covariance = cov(matrix(incidencem,
#                                     incidencek, 
#                                     nrow = n(), ncol = 2))[2,1],
#             prevse = sd(prevalence),
#             mse = sd(incidencem),
#             kse = sd(incidencek),
#             correlation = covariance/ (mse *kse),
#             weight = (kse^2 - covariance)/(mse^2 + kse^2 - 2 * covariance),
#             optse = sqrt((weight * mse)^2 + ((1 - weight) * kse)^2 + 
#                            (2 * covariance*(weight *(1 - weight))))) 
# 
# write.csv(sde, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/naivepooling_se_sd")
# 
# # sde <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013kznonemodel_sd")
# 
# weightdata <- sde %>% 
#   group_by(age, dates) %>% 
#   transmute(weight = weight)
# 
# 
# ######################################
# #  incidence estimates 
# ######################################
# 
# 
# incestimate <- naive_estimates(msfdata_summary = msfdata)
# 
# 
