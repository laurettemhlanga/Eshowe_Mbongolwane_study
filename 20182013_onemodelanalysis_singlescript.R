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
prevtaylororder = 3; epsilon = 0.01; rectaylororder = 1
previnclusion_dist = 20; recinclusion_dist = 20;  
anchorage <- 20; age_predict <- 15:40; oneparam = F; 
agetimetrans = F; len_anchorage <- length(anchorage);
reclink = "cloglog"; MDRI = 207; FRR = 0.0018
RSE_MDRI = 0.01 ; RSE_FRR = 0.25; 
epsi <- (MDRI/365.25)/2; BigT = 730
midpoint <- 2015.5; predicttimes <- c(2013, 2015.5,2018)

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

bootstrap_strategy = "both"

msfdata2013mod <- msfdata2013_original %>% 
  filter(sex == "Female")


msfdata2018mod <- msfdata2018_mod %>% 
  filter(sex == "Female")

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


bootstrap_incidence <- foreach(index = seq_along(1:surveyiterations), 
                               .combine = "rbind") %dopar% {
                                 
           analysis_data <- bootstrap_data(data_2013  = msfdata2013_split, 
                                           data_2018 = msfdata2018_split,
                                           bootstrap_strategy = bootstrap_strategy)
           
           incidence <- analyse_separate_recmodels(msfdata_summary  = analysis_data)
           
           incidence$index = index 
           
           incidence
        }

end_time <- Sys.time()
duration <- end_time - start_time

write.csv(bootstrap_incidence, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013bootsonemodel_sd")


sde <- bootstrap_incidence %>%
  group_by(age, dates) %>% 
  summarise(covariance = cov(matrix(incidencem,
                                    incidencek, 
                                    nrow = n(), ncol = 2))[2,1],
            prevse = sd(prevalence),
            mse = sd(incidencem),
            kse = sd(incidencek),
            correlation = covariance/ (mse *kse),
            weight = (kse^2 - covariance)/(mse^2 + kse^2 - 2 * covariance),
            optse = sqrt((weight * mse)^2 + ((1 - weight) * kse)^2 + 
                           (2 * covariance*(weight *(1 - weight))))) 

write.csv(sde, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013kznonemodel_sd")

# sde <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013kznonemodel_sd")

weightdata <- sde %>% 
  group_by(age, dates) %>% 
  transmute(weight = weight)


######################################
#  incidence estimates 
######################################


incestimate <- analyse_separate_recmodels(msfdata_summary = msfdata)


incestimate_prev <- incestimate %>% 
  select( age,  dates, prevalence)

write.csv(incestimate_prev, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/prev_est")


inc_estimate <- incestimate %>% 
  inner_join(weightdata) %>% 
  mutate(incidenceo = (incidencem*weight) + 
           ((1 - weight)*incidencek )) %>% 
  select(age, dates, incidencem,  incidencek, incidenceo)

write.csv(inc_estimate, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013kznonemodel_est")
inc_estimate <- read.csv()

incidence_long <- melt(inc_estimate, id.vars = c("dates","age")) 


sde_mod <- sde[,-c(3,4, 8,7)]
names(sde_mod) <- c( "age", "dates", "incidencem", "incidencek","incidenceo")
sde1 <- melt(sde_mod, id.vars = c("dates","age")) %>% 
  select(dates, age, variable, se = value )

incidencelong <- inner_join(x = incidence_long, y = sde1)




newdata <- incidencelong %>% 
  filter(dates != 2015.5)

ggplot()+
  geom_ribbon(data = newdata, aes(x = age, ymin = ifelse(value - (1.96 * se)<0, 0, value - (1.96 * se)),
                                  ymax = value + (1.96 * se), fill = variable), alpha = 0.3) +
  geom_line(data = newdata, aes(x = age, y = value , colour= variable)) +
  labs(x = "age", y = "estimate ", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  facet_grid(cols = vars(dates))+
  scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
  scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,40), breaks = seq(15,40, 5))+
  # scale_y_continuous(limits=c(0,0.08), breaks = seq(0,0.08,0.02))+
  theme(legend.position = "bottom")

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/maleplots/allplots_zoomedin30_20.png", w = 15, h = 7)



############################################
#difference
############################################


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

names(incidence_diff) <- c("age", "Mahiane", "Kassanjee", "Weighted")

diff_sde_long <- melt(incidence_diff, id.vars = "age")
names(diff_sde_long) <- c("age", "variable", "se")

diff_estimates <- inc_estimate %>%
  select(age,  dates, incidencem,  incidencek, incidenceo) %>% 
  filter(dates != 2015.5)%>%
  group_by(age, dates) %>%
  select( incidencem, incidencek, incidenceo) %>%
  pivot_wider(names_from = c(dates), values_from = c(incidencem, incidencek, incidenceo)) 


incidence_diff <- diff_estimates %>%
  transmute(Mahiane = incidencem_2018 - incidencem_2013, 
            Kassanjee = incidencek_2018 - incidencek_2013,
            Weighted = incidenceo_2018 - incidenceo_2013) %>% 
  melt(id.vars = "age") %>% 
  inner_join(diff_sde_long)

# write.csv(incidence_diff, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_est")


ggplot()+
  geom_ribbon(data = incidence_diff, aes(x = age, ymin = value - (1.96 * se),
                                         ymax = value + (1.96 * se), fill = variable), alpha = 0.3) +
  geom_line(data = incidence_diff, aes(x = age, y = value , colour= variable)) +
  labs(x = "age", y = "incidence difference ", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  geom_hline(aes(yintercept = 0), linetype ="dashed", color = "black")+
  scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
  scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,30), breaks = seq(16,30, 2))+
  scale_y_continuous(limits=c(-0.08,0.08), breaks = seq(-0.06,0.06,0.02))+
  theme(legend.position = "bottom")


ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/maleplots/incdifference_stimates230_20.png", 
       w=15,h= 10)






