

bootstrap_data <- function(data_2013  = msfdata2013_split, 
                            data_2018 = msfdata2018_split,
                            bootstrap_strategy ){
  
  
  if (bootstrap_strategy == "2013"){
  ###################################### SAMPLE 2013 ####################
  
  sampledclusters <- lapply(1:nward, function(x) sample(x = clusters_in_ward[[x]], 
                                                        size = length(clusters_in_ward[[x]]),
                                                        replace = T, prob = NULL))
  
  msfdata2013_splitcluster <- lapply(1:nward, function(x) 
    lapply(1:length(sampledclusters[[x]]), function(i) msfdata2013_split[[x]] %>%
             filter(cluster == sampledclusters[[x]][i])))
  
  bootsample <-  bind_rows(msfdata2013_splitcluster)
  
  
  } else if (bootstrap_strategy == "2018"){
  
  ############################################ SAMPLE 2018 ###############
  
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
  
  bootsample <- bind_rows(sampled_households)
  
  } else {
    
    ###################################### SAMPLE 2013 ####################
    
    sampledclusters <- lapply(1:nward, function(x) sample(x = clusters_in_ward[[x]], 
                                                          size = length(clusters_in_ward[[x]]),
                                                          replace = T, prob = NULL))
    
    msfdata2013_splitcluster <- lapply(1:nward, function(x) 
      lapply(1:length(sampledclusters[[x]]), function(i) msfdata2013_split[[x]] %>%
               filter(cluster == sampledclusters[[x]][i])))
    
    msf2013 <-  bind_rows(msfdata2013_splitcluster)
    
    
    ############################################ SAMPLE 2018 ###############
    
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
    
   bootsample <- bind_rows(msf2013, msf2018)
   }
  
  return(bootsample)
}


# analysis_data <- bootstrap_data(data_2013  = msfdata2013_split, 
#                                 data_2018 = msfdata2018_split,
#                                 bootstrap_strategy = "2013")


# analysis
# all parametric adjustments in the global environment

analyse_data <- function(msfdata_summary  = analysis_data){
  

  msfdata_onesvy <- msfdata_summary%>%
    mutate(age_years = as.integer(age_years)) %>%
    group_by(dates.x, age_years) %>%
    summarise(n_specimens = n(),
              positive = sum(hiv_status, na.rm = T),
              negative  = n_specimens - positive,
              recent  = sum(recent , na.rm = T),
              nonrecent = positive - recent) %>%
    mutate(age = age_years, overallweight = 1)


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
    extractprevmultitaylor3(prevalencefittedmodel = prevalencefittedmodel,
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
  

  incestimate <- prev_estimated_prevalencedata
}








analyse_separate_recmodels <- function(msfdata_summary = analysis_data){
  
  
  msfdata_onesvy <- msfdata_summary%>%
    mutate(age_years = as.integer(age_years)) %>%
    group_by(dates.x, age_years) %>%
    summarise(n_specimens = n(),
              positive = sum(hiv_status, na.rm = T),
              negative  = n_specimens - positive,
              recent  = sum(recent , na.rm = T),
              nonrecent = positive - recent) %>%
    mutate(age = age_years, overallweight = 1)
  
  
  filter_data <-   inclusionexclusion(ageradius = previnclusion_dist, age_0 = anchorage,
                                      agetimetrans =  agetimetrans,
                                      exponentterm = exponentterm, points_weights = points_weights,
                                      timeradius = timecutoff, recency = recency,
                                      time_0 = midpoint,
                                      inclusion_method = inclusionmethod,
                                      glmformatteddata = msfdata_onesvy)
  
  
  #separate the combined data
  filter_data2013 <- filter_data %>% 
    filter(time_param == -2.5)
  
  filter_data2018 <- filter_data %>% 
    filter(time_param == 2.5)
  
  
  
  # fit models prevalence and recency
  prevalencefittedmodel <-  tryCatch(fitprevmultitaylor(glmformatteddata = filter_data,
                                                        link = prevlink,
                                                        prevtaylororder = prevtaylororder),
                                     error = function(e) NA)
  
  
  filterdata <- list(filter_data2013, filter_data, filter_data2018)
  
  
  
  recentfittedmodel <-  lapply(seq_along(filterdata),
                               function(x) tryCatch(fitrecmultitaylor(glmformatteddata = filterdata[[x]],
                                                                      link = reclink,
                                                                      rectaylororder = rectaylororder),
                                 error = function(e) NA))
  
  predictdataset <- lapply(seq_along(filterdata),
                           function(x) newdata <- data.frame(age_param = age_predict - anchorage,
                                                             time_param = predicttimes[x] - midpoint))
  
  
  prevR_predictions <- lapply(seq_along(recentfittedmodel),
                                  function(x) stats::predict(recentfittedmodel[[x]], 
                                                             newdata = predictdataset[[x]], 
                                                             type = "response", se.fit = TRUE))
  
  recencydata <- data.frame(dates = sort(rep(predicttimes, length(age_predict))),
                            age = age_predict,
                            prevrecency = unlist(lapply(seq_along(prevR_predictions), 
                                                          function (x) prevR_predictions[[x]]$fit)),
                            prevrecency_stderror = unlist(lapply(seq_along(prevR_predictions), 
                                                          function (x) prevR_predictions[[x]]$se.fit)))
  
  
  # extract the prevalences for Mahiane
  prev_estimatedprevalencedata <- bind_rows(lapply(seq_along(predicttimes), function(x) tryCatch(
    extractprevmultitaylor3(prevalencefittedmodel = prevalencefittedmodel,
                            recencyfittedmodel = NULL,
                            agetimetrans = agetimetrans, epsilon = epsilon,
                            anchorage = anchorage, anchortime = midpoint,
                            prevtaylororder = prevtaylororder,
                            age_predict = age_predict,  #oneparam = oneparam,
                            time_predict = predicttimes[x]),
    error = function(e) NA)))
  
  
  prev_estimated_prevalencedata <- prev_estimatedprevalencedata %>%
    inner_join(excessmortality, by = c("age", "dates")) %>% 
    inner_join(recencydata, by = c("age", "dates"))
  
  
  prev_estimated_prevalencedata$incidencem <- mincidence_pe(prev = prev_estimated_prevalencedata$prevalence,
                                                            slope = prev_estimated_prevalencedata$prevslope,
                                                            excess_mortality = prev_estimated_prevalencedata$mortality)
  
  
  incidenceK <- tryCatch(incprops(PrevH = prev_estimated_prevalencedata$prevalence, PrevR = prev_estimated_prevalencedata$prevrecency,
                                  RSE_PrevH = prev_estimated_prevalencedata$prevalence_stderror/ prev_estimated_prevalencedata$prevalence,
                                  RSE_PrevR = prev_estimated_prevalencedata$prevrecency_stderror / prev_estimated_prevalencedata$prevrecency,
                                  MDRI = ifelse(prev_estimated_prevalencedata$dates ==2013, 217,  207),
                                  FRR = ifelse(prev_estimated_prevalencedata$dates ==2013, 0.001667, 0.0018),
                                  RSE_MDRI = RSE_MDRI, RSE_FRR = RSE_FRR, BigT = BigT), error = function(e) NULL)
  
  prev_estimated_prevalencedata$incidencek <- as.numeric(levels(incidenceK$Incidence.Statistics$Incidence))[incidenceK$Incidence.Statistics$Incidence]
  

  
  ###################################################### Estimates ######################################################
  
  incestimate <- prev_estimated_prevalencedata
}








naive_estimates <- function(msfdata_summary = analysis_data){
  
  msfdata_onesvy <- msfdata_summary%>%
    mutate(age_years = as.integer(age_years)) %>%
    group_by(dates.x, age_years) %>%
    summarise(n_specimens = n(),
              positive = sum(hiv_status, na.rm = T),
              negative  = n_specimens - positive,
              recent  = sum(recent , na.rm = T),
              nonrecent = positive - recent)
  
  filterdata <- msfdata_onesvy %>% 
    dplyr::select(age = age_years, dates = dates.x,
                  n_specimens, positive, recent) %>% 
    # filter(age>=20 & age <= 30) %>%
    summarise(n_specimens = sum(n_specimens),
              positive = sum(positive),
              recent = sum(recent))
  
  
  prev_binomestimates <- binom.confint(filterdata$positive, 
                                       filterdata$n_specimens,
                                       methods = "exact")
  
  names(prev_binomestimates)<- c("method ", "positive ", "n_specimens", 
                                 "prevalence", "prev_lower", "prev_upper")
  
  rec_binomestimates <- binom.confint(filterdata$recent, 
                                       filterdata$positive,
                                       methods = "exact")
  
  names(rec_binomestimates)<- c("method ", "recent", "positve",
                                "prev_recent", "rec_lower", 
                                "rec_upper")
  
  prevdata <- bind_cols(filterdata[,c(1, 2)],
                        prev_binomestimates[, 4:6],
                        rec_binomestimates[, c(4:6)]) %>% 
    mutate(prev_se = (prev_upper - prev_lower)/3.92,
           rec_se = (rec_upper - rec_lower)/3.92) 
  
  
  incidencedata <- inctools::incprops(PrevH = prevdata$prevalence, RSE_PrevH = (prevdata$prev_se)/prevdata$prevalence, 
                                      PrevR = prevdata$prev_recent, RSE_PrevR = (prevdata$rec_se)/prevdata$prev_recent, Boot = FALSE,
                                      BS_Count = 10000, alpha = 0.05, BMest = "same.test",
                                      MDRI = ifelse(prevdata$dates ==2013, 217,  207),
                                      FRR = ifelse(prevdata$dates ==2013, 0.001667, 0.0018),
                                      RSE_MDRI = 0, RSE_FRR = 0, BigT = 730, Covar_HR = 0)
  
  
 prevdata$incidence = as.numeric(as.vector(incidencedata$Incidence.Statistics[,2]))
  
   ########################### Estimates #################################
 
 incestimate <- prevdata
}
