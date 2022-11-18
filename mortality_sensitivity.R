sensitivity_mortality <- function(msfdata_summary = analysis_data){
  
  
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
  
  # separate the combined data
  # fit models prevalence and recency
  
  prevalencefittedmodel <-  tryCatch(fitprevmultitaylor(glmformatteddata = filter_data,
                                                        link = prevlink,
                                                        prevtaylororder = prevtaylororder),
                                     error = function(e) NA)

  
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
    inner_join(excessmortality, by = c("age", "dates")) 
  
  prev_estimated_prevalencedata$incidencem <- mincidence_pe(prev = prev_estimated_prevalencedata$prevalence,
                                                            slope = prev_estimated_prevalencedata$prevslope,
                                                            excess_mortality = prev_estimated_prevalencedata$mortality)
  
  ###################################################### Estimates ######################################################
  
  incestimate <- prev_estimated_prevalencedata
}
