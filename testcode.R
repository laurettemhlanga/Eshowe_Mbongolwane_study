# source('~/Desktop/Github/MSF_data/MSF_datacleaning.R')

# rec_estimated_prevalencedata %>% 
#   select(age, value = prevrecency, se = prevrecency_stderror) %>% 
#   ggplot()+
#   geom_ribbon( aes(x = age, ymin = value - (1.96 * se), 
#                                   ymax = value + (1.96 * se)), alpha = 0.3) +
#   geom_line( aes(x = age, y = value )) 


fitdata <- msfdata2018_mod %>% 
  filter(age_years <= 40)

form <- formula(recent ~ I(age_years) )

form0 <- formula(recent ~ I(age_years) + I(age_years^2))

form1 <- formula(recent ~ I(age_years) + I(age_years^2) +I(age_years^3))

form2 <- formula(recent ~ I(age_years) + I(age_years^2) +I(age_years^3) + I(age_years^4))

formulas <- list(form, form0,form1, form2)

fittedmodel <- list()


for (ii in seq_along(formulas)){



model <- glm2(formula = formulas[[ii]],
              data = fitdata,
              family = binomial(link = "log"))

predictdataset <- data.frame(age_years = 15:40)

# prevR_predictions <- stats::predict(model, newdata = predictdataset,
#                                     type = "response", se.fit = TRUE)


fittedmodel[[ii]] <- model

}

lapply(seq_along(fittedmodel), function(x) summary(fittedmodel[[x]]))

# data.frame(age = 15:40,
#            value = prevR_predictions$fit,
#            se = prevR_predictions$se.fit) %>% 
#   ggplot()+
#   geom_ribbon( aes(x = age, ymin = value - (1.96 * se), 
#                    ymax = value + (1.96 * se)), alpha = 0.3) +
#   geom_line( aes(x = age, y = value )) 
# 
# qplot(15:59,prevR_predictions$fit)+
#   ylim(0,0.2)+


msfdata2013mod <- msfdata2013_original %>% 
  filter(sex == "Female")


msfdata2018mod <- msfdata2018_mod %>% 
  filter(sex == "Female")

msfdata <- bind_rows(msfdata2013mod, msfdata2018mod)

fitdata <- msfdata %>% 
  filter(age_years <= 40)

form1 <- formula(hiv_status ~ I(age_years) + I(dates.x))

form2 <- formula(hiv_status ~ I(age_years) + I(age_years^2) +  I(dates.x) + I(age_years *dates.x))

form3 <- formula(hiv_status ~ I(age_years) + I(age_years^2) +I(age_years^3)+ I(dates.x) +
                   I(age_years *dates.x) + I(age_years^2 *dates.x))

form4 <- formula(hiv_status ~ I(age_years) + I(age_years^2) +I(age_years^3) + I(age_years^4)+
                   I(dates.x) + I(age_years *dates.x)+ I(age_years^2 *dates.x) + I(age_years^3 *dates.x))

form5 <- formula(hiv_status ~ I(age_years) + I(age_years^2) +I(age_years^3) + I(age_years^4)+
                   I(dates.x) + I(age_years *dates.x)+ I(age_years^2 *dates.x) + I(age_years^3 *dates.x)+
                   I(age_years^5)+ I(age_years^4 * dates.x))

formulas <- list(form1, form2,form3, form4, form5)

fittedmodel <- list()


for (ii in seq_along(formulas)){
  
  
  
  model <- glm2(formula = formulas[[ii]],
                data = fitdata,
                family = binomial(link = "logit"))
  
  predictdataset <- data.frame(age_years = 15:40)
  
  # prevR_predictions <- stats::predict(model, newdata = predictdataset,
  #                                     type = "response", se.fit = TRUE)
  
  
  fittedmodel[[ii]] <- model
  
}

lapply(seq_along(fittedmodel), function(x) summary(fittedmodel[[x]]))
