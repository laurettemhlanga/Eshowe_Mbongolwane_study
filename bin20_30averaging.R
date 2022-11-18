msfdata

source('~/Desktop/Github/MSF_data/MSF_datacleaning.R')

bootstrap_incidence <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013bootsonemodel20_30_sd")
  # read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013bootsonemodel20_30_sd")
   # read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013bootsonemodel_sd")

populationest <- msfdata %>%
  group_by(age_years, dates.x, sex) %>% 
  summarise(total = n()) %>%
  pivot_wider(names_from = c(sex), 
              values_from = c(total)) 


population_est <- populationest %>% 
  dplyr::select(age_years, dates = dates.x, Female) %>% 
  transmute(age = age_years,
            population = Female) %>% 
  filter(age >= 20 & age <= 30)

total_pop = sum(population_est$population)


bootstrap_pivoted <-  bootstrap_incidence  %>% 
  filter(dates !=2015.5) %>% 
  group_by(age, dates,  index) %>%
  dplyr::select( incidencem, incidencek) %>% 
  filter((age >=20 & age <= 30)) %>% 
  inner_join(population_est) %>% 
  group_by(index, dates) %>% 
  mutate(incidencem_se = sd(incidencem),
         incidencek_se = sd(incidencek),
         covariance = cov(matrix(incidencem, incidencek, 
                                 nrow = n(), ncol = 2))[2,1],
         weight = ifelse((incidencek_se^2 - covariance)/(incidencem_se^2 + 
                                                           incidencek_se^2 - (2 * covariance)) <0,0,
                         ifelse((incidencek_se^2 - covariance)/(incidencem_se^2 + 
                                                                  incidencek_se^2 - (2 * covariance)) >1 , 1,
                                (incidencek_se^2 - covariance)/(incidencem_se^2 + 
                                                                  incidencek_se^2 - (2 * covariance)))),
         weighted = (weight * incidencem) + (1 - weight) * incidencek) %>%  
  summarise(incidencem = sum((incidencem *population)/total_pop),
            incidencek = sum((incidencek *population)/total_pop),
            weighted = sum((weighted *population)/total_pop))#20182013bootsonemodel20_30_sd

estimates <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013kznonemodel20_30_est")

  data_est <- estimates %>% 
    filter(dates != 2015.5) %>% 
    inner_join(population_est) %>% 
    group_by(dates) %>% 
    summarise(incidencem = sum((incidencem *population)/total_pop),
              incidencek = sum((incidencek *population)/total_pop),
              incidenceo = sum((incidenceo *population)/total_pop))

  
  

bootstrap_pivoted%>% 
  group_by(dates) %>% 
  summarise(incidencem_se = sd(incidencem),
            incidencek_se = sd(incidencek),
            weight_se = sd(weighted)) %>% 
  inner_join(data_est)






