
mortalitydata <- read_feather("/home/laurette/Desktop/Github/MSF_data/incidence-combined-method-master/mortality_table.feather") 

excess_mortality <-  mortalitydata %>% 
  mutate(`2018` = `2014` + ((`2014` - `2015`)/(2014 - 2015)) * (2018 - 2014)) #%>%
  melt(id.vars = c("age","sex")) %>%
  filter(age >= 30)

  
  write.csv(excess_mortality, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/mortality2018")
  