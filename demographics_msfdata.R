source('~/Desktop/Github/msfanalysis20211110/MSFDataCleaning.R')


msf2013 = msfdata2013_original %>% 
  mutate(age = cut(age_years, c(15,20,25, 30, 35, 40, 45,50, 100), right = F)) %>% 
  group_by(age, sex) %>% 
  summarise(total = n(),
            negative = total - sum(hiv_status),
            positive = sum(hiv_status),
            recent = sum(recent, na.rm = T))
         


msf2018 = msfdata2018_mod %>%
  mutate(age = cut(age_years, c(15,20,25, 30, 35, 40, 45,50, 60, 100), right = F)) %>% 
  group_by(age, sex) %>% 
  summarise(total = n(),
            negative = total - sum(hiv_status, na.rm = T),
            positive = sum(hiv_status, na.rm = T),
            recent = sum(recent, na.rm = T))



msfdata2013_original %>% 
  group_by(age_years, sex) %>% 
  summarise(total = n()) %>% 
  pivot_wider()





