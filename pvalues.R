############################
###########estimates
############################
rm(list = ls())

bootstrap_incidence <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013bootsonemodel_sd")


bootstrap_pivoted <-  bootstrap_incidence  %>% 
  filter(dates !=2015.5) %>% 
  group_by(age, dates,  index) %>%
  dplyr::select( incidencem, incidencek) %>%
  pivot_wider(names_from = c(dates), 
              values_from = c(incidencem, incidencek)) 

Incidence_diff <- bootstrap_pivoted %>%
  dplyr::select(age, index, incidencem_2018, incidencem_2013,
                incidencek_2018, incidencek_2013) %>%
  transmute(age, index, Mahiane_diff = incidencem_2018 - incidencem_2013,
            Kassanjee_diff =  incidencek_2018 - incidencek_2013) %>% 
  group_by(age) %>% 
  mutate(Mahiane_diffse = sd(Mahiane_diff),
         Kassanjee_diffse = sd(Kassanjee_diff),
         covariance = cov(matrix(Mahiane_diff, 
                                 Kassanjee_diff, 
                                 nrow = n(), ncol = 2))[2,1],
         weight = (Kassanjee_diffse^2 - covariance)/(Mahiane_diffse^2 + 
                                                       Kassanjee_diffse^2 - (2 * covariance)),
         weighted_diff = (weight * Mahiane_diff) + (1 - weight) * Kassanjee_diff)


weight <- Incidence_diff %>% 
  dplyr::select(age,  weight, Kassanjee_diffse, 
                Mahiane_diffse, weighted_diff) %>%
  group_by(age) %>% 
  summarise(weight = mean(weight),
            Kassanjee_diffse = mean(Kassanjee_diffse), 
            Mahiane_diffse = mean(Mahiane_diffse),
            weighted_diffse = sd(weighted_diff))

estimate <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_est")

estimates <- estimate %>% 
  inner_join(weight) %>% 
  mutate(Weighted_diff  = (Mahiane *weight) +(1 - weight)*Kassanjee)

# write.csv(estimates, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_estse")

############################
###########estimates
############################


bootstrap_sd <- bootstrap_incidence %>%
  group_by(age, dates, index) %>% 
  filter(dates != 2015.5) %>% 
  dplyr::select(incidencem, incidencek) %>%
  pivot_wider(names_from = c(dates), 
              values_from = c(incidencem, incidencek)) %>%
  dplyr::select(age, index, incidencem_2018, incidencem_2013,
                incidencek_2018, incidencek_2013) %>%
  transmute(age, index, Mahiane_diff = incidencem_2018 - incidencem_2013,
            Kassanjee_diff =  incidencek_2018 - incidencek_2013) %>% 
  group_by(age) %>% 
  mutate(Mahiane_diffse = sd(Mahiane_diff),
         Kassanjee_diffse = sd(Kassanjee_diff),
         covariance = cov(matrix(Mahiane_diff, 
                                 Kassanjee_diff, 
                                 nrow = n(), ncol = 2))[2,1],
         weight = (Kassanjee_diffse^2 - covariance)/(Mahiane_diffse^2 + 
                                                       Kassanjee_diffse^2 - (2 * covariance)),
         weighted_diff = (weight * Mahiane_diff) + (1 - weight) * Kassanjee_diff)

# estimates <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_estse")

populationest <- read.csv("/home/laurette/Desktop/Github/MSF_data/poulation_est.csv")
###################################################################################################################
# delta t from one central age
###################################################################################################################

delta =0:5; central_age <- 15:40



iterationaverages1 <- data.frame(#taylororder = numeric(),inclusion_dist = numeric(),
  agebin = as.character(), delta = numeric(),
  agespecific = numeric(), M_ave  =numeric(), 
  centralage = numeric(), K_ave  =numeric(), O_ave  =numeric(), 
  M_agespec_ave = numeric(),K_agespec_ave = numeric(), 
  O_agespec_ave = numeric(), 
  O_agespec_ave_se = numeric(),  
  M_agespec_ave_se = numeric(), K_agespec_ave_se = numeric(),
  M_se = numeric(),  k_se = numeric(), 
  K_pvalue =  numeric(),M_pvalue =  numeric(), O_pvalue =  numeric())

for (ii in seq_along(delta)){
  
  iter_averages <- data.frame(#taylororder = numeric(),inclusion_dist = numeric(),
    agebin = as.character(), delta = numeric(),
    agespecific = numeric(), M_ave  =numeric(), 
    centralage = numeric(), K_ave  =numeric(), O_ave  =numeric(), 
    M_agespec_ave = numeric(),K_agespec_ave = numeric(), 
    O_agespec_ave = numeric(), 
    O_agespec_ave_se = numeric(),  
    M_agespec_ave_se = numeric(), K_agespec_ave_se = numeric(),
    M_se = numeric(),  k_se = numeric())
  
  for(jj in seq_along(central_age)){
    
    delt <-  delta[ii]
    
    centralage <- central_age[jj]
    
    diff_inc_estimates <- estimates%>% 
      filter(age >= centralage - delt & age <= centralage + delt) %>% 
      dplyr::select(age, Mahiane_diff_est = Mahiane, 
                    Kassanjee_diff_est = Kassanjee, 
                    weighted_diff_est = Weighted_diff)
    
    age25_estimates <-  bootstrap_sd %>%
      filter(age == central_age[jj]) %>%
      group_by(age, index) %>%
      dplyr::select(Mahiane_diff,Kassanjee_diff, weighted_diff) %>% 
      group_by(age) %>% 
      summarise(M_ave = mean(Mahiane_diff, na.rm = T),
                K_ave = mean(Kassanjee_diff, na.rm = T),
                O_ave = mean(weighted_diff, na.rm = T),
                M_avese = sd(Mahiane_diff, na.rm = T),
                K_avese = sd(Kassanjee_diff, na.rm = T),
                O_avese = sd(weighted_diff, na.rm = T))
    
    
    population_est <- populationest %>% 
      dplyr::select(age, year, females) %>% 
      filter(year == 2013) %>% 
      transmute(age = as.numeric(as.character(age)),
                population = females)%>% 
      filter(age >= centralage - delt & age <= centralage + delt)
    
    total_pop <- sum(population_est$population)
    
    bootstrap_pivoted <-  bootstrap_sd  %>% 
      filter(age >= centralage - delt & age <= centralage + delt) %>% 
      group_by(age,  index) %>%
      dplyr::select(Mahiane_diff,Kassanjee_diff, weighted_diff) %>% 
      inner_join(population_est) 
    
    inc_diff <- bootstrap_pivoted %>%
      inner_join(diff_inc_estimates) %>% 
      group_by(index) %>%
      transmute(age, index, Mahianeinc_diff = (Mahiane_diff*population)/total_pop,
                Kassanjeeinc_diff = (Kassanjee_diff*population)/total_pop,
                opt_diff = (weighted_diff*population)/total_pop,
                Mahiane_diff_est = (Mahiane_diff_est*population)/total_pop,
                Kassanjee_diff_est = (Kassanjee_diff_est*population)/total_pop,
                opt_diff_est = (weighted_diff_est*population)/total_pop)
    
    
    Mahianeeinc_diff_estimates <-  inc_diff %>%
      group_by(index) %>%
      summarise(Mave = sum(Mahiane_diff_est, na.rm = T),
                Kave = sum(Kassanjee_diff_est, na.rm = T),
                Oave = sum(opt_diff_est, na.rm = T),
                M_aveind = sum(Mahianeinc_diff, na.rm = T),
                K_aveind = sum(Kassanjeeinc_diff, na.rm = T),
                O_aveind = sum(opt_diff, na.rm = T)) 
    
    
    iteration_averages <- data.frame(age = centralage,
                                     delta = delta[ii],
                                     M_ave = mean(Mahianeeinc_diff_estimates$Mave),
                                     K_ave = mean(Mahianeeinc_diff_estimates$Kave),
                                     O_ave = mean(Mahianeeinc_diff_estimates$Oave),
                                     M_avese = sd(Mahianeeinc_diff_estimates$M_aveind, na.rm = T),
                                     K_avese = sd(Mahianeeinc_diff_estimates$K_aveind, na.rm = T),
                                     O_avese = sd(Mahianeeinc_diff_estimates$O_aveind, na.rm = T)) %>% 
      mutate(M_pvalue = 2*pnorm(-abs(M_ave/M_avese)), 
             K_pvalue = 2*pnorm(-abs(K_ave/K_avese)), 
             O_pvalue = 2*pnorm(-abs(O_ave/O_avese)) )
    
    iter_averages <- rbind(iteration_averages, iter_averages )
  }
  
  iterationaverages1 <- rbind(iterationaverages1, iter_averages)
}




pvalues <- iterationaverages1 %>% 
  dplyr::select(age, delta, 
                Mahiane = M_pvalue, Kassanjee = K_pvalue, 
                weighted = O_pvalue)%>%
  melt(id.vars = c("age", "delta")) %>% 
  dplyr::select(age, delta, variable,  value)




ggplot(pvalues)+
  geom_line(aes(x = age, y = value , colour= variable)) +
  labs(x = "age", y = "p-values", color = "")+
  facet_wrap(~delta)+
  theme_bw(base_size = 32, base_family = "") +
  geom_hline(aes(yintercept = 0.05), color = "black", linetype = "dashed")+
  scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
  scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,30), breaks = seq(16,30, 2))+
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2))+
  theme(legend.position = "bottom")

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/maleplots/pvalues.png", w=15,h=12)

