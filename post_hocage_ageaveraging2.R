############################
###########estimates
############################

# bootstrap_incidence_f <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013bootsonemodel_sd")
bootstrap_incidence <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013bootsonemodel_sd")


# bootstrap_incidencef <- bootstrap_incidence_f %>% 
#   mutate(sex = "female")
# 
# bootstrap_incidencem <- bootstrap_incidence_m %>% 
#   mutate(sex = "male")
# 
# bootstrap_incidence <- bind_rows(bootstrap_incidencef, bootstrap_incidencem)
# 
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

# estimates <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/incidence_diff_estse")

#############################
#Post hoc averaging into 5 year age bins
#############################
binslist <-  list(c(15,19), c(20, 24), c(25,29), 
                  c(30, 34), c(35, 39))
binlabels <- c( "15-19","20-24",
                "25-29", "30-34", 
                "35-39")
central_age <- seq(17,40,5)


populationest <- read.csv("/home/laurette/Desktop/Github/MSF_data/poulation_est.csv") %>% 
  filter(age != "90+") %>% 
  mutate(age = as.numeric(as.character(age)))


iterationaverages <- data.frame(#taylororder = numeric(),inclusion_dist = numeric(),
  agebin = as.character(), delta = numeric(),
  agespecific = numeric(), M_ave  =numeric(), 
  centralage = numeric(), K_ave  =numeric(), O_ave  =numeric(), 
  M_agespec_ave = numeric(),K_agespec_ave = numeric(), 
  O_agespec_ave = numeric(), 
  O_agespec_ave_se = numeric(),  
  M_agespec_ave_se = numeric(), K_agespec_ave_se = numeric(),
  M_se = numeric(),  k_se = numeric())


for(jj in seq_along(binslist)){
  
  
  lowerage <- binslist[[jj]][1]
  upperage <- binslist[[jj]][2]
  
  
  diff_inc_estimates <- estimates%>% 
    filter((age >= lowerage & age <= upperage)) %>% 
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
    dplyr::select(age, dates = year, males) %>% 
    filter(dates == 2013) %>% 
    transmute(age = age,
              population = males) %>% 
    filter(age >= lowerage & age <= upperage)
  
  
  total_pop <- sum(population_est$population)
  
  bootstrap_pivoted <-  bootstrap_sd  %>% 
    filter(age >= lowerage & age <= upperage) %>%  
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
              O_aveind = sum(opt_diff, na.rm = T)) %>%
    summarise(M_ave = mean(Mave, na.rm = T),
              K_ave = mean(Kave, na.rm = T),
              O_ave = mean(Oave, na.rm = T),
              M_avese = sd(M_aveind, na.rm = T),
              K_avese = sd(K_aveind, na.rm = T),
              O_avese = sd(O_aveind, na.rm = T))
  
  
  
  iteration_averages <- data.frame(#taylororder = prevtaylororder,
    # inclusion_dist = previnclusion_dist,
    centralage = age25_estimates$age,
    agebin = binlabels[jj],
    M_ave = Mahianeeinc_diff_estimates$M_ave,
    K_ave = Mahianeeinc_diff_estimates$K_ave,
    O_ave = Mahianeeinc_diff_estimates$O_ave,
    M_se = Mahianeeinc_diff_estimates$M_avese,
    K_se = Mahianeeinc_diff_estimates$K_avese,
    O_se = Mahianeeinc_diff_estimates$O_avese,
    M_agespec_ave = age25_estimates$M_ave,
    K_agespec_ave = age25_estimates$K_ave,
    O_agespec_ave =age25_estimates$O_ave,
    M_agespec_ave_se = age25_estimates$M_avese, 
    K_agespec_ave_se = age25_estimates$K_avese,
    O_agespec_ave_se = age25_estimates$O_avese)
  
  iterationaverages <- rbind(iterationaverages, iteration_averages)
}

 iteration_averages <- iterationaverages %>% 
   mutate(sex = "male")




# write.csv(iterationaverages, dataname)
# incidence difference averaging results 
# iterationaverages <- read.csv(dataname)


library(RColorBrewer)

newdatase <- iteration_averages %>%
  dplyr::select(agebin, sex, M_se, K_se, O_se,
                M_agespec_ave_se, K_agespec_ave_se,O_agespec_ave_se) %>% 
  transmute(agebin = agebin,
            sex = sex,
            M_ave = M_se,
            K_ave = K_se,
            O_ave = O_se,
            M_agespec_ave = M_agespec_ave_se,
            K_agespec_ave = K_agespec_ave_se,
            O_agespec_ave = O_agespec_ave_se) %>% 
  melt(id.vars = c("agebin", "sex"))
names(newdatase) <- c( "agebin", "sex","variable","se")


newdata <- iteration_averages %>%
  dplyr::select(agebin, sex, M_ave, M_agespec_ave,
                K_ave, K_agespec_ave, O_ave,
                O_agespec_ave) %>% 
  melt(id.vars = c("agebin", "sex"))  %>% 
  inner_join(newdatase) 

write.csv(newdata, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_5yrbins_estse")


  ggplot()+
  geom_point(data = newdata, aes(x = agebin, y = value, color = variable), size = 3, position = position_dodge(0.5))+
  geom_hline( aes(yintercept = 0), linetype = "dashed", color = "2")+
  geom_errorbar(data = newdata, aes(x = agebin, ymin = value - (qnorm(1 - 0.05/2) * se), 
                    ymax = value + (qnorm(1 - 0.05/2) * se), color = variable), position = position_dodge(0.5))+
  labs(x = "age bin", y = "incidence difference", color = "")+
  theme_bw(base_size = 24, base_family = "") +
  theme(axis.text.x = element_text())+
  scale_y_continuous(limits = c(-0.08, 0.08), breaks = seq(-0.08, 0.08, 0.02))+
  # scale_x_continuous(limits = c(14.5,30), breaks = seq(15,30, 0.01))+
   # scale_colour_manual(labels = c("Mahiane age range",
   #                               "Mahiane age specific",
   #                               "Kassanjee age range",
   #                               "Kassanjee age specific",
   #                               "Weighted age range",
   #                               "Weighted age specific"),
                      # values = brewer.pal(n = 6, name = "Dark2"))+
  theme( legend.position = "bottom", axis.text.x = element_text())

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/maleplots/incdifference_stimates_5agebin2.png", 
       w=15,h= 10)

# ggsave(paste0("agebintotalpopulation1000",".png"), w=7,h=4.5, dpi=2400)

# iterationaverages %>% 
#   dplyr::select(agebin,centralage ,M_pop_pvalue, K_pop_pvalue,
#                 M_agespec_pvalue, K_agespec_pvalue)
# 


###################################################################################################################
# delta t from one central age
###################################################################################################################

delta = 5; central_age <- 25



iterationaverages1 <- data.frame(#taylororder = numeric(),inclusion_dist = numeric(),
  agebin = as.character(), delta = numeric(),
  agespecific = numeric(), M_ave  =numeric(), 
  centralage = numeric(), K_ave  =numeric(), O_ave  =numeric(), 
  M_agespec_ave = numeric(),K_agespec_ave = numeric(), 
  O_agespec_ave = numeric(), 
  O_agespec_ave_se = numeric(),  
  M_agespec_ave_se = numeric(), K_agespec_ave_se = numeric(),
  M_se = numeric(),  k_se = numeric())

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
      dplyr::select(age, year,males) %>% 
      filter(year == 2013) %>% 
      transmute(age = as.numeric(as.character(age)),
                population = males)%>% 
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
                                     O_avese = sd(Mahianeeinc_diff_estimates$O_aveind, na.rm = T))
    
    iter_averages <- rbind(iteration_averages, iter_averages )
  }
  
  iterationaverages1 <- rbind(iterationaverages1, iter_averages)
}




sde <- iterationaverages1 %>% 
  dplyr::select(age, delta, 
                Mahiane = M_avese, Kassanjee = K_avese, 
                weighted = O_avese)%>%
  melt(id.vars = c("age", "delta")) %>% 
  dplyr::select(age, delta, variable, se = value)
  

newdata <- iterationaverages1 %>%
  dplyr::select(age, delta, Mahiane = M_ave, 
                Kassanjee = K_ave , 
                weighted = O_ave )%>%
  melt(id.vars = c("age", "delta")) %>% 
  inner_join(sde)

  ggplot(newdata)+
    geom_ribbon(aes(x = age, ymin = value - (1.96 * se),
                                         ymax = value + (1.96 * se), fill = variable), alpha = 0.3) +
    geom_line(aes(x = age, y = value , colour= variable)) +
    labs(x = "age", y = "incidence difference ", color = "", fill = "")+
    facet_wrap(~delta)+
    theme_bw(base_size = 32, base_family = "") +
    geom_hline(aes(yintercept = 0), linetype ="dashed", color = "black")+
    scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
    scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,30), breaks = seq(16,30, 2))+
  scale_y_continuous(limits=c(-0.07,0.08), breaks = seq(-0.06,0.08,0.02))+
  theme(legend.position = "bottom")

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/maleplots/facted_incdiff2.png", w=15,h=12)

