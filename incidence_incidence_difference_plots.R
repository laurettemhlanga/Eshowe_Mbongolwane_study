library(ggplot2); library(dplyr);
library(tidyverse);library(reshape2)

sde_fem <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013kznonemodel_sd")
sde_mal <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013kznonemodel_sd")

sdefem <- sde_fem %>% 
  mutate(sex = "female")

sdemal <- sde_mal %>% 
  mutate(sex = "male")

weightdata <- bind_rows(sdefem, sdemal) %>% 
  group_by(age, dates, sex) %>% 
  transmute(weight = weight)


######################################
#  incidence estimates 
######################################

incestimatef <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013kznonemodel_est")
incestimatem <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013kznonemodel_est")

# incestimate <- analyse_separate_recmodels(msfdata_summary = msfdata)

incestimate_f <- incestimatef %>% 
  mutate(sex = "female")

incestimate_m <- incestimatem %>% 
  mutate(sex = "male")

inc_estimate <- bind_rows(incestimate_f, incestimate_m) %>% 
  inner_join(weightdata) %>% 
  mutate(incidenceo = (incidencem*weight) + 
           ((1 - weight)*incidencek )) %>% 
  select(age, dates, sex,incidencem,  incidencek, incidenceo)

write.csv(inc_estimate, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/bothsexes20182013kznonemodel_est")


incidence_long <- melt(inc_estimate, id.vars = c("dates","age", "sex")) 


sde_mod <- bind_rows(sdefem, sdemal)[,-c(1,4,5, 8,9)]
names(sde_mod) <- c( "age", "dates", "incidencem", "incidencek","incidenceo", "sex")
sde1 <- melt(sde_mod, id.vars = c("dates","age", "sex")) %>% 
  select(dates, age, sex,variable, se = value )

incidencelong <- inner_join(x = incidence_long, y = sde1)




newdata <- incidencelong %>% 
  filter(dates != 2015.5 & age <= 30)

ggplot()+
  geom_ribbon(data = newdata, aes(x = age, ymin = ifelse(value - (1.96 * se)<0, 0, value - (1.96 * se)),
                                  ymax = ifelse(value + (1.96 * se)>0.1, 0.1, value + (1.96 * se)), fill = variable), alpha = 0.3) +
  geom_line(data = newdata, aes(x = age, y = value , colour= variable)) +
  labs(x = "age", y = "estimate ", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  facet_grid(cols = vars(sex), rows = vars(dates))+
  scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
  scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,30), breaks = seq(16,30, 2))+
  scale_y_continuous(limits=c(0,0.1), breaks = seq(0,0.1,0.02))+
  theme(legend.position = "bottom")

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/plots20211119/incidence_shared_prev.png", w = 18, h = 12)



############################################
#difference
############################################


bootstrap_incidence_f <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/20182013bootsonemodel_sd")
bootstrap_incidence_m <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/20182013bootsonemodel_sd")


bootstrap_incidencef <- bootstrap_incidence_f %>% 
  mutate(sex = "female")

bootstrap_incidencem <- bootstrap_incidence_m %>% 
  mutate(sex = "male")

bootstrap_incidence <- bind_rows(bootstrap_incidencef, bootstrap_incidencem)

bootstrap_estimates<- bootstrap_incidence %>%
  group_by(age, dates, sex) %>% 
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
  group_by(age, dates,  sex,index) %>%
  select( incidencem, incidencek, incidenceo) %>%
  pivot_wider(names_from = c(dates), values_from = c(incidencem, incidencek, incidenceo)) 


incidence_diff <- bootstrap_estimates %>%
  mutate(Mahiane = incidencem_2018 - incidencem_2013, 
         Kassanjee = incidencek_2018 - incidencek_2013,
         Weighted = incidenceo_2018 - incidenceo_2013) %>% 
  group_by(age, sex) %>% 
  summarise(mse = sd(Mahiane, na.rm = T),
            kse = sd(Kassanjee, na.rm = T),
            wse = sd(Weighted, na.rm = T))

names(incidence_diff) <- c("age", "sex","Mahiane", "Kassanjee", "Weighted")

diff_sde_long <- melt(incidence_diff, id.vars = c("age", "sex"))
names(diff_sde_long) <- c("age","sex","variable", "se")

diff_estimates <- inc_estimate %>%
  select(age,  dates, sex, incidencem,  incidencek, incidenceo) %>% 
  filter(dates != 2015.5)%>%
  group_by(age, dates, sex) %>%
  select( incidencem, incidencek, incidenceo) %>%
  pivot_wider(names_from = c(dates), values_from = c(incidencem, incidencek, incidenceo)) 


incidence_diff <- diff_estimates %>%
  transmute(Mahiane = incidencem_2018 - incidencem_2013, 
            Kassanjee = incidencek_2018 - incidencek_2013,
            Weighted = incidenceo_2018 - incidenceo_2013) %>% 
  melt(id.vars = c("age", "sex")) %>% 
  inner_join(diff_sde_long)

 write.csv(incidence_diff, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/bothsexesincidence_diff_est")

 dummy2 = data.frame(x = c("female", "male"),
                     y = c(0,0))

ggplot()+
  geom_ribbon(data = incidence_diff, aes(x = age, ymin = ifelse(value - (1.96 * se)< -0.1, -0.1, 
                                                                value - (1.96 * se)),
                                         ymax = ifelse(value + (1.96 * se)> 0.1, 0.1, 
                                                       value + (1.96 * se)), 
                                         fill = variable), alpha = 0.3) +
  geom_line(data = incidence_diff, aes(x = age, y = value , colour= variable)) +
  labs(x = "age", y = "incidence difference ", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  geom_hline(data = dummy2, aes(yintercept = y), linetype ="dashed", color = "black")+
  facet_grid(cols = vars(sex))+
  scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
  scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,30), breaks = seq(16,30, 2))+
  scale_y_continuous(limits=c(-0.07,0.08), breaks = seq(-0.06,0.08,0.02))+
  theme(legend.position = "bottom")


ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/plots20211119/bothsexesincdifference_stimates2.png", 
       w=20,h= 9)


######################
#prevalence
#####################
prev_se <- bind_rows(sdefem, sdemal) %>% 
  select(age,  dates,sex,prevse)

prev_estf <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/prev_est") %>% 
  mutate(sex = "female")


prev_estm<- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/prev_est") %>% 
  mutate(sex = "male")

prev_data <- bind_rows(prev_estf, prev_estm) 


prevdata <- inner_join(prev_data, prev_se) %>% 
  filter(dates != 2015.5)


ggplot()+
  geom_ribbon(data = prevdata, aes(x = age, ymin = ifelse(prevalence - (1.96 * prevse)< 0, 0, 
                                                          prevalence - (1.96 * prevse)),
                                         ymax = prevalence + (1.96 * prevse),
                                   fill  = as.factor(dates)), alpha = 0.3) +
  geom_line(data = prevdata, aes(x = age, y = prevalence, color = as.factor(dates))) +
  labs(x = "age", y = "estimated prevalence", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  facet_grid(cols = vars(sex))+
  scale_x_continuous(limits=c(15,35), breaks = seq(16,35, 2))+
  scale_y_continuous(limits=c(0,0.65), breaks = seq(0,0.6,0.1))+
  theme(legend.position = "bottom")

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/plots20211119/prevalence_estimates_bothdatsets_coloured.png", 
       w=15,h= 10)
