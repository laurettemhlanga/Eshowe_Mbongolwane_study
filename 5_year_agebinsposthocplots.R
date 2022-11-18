female_bin <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/females/incidence_diff_5yrbins_estse")
male_bin <- read.csv("/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_5yrbins_estse")

 newdata <- bind_rows(female_bin, male_bin) %>% 
   filter(agebin != "30-34" & agebin != "35-39")
   
  
 
 newdata %>%
   # filter(agebin != "35-39") %>% 
   ggplot( aes(x = agebin,   shape = variable, color = variable))+
   geom_point(aes(y = value), size = 5, position = position_dodge(0.75))+
   geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black")+
   facet_grid(cols = vars(sex))+
   geom_errorbar(aes(ymin = value - (qnorm(1 - 0.05/2) * se), 
                    ymax = value + (qnorm(1 - 0.05/2) * se)), 
                 position = position_dodge(0.75))+
   labs(x = "age bin", y = "incidence difference", color = "", shape = "")+
   theme(axis.text.x = element_text())+
   scale_shape_manual(labels = c("Kassanjee age specific", "Kassanjee age range",
                                 "Mahiane age specific", "Mahiane age range",
                                 "Optimal weighted age specific", "Optimal weighted age range"),
                      values = c (20, 18, 20, 18, 20, 18))+
   scale_colour_manual(labels = c("Kassanjee age specific", "Kassanjee age range",
                                  "Mahiane age specific", "Mahiane age range",
                                  "Optimal weighted age specific", "Optimal weighted age range"),
                      values = c("#DF536B", "#DF536B", "#2297E6", "#2297E6", "red", "red"))+
   theme_bw(base_size = 32, base_family = "")+
   scale_y_continuous(limits = c(-0.06, 0.035), breaks = round(seq(-0.06, 0.04, 0.02), digits = 3))+
   theme( legend.position = "bottom", axis.text.x = element_text())
 

ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/plots20211119/incdifference_stimates_5agebin2.png", 
       w=20,h= 10)
