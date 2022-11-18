############################################
# incidence difference weighted by the incidence difference
############################################
estimates <- estimate %>% 
  inner_join(weight) %>% 
  mutate(Weighted_diff  = (Mahiane *weight) +(1 - weight)*Kassanjee)

incidence_diff <-estimates %>% 
  dplyr::select(age, Mahiane, Kassanjee, Weighted = Weighted_diff) %>% 
  melt(id.vars = "age")
  

incidence_diffse <-estimates %>% 
  dplyr::select(age, Mahiane = Mahiane_diffse, 
                Kassanjee = Kassanjee_diffse, Weighted = weighted_diffse) %>% 
  melt(id.vars = "age") %>% 
  dplyr::select( age, variable, se = value)

incidence_diff <- inner_join(incidence_diff, incidence_diffse)

write.csv(incidence_diff, "/home/laurette/Desktop/Github/msfanalysis20211110/newoutputdata/males/incidence_diff_estse")

ggplot()+
  geom_ribbon(data = incidence_diff, aes(x = age, ymin = value - (1.96 * se),
                                         ymax = value + (1.96 * se), fill = variable), alpha = 0.3) +
  geom_line(data = incidence_diff, aes(x = age, y = value , colour= variable)) +
  labs(x = "age", y = "incidence difference ", color = "", fill = "")+
  theme_bw(base_size = 32, base_family = "") +
  geom_hline(aes(yintercept = 0), linetype ="dashed", color = "black")+
  scale_colour_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                      values = c("green4", "turquoise3", "red"))+
  scale_fill_manual(labels = c("Mahiane", "Kassanjee", "Optimal weighted"),
                    values = c("green4", "turquoise3", "red")) +
  scale_x_continuous(limits=c(15,34), breaks = seq(16,34, 2))+
  scale_y_continuous(limits=c(-0.08,0.1), breaks = seq(-0.06,0.1,0.02))+
  theme(legend.position = "bottom")


ggsave("/home/laurette/Dropbox/MhlangaWelteDeltaInc2020S2/MSFdataset/updatefigures20211020/maleplots/incdifference_stimates2.png", 
       w=15,h= 10)






