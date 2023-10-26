
set.seed(091020238)

#Model2
M <-1000
population_generations <- 20 #generations in the tree, not counting index case as generation

prob_clustering <- 0.5              #probability to add the extra edges
prob_transmission <- 0.326

lambda <- 2


epidemic_generations <- floor(population_generations / 2) +1


sumnoinf <- rep(0, epidemic_generations)
r0_list<- c()


for(j in 1:M){
  
  offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
  N <- number_of_nodes(offspring_distr_m2)
  while(N< 100){
    offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
    N <- number_of_nodes(offspring_distr_m2)
  }
  contacts_m2 <- creating_contacts_DK(offspring_distr_m2, N, prob_clustering)
  matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission)
  epidemic <- sim_epidemic(N, epidemic_generations, matrix_final_m2)
  
  sumnoinf<- sumnoinf + epidemic$n_infected  
  cat(j)
}

#calculating R0
r0 <- sumnoinf[epidemic_generations]/sumnoinf[epidemic_generations-1]


r0_epidemic_m2 <- c(0)
for (i in 1:(epidemic_generations-1)) {
  r0_epidemic_m2[i+1] <- sumnoinf[i+1]/sumnoinf[i]
}

r0_epidemic_m2




df_r0 <- data.frame(x=0:(length(r0_epidemic_m2)-1), y= r0_epidemic_m2)

ggplot(df_r0, aes(x=x, y=y)) +
  geom_point(color ='#0038A8') +
  geom_line(color ='#0038A8') +
  theme_light()+
  theme(legend.position = c(0.125, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  theme(
    #panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.y = element_blank()
  )+
  ylim(0 ,2.5)+
  scale_x_continuous(breaks = seq(0, 12, by = 1))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")

