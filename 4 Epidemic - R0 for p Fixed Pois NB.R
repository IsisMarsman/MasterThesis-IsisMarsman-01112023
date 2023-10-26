

library(ggplot2)
library(Matrix)


set.seed(111020234)


#Fixed 2

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

p_interval <- seq(from=0, to=1, length.out=11)
q_interval <- seq(from=0.25, to=1, length.out=4)

epidemic_generations <- floor(population_generations / 2) +1



r0_list_p<- c()


N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)


for(prob_transmission in q_interval){ 
  for(prob_clustering in p_interval){
    sumnoinf <- rep(0, epidemic_generations)
    for(j in 1:M){
      
      contacts_m2 <- creating_contacts_constant2(contacts_og, contacts_extra, prob_clustering)
      matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission)
      epidemic <- sim_epidemic(N, epidemic_generations, matrix_final_m2)
      sumnoinf<- sumnoinf + epidemic$n_infected  
      
    }
    
    
    
    r0 <- sumnoinf[epidemic_generations]/sumnoinf[epidemic_generations-1]
    if(is.na(r0)){r0 <- sumnoinf[5]/sumnoinf[4]
    }else{ 
      if(r0 < 1.1){r0 <- sumnoinf[5]/sumnoinf[4]}}
    
    r0_list_p <- c(r0_list_p,r0)
    cat("p = ", prob_clustering, "q=", prob_transmission, "R0 =", r0, "\n")
  }
  
}



r0_list_p <- c(rep(0,11), r0_list_p)

df_r0_p <- data.frame(x=p_interval, y= r0_list_p, q = c(rep("q = 0",11),rep("q = 0.25",11), rep("q = 0.5",11), rep("q = 0.75",11),rep("q = 1",11)))

ggplot(df_r0_p, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("q = 0" ='#D60270', "q = 0.25" = '#B92983',"q = 0.5" = '#9B4F96',"q = 0.75" = '#4E449F', "q = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.1, 0.81))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,4)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")




#Pois 2

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

lambda <-2

p_interval <- seq(from=0, to=1, length.out=11)
q_interval <- seq(from=0.25, to=1, length.out=4)

epidemic_generations <- floor(population_generations / 2) +1



r0_list_p_pois<- c()


for(prob_transmission in q_interval){ 
  for(prob_clustering in p_interval){
    sumnoinf <- rep(0, epidemic_generations)
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
      
    }
    
    
    
    r0 <- sumnoinf[epidemic_generations]/sumnoinf[epidemic_generations-1]
    if(is.na(r0)){r0 <- sumnoinf[5]/sumnoinf[4]
    }else{ 
      if(r0 < 1.1){r0 <- sumnoinf[5]/sumnoinf[4]}}
    
    r0_list_p_pois <- c(r0_list_p_pois,r0)
    cat("p = ", prob_clustering, "q=", prob_transmission, "R0 =", r0, "\n")
  }
  
}

r0_list_p_pois <- c(rep(0,11), r0_list_p_pois)

df_r0_p_pois <- data.frame(x=p_interval, y= r0_list_p_pois, q = c(rep("q = 0",11),rep("q = 0.25",11), rep("q = 0.5",11), rep("q = 0.75",11),rep("q = 1",11)))

ggplot(df_r0_p_pois, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("q = 0" ='#D60270', "q = 0.25" = '#B92983',"q = 0.5" = '#9B4F96',"q = 0.75" = '#4E449F', "q = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,4.01)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")



#NB(2,0.5)


M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

nb_size <-2
nb_prob <- 0.5

p_interval <- seq(from=0, to=1, length.out=11)
q_interval <- seq(from=0.25, to=1, length.out=4)

epidemic_generations <- floor(population_generations / 2) +1



r0_list_p_nb<- c()


for(prob_transmission in q_interval){ 
  for(prob_clustering in p_interval){
    sumnoinf <- rep(0, epidemic_generations)
    for(j in 1:M){
      
      offspring_distr_m2 <- creating_offspringdistr_nb(nb_size, nb_prob, population_generations)
      N <- number_of_nodes(offspring_distr_m2)
      while(N< 100){
        offspring_distr_m2 <- creating_offspringdistr_nb(nb_size, nb_prob, population_generations)
        N <- number_of_nodes(offspring_distr_m2)
      }
      contacts_m2 <- creating_contacts_DK(offspring_distr_m2, N, prob_clustering)
      matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission)
      
      
      epidemic <- sim_epidemic(N, epidemic_generations, matrix_final_m2)
      sumnoinf<- sumnoinf + epidemic$n_infected  
      
    }
    
    
    
    r0 <- sumnoinf[epidemic_generations]/sumnoinf[epidemic_generations-1]
    if(is.na(r0)){r0 <- sumnoinf[5]/sumnoinf[4]
    }else{ 
      if(r0 < 1.1){r0 <- sumnoinf[5]/sumnoinf[4]}}
    
    r0_list_p_nb <- c(r0_list_p_nb,r0)
    cat("p = ", prob_clustering, "q=", prob_transmission, "R0 =", r0, "\n")
  }
  
}


r0_list_p_nb <- c(rep(0,11), r0_list_p_nb)

df_r0_p_nb <- data.frame(x=p_interval, y= r0_list_p_nb, q = c(rep("q = 0",11),rep("q = 0.25",11), rep("q = 0.5",11), rep("q = 0.75",11),rep("q = 1",11)))

ggplot(df_r0_p_nb, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("q = 0" ='#D60270', "q = 0.25" = '#B92983',"q = 0.5" = '#9B4F96',"q = 0.75" = '#4E449F', "q = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,4)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")



