
library(ggplot2)
library(Matrix)


set.seed(15092023)


#Fixed 2

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

p_interval <- seq(from=0, to=1, length.out=5)
q_interval <- seq(from=0.2, to=1, length.out=9)

epidemic_generations <- floor(population_generations / 2) +1



r0_list_q<- c()


N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)

for(prob_clustering in p_interval){
  for(prob_transmission in q_interval){ 
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
    
    
    r0_list_q <- c(r0_list_q, r0)
    cat("p = ", prob_clustering, "q=", prob_transmission, "R0 =", r0, "\n")
  }
  
}



r0_list_q_plot <- c(0 , r0_list_q[1:9], 0, r0_list_q[10:18], 0 ,r0_list_q[19:27], 0 ,r0_list_q[28:36], 0 ,r0_list_q[37:45])

q_interval_plot <- c(0, q_interval)

df_r0_q_plot <- data.frame(x=q_interval_plot, y= r0_list_q_plot, q = c(rep("p = 0",10),rep("p = 0.25",10), rep("p = 0.5",10), rep("p = 0.75",10),rep("p = 1",10)))

ggplot(df_r0_q_plot, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("p = 0" ='#D60270', "p = 0.25" = '#B92983',"p = 0.5" = '#9B4F96',"p = 0.75" = '#4E449F', "p = 1" =  '#0038A8'))+ 
  theme_light()+
  theme(legend.position = c(0.11, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,4.01)+
  labs(x = "Transmission Probability (q)", y = "Reproduction Number")




#Pois(2)

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

lambda <- 2

p_interval <- seq(from=0, to=1, length.out=5)
q_interval <- seq(from=0.2, to=1, length.out=9)

epidemic_generations <- floor(population_generations / 2) +1



r0_list_q<- c()


for(prob_clustering in p_interval){
  for(prob_transmission in q_interval){ 
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
    
    
    r0_list_q <- c(r0_list_q, r0)
    cat("p = ", prob_clustering, "q=", prob_transmission, "R0 =", r0, "\n")
  }
  
}


r0_list_q_pois_plot <- c(0 , r0_list_q[1:9], 0, r0_list_q[10:18], 0 ,r0_list_q[19:27], 0 ,r0_list_q[28:36], 0 ,r0_list_q[37:45])

q_interval_plot <- c(0, q_interval)

df_r0_q_pois_plot <- data.frame(x=q_interval_plot, y= r0_list_q_pois_plot, q = c(rep("p = 0",10),rep("p = 0.25",10), rep("p = 0.5",10), rep("p = 0.75",10),rep("p = 1",10)))

ggplot(df_r0_q_pois_plot, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("p = 0" ='#D60270', "p = 0.25" = '#B92983',"p = 0.5" = '#9B4F96',"p = 0.75" = '#4E449F', "p = 1" =  '#0038A8'))+ 
  theme_light()+
  theme(legend.position = c(0.15, 0.75))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,4.01)+
  labs(x = "Transmission Probability (q)", y = "Reproduction Number")





#NB(2,0.5)



M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

nb_size <-2
nb_prob <- 0.5

p_interval <- seq(from=0, to=1, length.out=5)
q_interval <- seq(from=0.2, to=1, length.out=9)

epidemic_generations <- floor(population_generations / 2) +1



r0_list_q<- c()


for(prob_clustering in p_interval){
  for(prob_transmission in q_interval){ 
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
    
    
    r0_list_q <- c(r0_list_q, r0)
    cat("p = ", prob_clustering, "q=", prob_transmission, "R0 =", r0, "\n")
  }
  
}


r0_list_q_nb_plot <- c(0 , r0_list_q[1:9], 0, r0_list_q[10:18], 0 ,r0_list_q[19:27], 0 ,r0_list_q[28:36], 0 ,r0_list_q[37:45])

q_interval_plot <- c(0, q_interval)

df_r0_q_nb_plot <- data.frame(x=q_interval_plot, y= r0_list_q_nb_plot, q = c(rep("p = 0",10),rep("p = 0.25",10), rep("p = 0.5",10), rep("p = 0.75",10),rep("p = 1",10)))

ggplot(df_r0_q_nb_plot, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("p = 0" ='#D60270', "p = 0.25" = '#B92983',"p = 0.5" = '#9B4F96',"p = 0.75" = '#4E449F', "p = 1" =  '#0038A8'))+ 
  theme_light()+
  theme(legend.position = c(0.11, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,4.05)+
  labs(x = "Transmission Probability (q)", y = "Reproduction Number")







