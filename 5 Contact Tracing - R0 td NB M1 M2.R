
library(ggplot2)
library(Matrix)

set.seed(041020232)

#NB(2,0.5)

#model 2 

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

prob_transmission <- 0.5
prob_clustering <-0.5

nb_size_m2 <- 2
nb_prob_m2 <- 0.5

r0_list_t_nb<- c()
r0_tracing_list_t_nb<- c()
sumnoinf <- rep(0, epidemic_generations)
sumnoinf_tracing <- rep(0, epidemic_generations)


for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    sumnoinf_tracing <- rep(0, epidemic_generations)
    
    for(j in 1:M){
      
      offspring_distr_m2 <- creating_offspringdistr_nb(nb_size_m2, nb_prob_m2, population_generations)
      N <- number_of_nodes(offspring_distr_m2)
      while(N< 1000){
        offspring_distr_m2 <- creating_offspringdistr_nb(nb_size_m2, nb_prob_m2, population_generations)
        N <- number_of_nodes(offspring_distr_m2)
      }
      contacts_m2 <- creating_contacts_DK(offspring_distr_m2, N, prob_clustering)
      matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission)
      
      epidemic_tracing <- sim_epidemic_tracing_m2_v2(N, matrix_final_m2, contacts_m2, prob_tracing, prob_detection)
      sumnoinf_tracing <- sumnoinf_tracing + epidemic_tracing
      
      
    }
    
    
    
    r0_tracing <- sumnoinf_tracing[epidemic_generations]/sumnoinf_tracing[epidemic_generations-1]
    
    if(is.na(r0_tracing)){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]
    }else{ 
      if(r0_tracing < 1.1){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]}}
    
    r0_tracing_list_t_nb <- c(r0_tracing_list_t_nb,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}




prob_detection<-1

  for(prob_tracing in c(0.5,0.6,0.7,0.8,0.9,1)){
    sumnoinf_tracing <- rep(0, epidemic_generations)
    
    for(j in 1:M){
      
      offspring_distr_m2 <- creating_offspringdistr_nb(nb_size_m2, nb_prob_m2, population_generations)
      N <- number_of_nodes(offspring_distr_m2)
      while(N< 1000){
        offspring_distr_m2 <- creating_offspringdistr_nb(nb_size_m2, nb_prob_m2, population_generations)
        N <- number_of_nodes(offspring_distr_m2)
      }
      contacts_m2 <- creating_contacts_DK(offspring_distr_m2, N, prob_clustering)
      matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission)
      
      epidemic_tracing <- sim_epidemic_tracing_m2_v2(N, matrix_final_m2, contacts_m2, prob_tracing, prob_detection)
      sumnoinf_tracing <- sumnoinf_tracing + epidemic_tracing
      
      
    }
    
    
    
    r0_tracing <- sumnoinf_tracing[epidemic_generations]/sumnoinf_tracing[epidemic_generations-1]
    
    if(is.na(r0_tracing)){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]
    }else{ 
      if(r0_tracing < 1.1){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]}}
    
    r0_tracing_list_t_nb <- c(r0_tracing_list_t_nb,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  



r0_p0_list_nb<-c()
for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    r0_p0 <- r0_tracing_list_t_nb[1]*(1-prob_detection*prob_tracing)
    
    r0_p0_list_nb <- c(r0_p0_list_nb, r0_p0)
    
  }
}


r0_tracing_list_t_nb[length(r0_tracing_list_t_nb)] <- 0
df_r0_t_nb <- data.frame(x=t_interval, y= r0_tracing_list_t_nb, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0p0_t_nb <- data.frame(x=t_interval, y= r0_p0_list_nb, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0_t_nb$dataset <- "p=0.5 q=0.5"
df_r0p0_t_nb$dataset <- "no clustering"

combined_df_nb <- rbind( df_r0_t_nb, df_r0p0_t_nb)

# Create the ggplot object
ggplot(data = combined_df_nb, aes(x = x, y = y, color = q, linetype = dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(   "dashed", "solid")) +
  geom_point(data = df_r0_t_nb, aes(x = x, y = y))+
  scale_color_manual(values = c("d = 0" ='#D60270', "d = 0.25" = '#B92983',"d = 0.5" = '#9B4F96',"d = 0.75" = '#4E449F', "d = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.24))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,2.2)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Tracing Probability (t)", y = "Reproduction Number")






#model 1


M <-10000
population_generations <- 10 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

prob_clustering <-0.5

nb_size_m1 <- 2
nb_prob_m1 <- 0.5

r0_list_t_nb_m1<- c()
r0_tracing_list_t_nb_m1<- c()



for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    infection_gens <- c()
    
    for(j in 1:M){
      
      offspring_distr_m1 <- creating_offspringdistr_nb_m1(nb_size_m1, nb_prob_m1, population_generations)
      N <- number_of_nodes(offspring_distr_m1$offspring_distr)
      
      while(N< 400){
        offspring_distr_m1 <- creating_offspringdistr_nb_m1(nb_size_m1, nb_prob_m1, population_generations)
        N <- number_of_nodes(offspring_distr_m1$offspring_distr)
      }
      
      contacts_transmission_m1 <- create_transmission_network_m1(offspring_distr_m1$offspring_distr)
      
      contacts_tracing_m1 <- creating_contacts_DK(offspring_distr_m1$offspring_distr, N, prob_clustering)
      matrix_tracing_m1 <- Matrix::sparseMatrix(i = contacts_tracing_m1[,1] , j = contacts_tracing_m1[,2], x=1, dims = c(N,N), symmetric = TRUE) 
      
      generation_per_node<- list_generation_per_node_distr(population_generations, offspring_distr_m1$no_nodes)
      
      gen_infected_nodes <- sim_epidemic_tracing_m1(N, contacts_tracing_m1, prob_tracing, prob_detection, generation_per_node, contacts_transmission_m1)
      infection_gens <- c(infection_gens, gen_infected_nodes)
      
    }
    
    
    no_infected_tracing <- table(infection_gens)
    r0_tracing <- no_infected_tracing[epidemic_generations]/no_infected_tracing[epidemic_generations-1]
    
    if(is.na(r0_tracing)){r0_tracing <- no_infected_tracing[4]/no_infected_tracing[3]
    }else{ 
      if(r0_tracing < 1.1){r0_tracing <- no_infected_tracing[4]/no_infected_tracing[3]}}
    
    r0_tracing_list_t_nb_m1 <- c(r0_tracing_list_t_nb_m1,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}



r0_p0_list_nb_m1<-c()
for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    r0_p0 <- r0_tracing_list_t_nb_m1[1]*(1-prob_detection*prob_tracing)
    
    r0_p0_list_nb_m1 <- c(r0_p0_list_nb_m1, r0_p0)
    
  }
}


r0_tracing_list_t_nb_m1[length(r0_tracing_list_t_nb_m1)] <- 0
df_r0_t_nb_m1 <- data.frame(x=t_interval, y= r0_tracing_list_t_nb_m1, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0p0_t_nb_m1 <- data.frame(x=t_interval, y= r0_p0_list_nb_m1, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0_t_nb_m1$dataset <- "p=0.5"
df_r0p0_t_nb_m1$dataset <- "no clustering"

combined_df_nb_m1 <- rbind( df_r0_t_nb_m1, df_r0p0_t_nb_m1)

# Create the ggplot object
ggplot(data = combined_df_nb_m1, aes(x = x, y = y, color = q, linetype = dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(   "dashed", "solid")) +
  geom_point(data = df_r0_t_nb_m1, aes(x = x, y = y))+
  scale_color_manual(values = c("d = 0" ='#D60270', "d = 0.25" = '#B92983',"d = 0.5" = '#9B4F96',"d = 0.75" = '#4E449F', "d = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.24))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,2.2)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Tracing Probability (t)", y = "Reproduction Number")

