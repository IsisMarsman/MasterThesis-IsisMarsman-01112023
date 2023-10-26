library(ggplot2)
library(Matrix)

set.seed(091020235)

#MODEL 1
#Fixed 2

M <-10000
population_generations <- 10 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=5)
p_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_detection <- 0.5
prob_transmission <-0.5

r0_list_p_fixed_m1<- c()
r0_tracing_list_p_fixed_m1<- c()

N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)
generation_per_node<- list_generation_per_node_const2(population_generations)

contacts_transmission_m1 <- contacts_og


for(prob_tracing in t_interval){
  for(prob_clustering in p_interval){  
    infection_gens <- c()
  
    for(j in 1:M){
      
      contacts_tracing_m1 <- creating_contacts_constant2(contacts_og, contacts_extra, prob_clustering)
      matrix_final_m1 <- Matrix::sparseMatrix(i = contacts_tracing_m1[,1] , j = contacts_tracing_m1[,2], dims = c(N,N), symmetric = TRUE)
      
      gen_infected_nodes <- sim_epidemic_tracing_m1(N, contacts_tracing_m1, prob_tracing, prob_detection, generation_per_node, contacts_transmission_m1)
      infection_gens <- c(infection_gens, gen_infected_nodes)
      
    }
    
    
    no_infected_tracing <- table(infection_gens)
    r0_tracing <- no_infected_tracing[epidemic_generations]/no_infected_tracing[epidemic_generations-1]
    
    r0_tracing_list_p_fixed_m1 <- c(r0_tracing_list_p_fixed_m1,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
}


df_r0_p_fixed_m1 <- data.frame(x=p_interval, y= r0_tracing_list_p_fixed_m1, q = c(rep("t = 0",11),rep("t = 0.25",11), rep("t = 0.5",11), rep("t = 0.75",11),rep("t = 1",11)))


ggplot(df_r0_p_fixed_m1, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("t = 0" ='#D60270', "t = 0.25" = '#B92983',"t = 0.5" = '#9B4F96',"t = 0.75" = '#4E449F', "t = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.84))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")





#Pois(2)

M <-10000
population_generations <- 10 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=5)
p_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_detection <- 0.5
prob_transmission <-0.5

lambda_m1 <- 2

r0_list_p_pois_m1<- c()
r0_tracing_list_p_pois_m1<- c()


for(prob_tracing in t_interval){
  for(prob_clustering in p_interval){ 
    infection_gens <- c()
    
    for(j in 1:M){
      
      offspring_distr_m1 <- creating_offspringdistr_pois_m1(lambda_m1, population_generations)
      N <- number_of_nodes(offspring_distr_m1$offspring_distr)
      
      while(N< 400){
        offspring_distr_m1 <- creating_offspringdistr_pois_m1(lambda_m1, population_generations)
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
    
    r0_tracing_list_p_pois_m1 <- c(r0_tracing_list_p_pois_m1,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
}


df_r0_p_pois_m1 <- data.frame(x=p_interval, y= r0_tracing_list_p_pois_m1, q = c(rep("t = 0",11),rep("t = 0.25",11), rep("t = 0.5",11), rep("t = 0.75",11),rep("t = 1",11)))



ggplot(df_r0_p_pois_m1, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("t = 0" ='#D60270', "t = 0.25" = '#B92983',"t = 0.5" = '#9B4F96',"t = 0.75" = '#4E449F', "t = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.84))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")








#NB(2,0.5)

M <-10000
population_generations <- 10 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=5)
p_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_detection <- 0.5
prob_transmission <-0.5

nb_size_m1 <- 2
nb_prob_m1 <- 0.5


r0_list_p_nb_m1<- c()
r0_tracing_list_p_nb_m1<- c()


for(prob_tracing in t_interval){
  for(prob_clustering in p_interval){ 
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
    
    r0_tracing_list_p_nb_m1 <- c(r0_tracing_list_p_nb_m1,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
}


df_r0_p_nb_m1 <- data.frame(x=p_interval, y= r0_tracing_list_p_nb_m1, q = c(rep("t = 0",11),rep("t = 0.25",11), rep("t = 0.5",11), rep("t = 0.75",11),rep("t = 1",11)))


ggplot(df_r0_p_nb_m1, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("t = 0" ='#D60270', "t = 0.25" = '#B92983',"t = 0.5" = '#9B4F96',"t = 0.75" = '#4E449F', "t = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.84))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")





