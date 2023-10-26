
library(ggplot2)
library(Matrix)


set.seed(021020232)


#Fixed 2

M <-10000
population_generations <- 10 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

prob_clustering <-0.5

r0_list_t<- c()
r0_tracing_list_t<- c()


N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)
generation_per_node<- list_generation_per_node_const2(population_generations)

contacts_transmission_m1 <- contacts_og


for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    infection_gens <- c()
    
    for(j in 1:M){
      
      contacts_tracing_m1 <- creating_contacts_constant2(contacts_og, contacts_extra, prob_clustering)
      matrix_final_m1 <- Matrix::sparseMatrix(i = contacts_tracing_m1[,1] , j = contacts_tracing_m1[,2], dims = c(N,N), symmetric = TRUE)
      
      gen_infected_nodes <- sim_epidemic_tracing_m1(N, contacts_tracing_m1, prob_tracing, prob_detection, generation_per_node, contacts_transmission_m1)
      infection_gens <- c(infection_gens, gen_infected_nodes)
      
    }
    
    
    no_infected_tracing <- table(infection_gens)
    r0_tracing <- no_infected_tracing[epidemic_generations]/no_infected_tracing[epidemic_generations-1]
    
    if(is.na(r0_tracing)){r0_tracing <- no_infected_tracing[4]/no_infected_tracing[3]
    }else{ 
      if(r0_tracing < 1.1){r0_tracing <- no_infected_tracing[4]/no_infected_tracing[3]}}
    
    r0_tracing_list_t <- c(r0_tracing_list_t,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}


r0_p0_list<-c()
for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    r0_p0 <- r0_tracing_list_t[1]*(1-prob_detection*prob_tracing)
    
    r0_p0_list <- c(r0_p0_list, r0_p0)
    
  }
}


r0_tracing_list_t[length(r0_tracing_list_t)] <- 0
df_r0_t <- data.frame(x=t_interval, y= r0_tracing_list_t, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0p0_t <- data.frame(x=t_interval, y= r0_p0_list, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0_t$dataset <- "p=0.5"
df_r0p0_t$dataset <- "no clustering"

combined_df <- rbind( df_r0_t, df_r0p0_t)

# Create the ggplot object
ggplot(data = combined_df, aes(x = x, y = y, color = q, linetype = dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(   "dashed", "solid")) +
  geom_point(data = df_r0_t, aes(x = x, y = y))+
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











#Pois(2)


M <-10000
population_generations <- 10 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

prob_clustering <-0.5
lambda<-2

r0_list_t_pois<- c()
r0_tracing_list_t_pois<- c()



for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    infection_gens <- c()
    
    for(j in 1:M){
      
      offspring_distr_m1 <- creating_offspringdistr_pois_m1(lambda, population_generations)
      N <- number_of_nodes(offspring_distr_m1$offspring_distr)
      
      while(N< 600){
        offspring_distr_m1 <- creating_offspringdistr_pois_m1(lambda, population_generations)
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
    
    r0_tracing_list_t_pois <- c(r0_tracing_list_t_pois,r0_tracing)
    cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}



r0_p0_list_pois<-c()
for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
    r0_p0 <- r0_tracing_list_t_pois[1]*(1-prob_detection*prob_tracing)
    
    r0_p0_list_pois <- c(r0_p0_list_pois, r0_p0)
    
  }
}


r0_tracing_list_t_pois[length(r0_tracing_list_t_pois)] <- 0
df_r0_t_pois <- data.frame(x=t_interval, y= r0_tracing_list_t_pois, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0p0_t_pois <- data.frame(x=t_interval, y= r0_p0_list_pois, q = c(rep("d = 0",11),rep("d = 0.25",11), rep("d = 0.5",11), rep("d = 0.75",11),rep("d = 1",11)))


df_r0_t_pois$dataset <- "p=0.5"
df_r0p0_t_pois$dataset <- "no clustering"


combined_df_pois <- rbind( df_r0_t_pois, df_r0p0_t_pois)

# Create the ggplot object
ggplot(data = combined_df_pois, aes(x = x, y = y, color = q, linetype = dataset)) +
  geom_line() +
  scale_linetype_manual(values = c(   "dashed", "solid")) +
  geom_point(data = df_r0_t_pois, aes(x = x, y = y))+
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


