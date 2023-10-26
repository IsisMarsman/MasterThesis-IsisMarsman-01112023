library(Matrix)
library(ggplot2)

set.seed(071020231)

#model 2
prob_detection <- 0.5

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_transmission_m2 <- 0.201
prob_clustering_m2 <-1

lambda <- 2

r0_tracing_list_t_m2<- c()


for(prob_tracing in t_interval){
  sumnoinf_tracing <- rep(0, epidemic_generations)
  
  for(j in 1:M){
    offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
    N <- number_of_nodes(offspring_distr_m2)
    while(N< 1000){
      offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
      N <- number_of_nodes(offspring_distr_m2)
    }
    contacts_m2 <- creating_contacts_DK(offspring_distr_m2, N, prob_clustering_m2)
    matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission_m2)
    
    epidemic_tracing <- sim_epidemic_tracing_m2_v2(N, matrix_final_m2, contacts_m2, prob_tracing, prob_detection)
    sumnoinf_tracing <- sumnoinf_tracing + epidemic_tracing
    
    
  }
  
  
  
  r0_tracing <- sumnoinf_tracing[epidemic_generations]/sumnoinf_tracing[epidemic_generations-1]
  
  if(is.na(r0_tracing)){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]
  }else{ 
    if(r0_tracing < 1.1){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]}}
  
  r0_tracing_list_t_m2 <- c(r0_tracing_list_t_m2,r0_tracing)
  cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
}






#model 1 

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

#results from table hfst 7
prob_clustering_m1 <-  0.6319
nb_size_m1 <- (1.5061^2)/(2.1135^2 -1.5061)
nb_prob_m1 <- 1.5061/(2.1135^2)

r0_tracing_list_t_m1<- c()

for(prob_tracing in t_interval){
  infection_gens <- c()
  
  for(j in 1:M){
    offspring_distr_m1 <- creating_offspringdistr_nb_m1(nb_size_m1, nb_prob_m1, population_generations)
    N <- number_of_nodes(offspring_distr_m1$offspring_distr)
    
    while(N< 600){
      offspring_distr_m1 <- creating_offspringdistr_nb_m1(nb_size_m1, nb_prob_m1, population_generations)
      N <- number_of_nodes(offspring_distr_m1$offspring_distr)
    }
    
    contacts_transmission_m1 <- create_transmission_network_m1(offspring_distr_m1$offspring_distr)
    
    contacts_tracing_m1 <- creating_contacts_DK(offspring_distr_m1$offspring_distr, N, prob_clustering_m1)
    matrix_tracing_m1 <- Matrix::sparseMatrix(i = contacts_tracing_m1[,1] , j = contacts_tracing_m1[,2], x=1, dims = c(N,N), symmetric = TRUE) 
    
    generation_per_node<- list_generation_per_node_distr(population_generations, offspring_distr_m1$no_nodes)
    
    gen_infected_nodes <- sim_epidemic_tracing_m1(N, contacts_tracing_m1, prob_tracing, prob_detection, generation_per_node, contacts_transmission_m1)
    infection_gens <- c(infection_gens, gen_infected_nodes)
    
  }
  
  
  no_infected_tracing <- table(infection_gens)
  r0_tracing <- no_infected_tracing[epidemic_generations]/no_infected_tracing[epidemic_generations-1]
  
  r0_tracing_list_t_m1 <- c(r0_tracing_list_t_m1,r0_tracing)
  cat("d = ", prob_detection, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
}




r0_p0_list<-c()
for(prob_tracing in t_interval){
  r0_p0 <- r0_tracing_list_t_m1[1]*(1-prob_detection*prob_tracing)
  
  r0_p0_list<- c(r0_p0_list, r0_p0)
  
}


r0_tracing_list_t_m1[length(r0_tracing_list_t_m1)]<-r0_tracing

df_r0_t_m1 <- data.frame(x=t_interval, y= r0_tracing_list_t_m1)

df_r0_t_m2 <- data.frame(x=t_interval, y= r0_tracing_list_t_m2)

df_r0p0_t <- data.frame(x=t_interval, y= r0_p0_list)

df_r0_t_m2$dataset <- "Model 2"
df_r0_t_m1$dataset <- "Model 1"
df_r0p0_t$dataset <- "no clustering"

combined_df<- rbind( df_r0_t_m1, df_r0_t_m2, df_r0p0_t)

# Create the ggplot object
ggplot(data = combined_df, aes(x = x, y = y, color = dataset)) +
  geom_line() +
  geom_point(data = df_r0_t_m1, aes(x = x, y = y))+
  geom_point(data = df_r0_t_m2, aes(x = x, y = y))+
  scale_color_manual(values = c("Model 1" ='#D60270', "Model 2" = '#9B4F96', "no clustering" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.11, 0.24))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0.5 ,1.75)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Tracing Probability (t)", y = "Reproduction Number")


