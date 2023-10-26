library(ggplot2)
library(Matrix)

set.seed(041020235)

#MODEL 2
#Fixed 2

M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=5)
p_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_detection <- 0.5
prob_transmission <-0.5


r0_tracing_list_pt_fixed<- c()
sumnoinf_tracing <- rep(0, epidemic_generations)

N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)

for(prob_tracing in t_interval){
  for(prob_clustering in p_interval){ 
    sumnoinf_tracing <- rep(0, epidemic_generations)
    
    
    for(j in 1:M){
      
      contacts_m2 <- creating_contacts_constant2(contacts_og, contacts_extra, prob_clustering)
      matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission)
      epidemic_tracing <- sim_epidemic_tracing_m2_v2(N, matrix_final_m2, contacts_m2, prob_tracing, prob_detection)
      sumnoinf_tracing <- sumnoinf_tracing + epidemic_tracing
      
      
    }
    
    
    
    r0_tracing <- sumnoinf_tracing[epidemic_generations]/sumnoinf_tracing[epidemic_generations-1]
    
    if(is.na(r0_tracing)){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]
    }else{ 
      if(r0_tracing < 1.1){r0_tracing <- sumnoinf_tracing[4]/sumnoinf_tracing[3]}}
    
    r0_tracing_list_pt_fixed <- c(r0_tracing_list_pt_fixed,r0_tracing)
    cat("p = ", prob_clustering, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}


df_r0_pt_fixed <- data.frame(x=p_interval, y= r0_tracing_list_pt_fixed, q = c(rep("t = 0",11),rep("t = 0.25",11), rep("t = 0.5",11), rep("t = 0.75",11),rep("t = 1",11)))


ggplot(df_r0_pt_fixed, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("t = 0" ='#D60270', "t = 0.25" = '#B92983',"t = 0.5" = '#9B4F96',"t = 0.75" = '#4E449F', "t = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.125, 0.8))+
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
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=5)
p_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_detection <- 0.5
prob_transmission <-0.5

lambda <- 2

r0_tracing_list_pt_pois<- c()
sumnoinf_tracing <- rep(0, epidemic_generations)


for(prob_tracing in t_interval){
  for(prob_clustering in p_interval){ 
    sumnoinf_tracing <- rep(0, epidemic_generations)
    
    
    for(j in 1:M){
      
      offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
      N <- number_of_nodes(offspring_distr_m2)
      while(N< 100){
        offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
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
    
    r0_tracing_list_pt_pois <- c(r0_tracing_list_pt_pois,r0_tracing)
    cat("p = ", prob_clustering, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}


df_r0_pt_pois <- data.frame(x=p_interval, y= r0_tracing_list_pt_pois, q = c(rep("t = 0",11),rep("t = 0.25",11), rep("t = 0.5",11), rep("t = 0.75",11),rep("t = 1",11)))


ggplot(df_r0_pt_pois, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("t = 0" ='#D60270', "t = 0.25" = '#B92983',"t = 0.5" = '#9B4F96',"t = 0.75" = '#4E449F', "t = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.125, 0.8))+
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
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=5)
p_interval <- seq(from=0, to=1, length.out=11)

epidemic_generations <- floor(population_generations / 2) +1

prob_detection <- 0.5
prob_transmission <-0.5

nb_size <- 2
nb_prob <- 0.5

r0_tracing_list_pt_nb<- c()
sumnoinf_tracing <- rep(0, epidemic_generations)

N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)

for(prob_tracing in t_interval){
  for(prob_clustering in p_interval){ 
    sumnoinf_tracing <- rep(0, epidemic_generations)
    
    
    for(j in 1:M){
      
      offspring_distr_m2 <- creating_offspringdistr_nb(nb_size, nb_prob, population_generations)
      N <- number_of_nodes(offspring_distr_m2)
      while(N< 100){
        offspring_distr_m2 <- creating_offspringdistr_nb(nb_size, nb_prob, population_generations)
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
    
    r0_tracing_list_pt_nb <- c(r0_tracing_list_pt_nb,r0_tracing)
    cat("p = ", prob_clustering, "t=", prob_tracing, "R0_tracing =", r0_tracing,"\n")
  }
  
}


df_r0_pt_nb <- data.frame(x=p_interval, y= r0_tracing_list_pt_nb, q = c(rep("t = 0",11),rep("t = 0.25",11), rep("t = 0.5",11), rep("t = 0.75",11),rep("t = 1",11)))


ggplot(df_r0_pt_nb, aes(x=x, y=y, color=q)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("t = 0" ='#D60270', "t = 0.25" = '#B92983',"t = 0.5" = '#9B4F96',"t = 0.75" = '#4E449F', "t = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.125, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  ylim(0 ,3)+
  scale_x_continuous(breaks = seq(0, 1, by = 0.2))+
  labs(x = "Clustering Probability (p)", y = "Reproduction Number")





