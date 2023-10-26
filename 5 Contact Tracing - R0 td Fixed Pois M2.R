
library(ggplot2)
library(Matrix)


set.seed(29092023)

#Fixed 2


M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

prob_transmission <- 0.5
prob_clustering <-0.5

r0_list_t<- c()
r0_tracing_list_t<- c()
sumnoinf <- rep(0, epidemic_generations)
sumnoinf_tracing <- rep(0, epidemic_generations)

N <- number_of_nodes_constant2(population_generations)
contacts_og <- creating_fixed_contacts_constant2(population_generations)
contacts_extra <- creating_extra_contacts_constant2(population_generations)

for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
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


df_r0_t$dataset <- "p=0.5 q=0.5"
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


#plot for t
M <-10000
population_generations <- 14 #generations in the tree, not counting index case as generation

t_interval <- seq(from=0, to=1, length.out=11)
d_interval <- seq(from=0, to=1, length.out=5)

epidemic_generations <- floor(population_generations / 2) +1

prob_transmission <- 0.5
prob_clustering <-0.5
lambda<-2

r0_list_t_pois<- c()
r0_tracing_list_t_pois<- c()
sumnoinf <- rep(0, epidemic_generations)
sumnoinf_tracing <- rep(0, epidemic_generations)


for(prob_detection in d_interval){ 
  for(prob_tracing in t_interval){
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


df_r0_t_pois$dataset <- "p=0.5 q=0.5"
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


