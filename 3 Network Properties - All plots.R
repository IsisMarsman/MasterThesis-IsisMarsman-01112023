#Analytical Figures

library(ggplot2)


get_new_degree_pois <- function(lambda, p){
  newdeg <- c(0)
  N<- 40
  for (k in 1: N){
    sum <- 0
    
    for (i in 1: k)         #current degree
    { 
      M <- i * (N -1)
      for(l in (k-i):M)                          #number of grandchildren
      { 
        sum <- sum + dpois(i-1, lambda)* dpois(l-1, i*lambda) *dbinom(k-i, l, p)
      } 
      
    }
    newdeg[k+1] <- sum
    cat(k)
    
  }
  
  return(newdeg)
}




get_new_degree_constant <- function(lambda, p){
  newdeg <- c(0)
  
  i <- lambda +1
  l <- lambda *(lambda +1)
  
  N <- i + l +1
  
  for (k in 1: N){
    
    newdeg[k+1] <- dbinom(k-i, l, p)
    cat(k)
    
  }
  
  return(newdeg)
}




get_new_degree_nb <- function(r, q, p){
  newdeg <- c(0)
  maxdegree <- 100                               #say that max degree of poison original distr is 30
  #N <- maxdegree + maxdegree*(maxdegree-1)     #maximal possible new degree
  N <- 100                                       #setting maxdegree and N higher has only effect at 1e-8
  
  for (k in 1: maxdegree){
    sum <- 0
    
    for (i in 1: k)                              #current degree
    { 
      M <- min(N, i * (maxdegree -1))
      for(l in (k-i):N)                          #number of grandchildren
      { 
        sum <- sum + dbinom(k-i, l, p)*  dnbinom( l-1, size= i*r+1, prob=q) * dnbinom(i-1, size = r, prob= q)
        
      } 
      
    }
    newdeg[k+1] <- sum
    cat(k)
    
  }
  
  return(newdeg)
}






#constant

df_const <- data.frame(x=seq(0,11, by=1), y= c( c(get_new_degree_constant(2,0.25),0), c(get_new_degree_constant(2,0.5),0), c(get_new_degree_constant(2,0.75),0)), 
                 ClusteringProbability = c( rep("p = 0.25", 12), rep("p = 0.5",12),rep("p = 0.75", 12)))




ggplot(df_const, aes(x=x, y=y,fill=ClusteringProbability)) +
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("p = 0.25" ='#FFA1D1', "p = 0.5" = '#EB52A1',"p = 0.75" =  '#D60270'))+ 
  theme_light()+
  theme(legend.position = c(0.1, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  scale_x_continuous(breaks = seq(0, 12, by = 1))+
  theme(
    panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    #panel.grid.minor.y = element_blank()
  )+
  labs(x = "k", y = "P(Y=k)")
  


  






#Pois



df_pois <- data.frame(x=seq(0,30, by=1), y= c(get_new_degree_pois(2,0)[1:31], get_new_degree_pois(2,0.5)[1:31], get_new_degree_pois(2,1)[1:31]), 
                      ClusteringProbability = c( rep("p = 0", 31), rep("p = 0.5",31), rep("p = 1",31)))






ggplot(df_pois, aes(x=x, y=y,fill=ClusteringProbability)) +
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("p = 0" ='#FFBAFA', "p = 0.5" = '#CD85C8', "p = 1" =  '#9B4F96'))+
  theme_light()+
  theme(legend.position = c(0.85, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  scale_x_continuous(breaks = seq(0, 30, by = 5))+
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05))+
  #xlim(0,30)+
  theme(
    panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )+
  labs(x = "k", y = "P(Y=k)")










#NB



df_nb <- data.frame(x=seq(0,40, by=1), y= c(get_new_degree_nb(2, 0.5 ,1)[1:41], get_new_degree_nb(2, 0.5, 0.5)[1:41], get_new_degree_nb(2, 0.5, 0)[1:41]), 
                    ClusteringProbability = c( rep("p = 1", 41), rep("p = 0.5", 41), rep("p = 0", 41)))




ggplot(df_nb, aes(x=x, y=y,fill=ClusteringProbability)) +
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("p = 0" ='#91B6FF', "p = 0.5" = '#4977D4', "p = 1" =  '#0038A8'))+
  theme_light()+
  theme(legend.position = c(0.85, 0.8))+
  theme(legend.title = element_blank(),
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  scale_x_continuous(breaks = seq(0, 40, by = 5))+
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.05))+
  #xlim(0,40)+
  theme(
    panel.grid.major.x = element_blank(),
    #panel.grid.minor.x = element_blank(),
    #panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )+
  labs(x = "k", y = "P(Y=k)")







#Expected degree



#expected degree
mean_pois <- c()
mean_constant <- c()
mean_nb <- c()

for(p in seq(0,1, by=0.05 )){
  newdeg_pois <- get_new_degree_pois(2, p)
  value_pois <- 0 : (length(newdeg_pois)-1)
  
  newdeg_constant <- get_new_degree_constant(2, p)
  value_constant <- 0 : (length(newdeg_constant)-1)
  
  newdeg_nb <- get_new_degree_nb(2, 0.5, p)
  value_nb <- 0 : (length(newdeg_nb)-1)
  
  mean_pois <- c(mean_pois, weighted.mean(value_pois, newdeg_pois))
  mean_constant <- c(mean_constant, weighted.mean(value_constant, newdeg_constant))
  mean_nb <- c(mean_nb, weighted.mean(value_nb, newdeg_nb))
  
}



df_expdeg <- data.frame(x=seq(0,1, by=0.05), y= c( mean_constant, mean_pois, mean_nb), Offspring = c(rep("Fixed 2",21), rep("Pois(2)", 21), rep("NB(2,0.5)",21)))



ggplot(df_expdeg, aes(x=x, y=y, color=Offspring)) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = c("Fixed 2" =  '#D60270', "Pois(2)" ='#9B4F96', "NB(2,0.5)" =  '#0038A8' ))+ 
  theme_light()+
  theme(legend.position = c(0.85, 0.2))+
  theme(
        legend.spacing.y = unit(0, "mm"), 
        panel.border = element_rect(colour = "black", fill=NA),
        #aspect.ratio = 1, axis.text = element_text(colour = 1, size = 12),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))+
  #ylim(3 ,11.5)+
  scale_y_continuous(breaks = seq(3, 12, by = 2))+
  labs(x = "Clustering Probability (p)", y = "Expected Degree (E(Y))")







#clustercoeff voor p

clusteringcoeff_pois <- c()
clusteringcoeff_constant <- c()

clusteringcoeff_pois3 <- c()
clusteringcoeff_constant3 <- c()

clusteringcoeff_pois4 <- c()
clusteringcoeff_constant4 <- c()

clusteringcoeff_pois5 <- c()
clusteringcoeff_constant5 <- c()



p <- seq(0,1, by=0.01)
lambda<-2
cluster_pois2 <- (p*(lambda*(3*lambda + (lambda^2 + lambda + 2)*p^2 +4)+2))/ 
  (lambda*(lambda+2)+lambda*(lambda+2)*(lambda^2 + lambda+1)*p^2 +2*(lambda*(lambda+1)*(lambda+2)+1)*p)
cluster_constant2 <- (p* ((lambda - 1)* p^2 + 3))/((lambda^2 + lambda - 1)*p^2 + 2*(lambda + 1)*p + 1) 

lambda<-3
cluster_pois3 <- (p*(lambda*(3*lambda + (lambda^2 + lambda + 2)*p^2 +4)+2))/ 
  (lambda*(lambda+2)+lambda*(lambda+2)*(lambda^2 + lambda+1)*p^2 +2*(lambda*(lambda+1)*(lambda+2)+1)*p)
cluster_constant3 <- (p* ((lambda - 1)* p^2 + 3))/((lambda^2 + lambda - 1)*p^2 + 2*(lambda + 1)*p + 1) 

lambda<-5
cluster_pois5 <- (p*(lambda*(3*lambda + (lambda^2 + lambda + 2)*p^2 +4)+2))/ 
  (lambda*(lambda+2)+lambda*(lambda+2)*(lambda^2 + lambda+1)*p^2 +2*(lambda*(lambda+1)*(lambda+2)+1)*p)
cluster_constant5 <- (p* ((lambda - 1)* p^2 + 3))/((lambda^2 + lambda - 1)*p^2 + 2*(lambda + 1)*p + 1) 




df_cluster <- data.frame(x=p, y= c(cluster_pois2, cluster_constant2, cluster_pois3, cluster_constant3,  cluster_pois5, cluster_constant5), 
                         Offspring = c( rep("Pois(2)", length(p)), rep("Fixed 2", length(p)), rep("Pois(3)", length(p)), rep("Fixed 3", length(p)),  rep("Pois(5)", length(p)), rep("Fixed 5",length(p))))



ggplot(df_cluster, aes(x = x, y = y, color = Offspring)) +
  geom_line() +
  scale_color_manual(values = c("Pois(2)" = '#147BD1', "Fixed 2" = '#753BBD', "Pois(3)" ='#F7EA48' , "Fixed 3" =  '#2DC84D',  
                                "Pois(5)" = '#E03C31', "Fixed 5" ='#FF7F41'             )) +
  theme_light() +
  theme(legend.position = c(0.85, 0.2),
        legend.spacing.y = unit(0, "mm"),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  labs(x = "Clustering Probability (p)", y = "Clustering Coefficient (C)") +
  guides(color = guide_legend(ncol = 2))  



