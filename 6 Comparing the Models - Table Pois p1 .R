library(Matrix)

#Model 2

set.seed(031020233)

M <- 10000                       #number of iterations
population_generations <-14  #generations in the tree, not counting index case as generation
epidemic_generations <- floor(population_generations / 2) +1  #generations in the epidemic (+1 because index case is first generation in the list)

lambda <- 2
prob_clustering_m2 <- 1              #probability to add the extra edges
prob_transmission_m2 <-0.201

prob_tracing_m2 <- 0.5
prob_detection_m2 <- 0.5

sumnoinf <- rep(0, epidemic_generations)
sumnoinf_tracing <- rep(0, epidemic_generations)
list_clustercoefficient_tracing_m2 <- c()
list_clustercoefficient_contact_m2 <- c()
list_degree_tracing_m2 <- c()
list_offspring_transmission_m2 <-c()


for(j in 1:M){
  
  offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
  N <- number_of_nodes(offspring_distr_m2)
  while(N< 1000){
    offspring_distr_m2 <- creating_offspringdistr_pois(lambda, population_generations)
    N <- number_of_nodes(offspring_distr_m2)
  }
  contacts_m2 <- creating_contacts_DK(offspring_distr_m2, N, prob_clustering_m2)
  matrix_final_m2 <- creating_final_matrix(contacts_m2, N, prob_transmission_m2)
  
  epidemic <- sim_epidemic(N, epidemic_generations, matrix_final_m2)
  epidemic_tracing <- sim_epidemic_tracing_m2_v2(N, matrix_final_m2, contacts_m2, prob_tracing_m2, prob_detection_m2)
  
  
  
  sumnoinf<- sumnoinf + epidemic$n_infected  
  sumnoinf_tracing <- sumnoinf_tracing + epidemic_tracing
  
  
  
  contact_network_m2 <- get_contact_network_m2(contacts_m2, N)
  
  
  
  #tracing
  clustercoefficient_tracing_m2 <- get_clusteringcoefficient_tracing_m2_v2(contact_network_m2 , epidemic$infected, epidemic$n_infected, epidemic_generations-3)
  clustercoefficient_contact_m2 <- get_clusteringcoefficient_contact_m2(contact_network_m2, N)
  
  degree_tracing_m2 <- get_degreedistr_tracing_m2(contact_network_m2,epidemic_generations, epidemic$susceptibles, epidemic$infected )
  
  list_clustercoefficient_contact_m2 <- c(list_clustercoefficient_contact_m2, clustercoefficient_contact_m2)
  list_clustercoefficient_tracing_m2 <- c(list_clustercoefficient_tracing_m2, clustercoefficient_tracing_m2)
  list_degree_tracing_m2 <- c(list_degree_tracing_m2, degree_tracing_m2)
  
  
  #transmission
  offspring_transmission_m2 <- get_offspring_transmission_network_m2(contact_network_m2, epidemic$infected, epidemic_generations)
  list_offspring_transmission_m2 <- c(list_offspring_transmission_m2, offspring_transmission_m2)
  
  
  
  cat(j)
  
  
}



#calculating R0
r0 <- sumnoinf[epidemic_generations]/sumnoinf[epidemic_generations-1]
r0_tracing <- sumnoinf_tracing[epidemic_generations]/sumnoinf_tracing[epidemic_generations -1]



r0_tracing_m2<- c(0)
r0_epidemic_m2 <- c(0)
for (i in 1:(epidemic_generations-1)) {
  r0_tracing_m2[i+1] <- sumnoinf_tracing[i+1]/sumnoinf_tracing[i]
  r0_epidemic_m2[i+1] <- sumnoinf[i+1]/sumnoinf[i]
}

r0_tracing_m2
r0_epidemic_m2



expected_deg_transmission_m2 <- r0 +1
avg_cluster_tracing_m2 <- mean(list_clustercoefficient_tracing_m2, na.rm = TRUE)
avg_cluster_contact_m2 <- mean(list_clustercoefficient_contact_m2, na.rm = TRUE)
avg_deg_tracing_m2 <- mean(list_degree_tracing_m2 ,na.rm = TRUE )
var_offspring_transmission_m2 <- var(list_offspring_transmission_m2, na.rm = TRUE )





#Model 1

prob_tracing_m1 <- 0.5
prob_detection_m1 <- 0.5

prob_clustering_m1 <- (avg_deg_tracing_m2 - expected_deg_transmission_m2) / ((r0+1)*(r0) + var_offspring_transmission_m2/r0)

nb_size_m1 <- r0^2/(var_offspring_transmission_m2 - r0)
nb_prob_m1 <- r0/var_offspring_transmission_m2


list_clustercoefficient_tracing_m1 <- c()
list_degree_tracing_m1 <- c()
infection_gens <- c()



for (j in 1:M) {
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
  gen_infected_nodes <- sim_epidemic_tracing_m1(N, contacts_tracing_m1, prob_tracing_m1, prob_detection_m1, generation_per_node, contacts_transmission_m1)
  
  infection_gens <- c(infection_gens, gen_infected_nodes)
  
  clustercoefficient_tracing_m1 <- get_clusteringcoefficient_tracing_m1(matrix_tracing_m1, offspring_distr_m1$no_nodes, population_generations-3 )
  degreedistr_tracing_m1 <- get_degreedistr_tracing_m1(matrix_tracing_m1, population_generations+1 , offspring_distr_m1$no_nodes) 
  
  list_clustercoefficient_tracing_m1 <- c(list_clustercoefficient_tracing_m1, clustercoefficient_tracing_m1)
  list_degree_tracing_m1 <- c(list_degree_tracing_m1, degreedistr_tracing_m1)
  
  cat( j )
  
}

no_infected_tracing <- table(infection_gens)

r0_tracing_m1 <- c(0)
for (i in 1:(epidemic_generations-1)) {
  r0_tracing_m1[i+1] <- no_infected_tracing[i+1]/no_infected_tracing[i]
}



avg_clustercoefficient_tracing_m1 <- mean(list_clustercoefficient_tracing_m1 , na.rm = TRUE)
avg_degree_tracing_m1 <- mean(list_degree_tracing_m1, na.rm = TRUE)

prob_clustering_m1
avg_clustercoefficient_tracing_m1
avg_degree_tracing_m1
r0_tracing_m1


expected_deg_transmission_m2
r0
mean(list_offspring_transmission_m2)
var_offspring_transmission_m2
avg_deg_tracing_m2
r0_tracing_m2
r0_epidemic_m2
avg_cluster_tracing_m2
avg_cluster_contact_m2

