
#model 2

#creating contact network



#Model 2 constant 2 offspring

creating_fixed_contacts_constant2 <- function(population_generations){ 
  
  N <- 3*2^(population_generations) - 2           #number of nodes
  
  #contacts in the branching process with offspring degree 2
  i=2:(3*2^(population_generations - 1) - 2)
  contacts_fixed <- matrix(c(rep(1, 3),i, i, c(2,3,4), 2*i+1,  2*i+2), ncol = 2)
  return(contacts_fixed)
}


creating_extra_contacts_constant2<- function(population_generations){ 
  #contacts with distance 2
  i=2:(3*2^(population_generations - 2) - 2)
  k = seq(3, (3*2^(population_generations - 1) - 2), by =2)
  contacts_extra<- matrix( c( c(2,2), rep(1, 6), i, i, i, i, k , c(3,4), c(5,6,7,8,9,10), 4*i+3, 4*i+4, 4*i+5, 4*i+6, k+1), ncol = 2 )
  return(contacts_extra)
}


creating_contacts_constant2<- function(contacts_fixed, contacts_extra, prob_clustering){
  #apply cluster probability 
  prop_p <-  sample(c(TRUE, FALSE), nrow(contacts_extra), prob = c(prob_clustering, 1 - prob_clustering), replace = TRUE)
  contacts_extra_p <- contacts_extra[prop_p,]
  
  contacts <- rbind(contacts_fixed, contacts_extra_p) #join the two matrices with original contacts and the extra contacts
  
  return(contacts)
}


number_of_nodes_constant2 <- function(population_generations){
  N <- 3*2^(population_generations) - 2
  return(N)
}



#model 2 poison & nb offspring

creating_offspringdistr_pois <- function(lambda, population_generations){
  offspring_distr<- c( 1+rpois(1, lambda) )
  n<- offspring_distr[1]
  
  for (l in 1:population_generations-1) {
    offspring_distr_generation <- rpois(n, lambda)
    n<- sum(offspring_distr_generation)
    offspring_distr<- c(offspring_distr, offspring_distr_generation)
  }
  return(offspring_distr)
}

creating_offspringdistr_nb <- function(size, prob, population_generations){
  offspring_distr<- c( 1+rnbinom(1, size, prob) )
  n<- offspring_distr[1]
  
  for (l in 1:population_generations-1) {
    offspring_distr_generation <- rnbinom(n, size, prob)
    n<- sum(offspring_distr_generation)
    offspring_distr<- c(offspring_distr, offspring_distr_generation)
  }
  return(offspring_distr)
}

number_of_nodes <- function(offspring_distr){
  N<- sum(offspring_distr) + 1
  return(N)
}

creating_contacts_DK <- function(offspring_distr, N, prob_clustering){
  parent <- c(0L, rep(1:length(offspring_distr), offspring_distr))
  
  child <- seq(1, N)
  
  contacts_fixed <- cbind(parent[-1],child[-1])
  
  grandchild <- seq(2, N)
  grandparent <- parent[parent[grandchild]]
  
  #siblings
  firstsibling_nodes <- which(!duplicated(parent))[-1] #unique ivm offspring 0
  familysizes <- offspring_distr[offspring_distr > 0]
  firstsibling_nodes <- firstsibling_nodes[familysizes > 1] #family size 1 has no siblings
  familysizes <- familysizes[familysizes > 1]
  
  siblingpairs <- sapply(2:max(familysizes), function(x)
    outer(combn(x, 2), firstsibling_nodes[familysizes == x], '+') - 1
  )
  
  siblingpairs_1 <- do.call(c, sapply(siblingpairs, function(x) as.vector(x[1, , ])))
  siblingpairs_2 <- do.call(c, sapply(siblingpairs, function(x) as.vector(x[2, , ])))
  
  contacts_extra <- cbind(c(tail(grandparent, -offspring_distr[1]), siblingpairs_1) , 
                          c(tail(grandchild, -offspring_distr[1]), siblingpairs_2))
  
  #apply probability p
  prob_p <-  sample(c(TRUE, FALSE), nrow(contacts_extra), prob = c(prob_clustering, 1 - prob_clustering), replace = TRUE)
  contacts_extra_p <- contacts_extra[prob_p,]
  
  contacts <- rbind(contacts_fixed, contacts_extra_p) #join the two matrices with original contacts and the extra contacts
  
  return(contacts)
}  




#epidemic
creating_final_matrix <- function(contacts, N, prob_transmission){
  
  #apply transmission probability q
  prob_q <- sample(c(TRUE, FALSE), nrow(contacts), prob = c(prob_transmission, 1 - prob_transmission), replace = TRUE)
  contacts_q <- contacts[prob_q,]
  
  #make sparse adjacency matrix
  matrix_final <- Matrix::sparseMatrix(i = contacts_q[,1] , j = contacts_q[,2], dims = c(N,N), symmetric = TRUE) 
  
  return(matrix_final)
}  


sim_epidemic <- function( N, epidemic_generations, m_final ){
  
  infected <- list()                         #list to keep track of the infected per generation
  noinfected <- c()                          #vector with number of infected per generation
  infected[[1]]<- c(TRUE, rep(FALSE, N-1))   #initial case
  susc <- c(FALSE, rep( TRUE, N-1))          #keeping track of the susceptible nodes
  
  
  for (gen in 2:epidemic_generations) {
    infected[[gen]]<- m_final%*%infected[[gen-1]]  & susc
    susc <- susc - infected[[gen]]
  }
  
  
  for (gen in 1:epidemic_generations) {
    noinfected[gen]<- sum(infected[[gen]],na.rm = TRUE)
  }
  
  return(list(infected = infected, susceptibles = susc, n_infected = noinfected))
}



#resulting tracing and transmission networks


get_contact_network_m2 <- function(contacts_m2, N){
  contact_network_m2 <- Matrix::sparseMatrix(i = contacts_m2[,1] , j = contacts_m2[,2], x=1, dims = c(N,N), symmetric = TRUE)
  return(contact_network_m2)
}

get_tracing_network_m2 <- function(contact_network_m2, susceptibles){ 
 
  allinfected <- !susceptibles    
  inf_nodes <- which(as.logical(allinfected)) 
  m_tracing <- contact_network_m2[inf_nodes,inf_nodes]
  
  return(m_tracing)
}



get_degreedistr_tracing_m2 <- function(contact_network_m2, epidemic_generations, susceptibles, infected){
  
  inf_nodes <- which(!susceptibles)
  nodes_interested <- which(infected[[epidemic_generations-3]])
  if(length(nodes_interested) == 0){
    degreedistr_tracing <- c()}else{
      if(length(nodes_interested) == 1){
        m_tracing <- contact_network_m2[inf_nodes,nodes_interested]
        
        degreedistr_tracing <- sum(m_tracing)
      }else{
        m_tracing <- contact_network_m2[inf_nodes,nodes_interested]
        
        degreedistr_tracing <- Matrix::colSums(m_tracing)
      }
     }  
  
  return(degreedistr_tracing)
}


get_expected_degreedistr_tracing_m2 <- function(contact_network_m2, epidemic_generations, susceptibles, infected){
  
  inf_nodes <- which(!susceptibles)
  nodes_interested <- which(infected[[epidemic_generations-3]])
  m_tracing <- contact_network_m2[inf_nodes,nodes_interested]
  
  degreedistr_tracing <- Matrix::colSums(m_tracing)
  expected_deg_tracing <- mean(degreedistr_tracing, na.rm = TRUE)     
  
  return(expected_deg_tracing)
}




get_clusteringcoefficient_tracing_m2<- function(contact_network_m2, infected, susceptibles, generation){
  
  
  susc <- which(as.logical(susceptibles))
  matrix_tracing_m2 <- contact_network_m2
  if(length(susc) !=0){
      matrix_tracing_m2[susc,] <- FALSE
      matrix_tracing_m2[,susc] <- FALSE
    }
  
  A2<- matrix_tracing_m2%*%matrix_tracing_m2
  A3<- A2%*%matrix_tracing_m2
  
  
  inf_nodes <- which(!susceptibles)
  nodes_interested <- which(infected[[generation]])
  m_tracing <- contact_network_m2[inf_nodes,nodes_interested]
  if(length(nodes_interested) == 0){  cluster<- c()
  }else{
    if(length(nodes_interested) == 1){
      degreedistr_tracing <- sum(m_tracing)
      }else{
      degreedistr_tracing <- Matrix::colSums(m_tracing)
      }
    
    triplets <- choose(degreedistr_tracing, 2)
    all_triplets <- sum(triplets ,na.rm = TRUE)
    all_triangles <- sum(Matrix::diag(A3[nodes_interested,nodes_interested]))/2
    cluster <- all_triangles / all_triplets
    }
  
  return(cluster)
}





get_clusteringcoefficient_tracing_m2_v2<- function(contact_network_m2, infected, n_infected, generation){
  
  
  if(n_infected[[generation]] ==0){cluster<- c()}else{
  
  inf<-c()
  for(i in 1:epidemic_generations){
    inf <- c(inf, which(infected[[i]]))
  }
  if(length(inf) !=0){
    matrix_tracing_m2 <- contact_network_m2[inf,inf]
  }
  
  A2<- matrix_tracing_m2%*%matrix_tracing_m2
  A3<- A2%*%matrix_tracing_m2
  
  
  nodes_interested <- (sum(n_infected[1:(generation)-1])+1):(sum(n_infected[1:generation])) 
  
  m_tracing <- matrix_tracing_m2[,nodes_interested]
  if(length(nodes_interested) == 0){  cluster<- c()
  }else{
    if(length(nodes_interested) == 1){
      degreedistr_tracing <- sum(m_tracing)
    }else{
      degreedistr_tracing <- Matrix::colSums(m_tracing)
    }
    
    triplets <- choose(degreedistr_tracing, 2)
    all_triplets <- sum(triplets ,na.rm = TRUE)
    all_triangles <- sum(Matrix::diag(A3[nodes_interested,nodes_interested]))/2
    cluster <- all_triangles / all_triplets
  }
  }
  return(cluster)
}





get_clusteringcoefficient_contact_m2<- function(contact_network_m2, N){
  
  
    A2<- contact_network_m2%*%contact_network_m2
    A3<- A2%*%contact_network_m2
    
    nodes_interested <- 100:150 
    
    m_contact <- contact_network_m2[,nodes_interested]
    if(length(nodes_interested) == 0){  cluster<- c()
    }else{
      if(length(nodes_interested) == 1){
        degreedistr_contact <- sum(m_contact)
      }else{
        degreedistr_contact <- Matrix::colSums(m_contact)
      }
      
      triplets <- choose(degreedistr_contact, 2)
      all_triplets <- sum(triplets ,na.rm = TRUE)
      all_triangles <- sum(Matrix::diag(A3[nodes_interested,nodes_interested]))/2
      cluster <- all_triangles / all_triplets
    }

  return(cluster)
}




create_transmission_network_gen_m2 <- function(contact_network_m2, infected, susceptibles, epidemic_generations){
  infector <- function(infected_node, epi_gen ){
    x<- which(contact_network_m2[,infected_node] & infected[[epi_gen-1]])
    if(length(x)==1){
      return(c(x,infected_node))
    }
    else{
      return(c(sample(x,1),infected_node))
    }
  }
  
  trans_contacts2<- unlist( lapply(which(infected[[epidemic_generations]] == TRUE), infector, epi_gen=epidemic_generations))
  
  
  if(is.null(trans_contacts2)){
    degdistr_transmission_m2 <- c()
  }else{
    trans_contacts2 <- matrix(trans_contacts2, ncol = 2, byrow = TRUE)
    
    matrix_transmission <- Matrix::sparseMatrix(i = trans_contacts2[,1] , j = trans_contacts2[,2],x=1, dims = c(N,N)) 
    matrix_transmission <- matrix_transmission + Matrix::t(matrix_transmission)
    
    degdistr_transmission2 <- Matrix::rowSums(matrix_transmission)
    nodes_interested <- which(infected[[epidemic_generations-1]])
    degdistr_transmission_m2 <- degdistr_transmission2[nodes_interested]
  }
  
  return(degdistr_transmission_m2)
}



get_offspring_transmission_network_m2 <- function(contact_network_m2, infected, epidemic_generations){
  transmission_contacts <- c()
  for (infected_node in which(infected[[epidemic_generations]])){
    infectors <- which(contact_network_m2[,infected_node] & infected[[epidemic_generations-1]])
    if(length(infectors)==1){
      transmission_contacts <- c(transmission_contacts, c(infectors ,infected_node))
    }
    else{
      transmission_contacts <- c(transmission_contacts, c(sample(infectors,1),infected_node))
    }
  }
  
  if(is.null(transmission_contacts)){
    offspring_transmission <- c()
  }else{
    transmission_contacts <- matrix(transmission_contacts, ncol = 2, byrow = TRUE)
    
    matrix_transmission <- Matrix::sparseMatrix(i = transmission_contacts[,1] , j = transmission_contacts[,2],x=1, dims = c(N,N)) 
    matrix_transmission <- matrix_transmission + Matrix::t(matrix_transmission)
    
    offspring_transmission <- Matrix::rowSums(matrix_transmission)
    nodes_interested <- which(infected[[epidemic_generations-1]])
    offspring_transmission <- offspring_transmission[nodes_interested]
  }
  
  return(offspring_transmission)
}






create_transmission_network_m2 <- function(contact_network_m2, infected, susceptibles, epidemic_generations){
  infector <- function(infected_node, epi_gen ){
    x<- which(contact_network_m2[,infected_node] & infected[[epi_gen-1]])
    if(length(x)==1){
      return(c(x,infected_node))
    }
    else{
      return(c(sample(x,1),infected_node))
    }
  }
  
  trans_contacts2 <-c()
  for(i in 2:epidemic_generations){
    trans_contacts2<- c(trans_contacts2, unlist( lapply(which(infected[[i]] == TRUE), infector, epi_gen=i)))
    
  }
  
  if(is.null(trans_contacts2)){
    matrix_transmission <- Matrix::sparseMatrix(i =c() , j = c(), dims = c(N,N))
  }else{
    trans_contacts2 <- matrix(trans_contacts2, ncol = 2, byrow = TRUE)
    
    matrix_transmission <- Matrix::sparseMatrix(i = trans_contacts2[,1] , j = trans_contacts2[,2],x=1, dims = c(N,N)) 
    matrix_transmission <- matrix_transmission + Matrix::t(matrix_transmission)
  }
  
  return(matrix_transmission)
}



get_degreedistribution_transmission_m2 <- function(matrix_transmission, infected, epidemic_generations){
  if(sum(matrix_transmission[1,])==0){degdistr_transmission_m2 <- c()} else{ 
    degdistr_transmission2 <- Matrix::rowSums(matrix_transmission)
    nodes_interested <- which(infected[[epidemic_generations-1]])
    degdistr_transmission_m2 <- degdistr_transmission2[nodes_interested]
  }
  
  return(degdistr_transmission_m2)
}









#Model 1

creating_offspringdistr_degdistr_m1 <- function(degdistr_transmission2, population_generations){
  offspring_distr<- 1+ sample(degdistr_transmission2, 1)
  n<- offspring_distr[1]
  no_nodes<- c(n)
  
  for (l in 1:population_generations-1) {
    offspring_distr_gen <- sample(degdistr_transmission2 -1 , n, replace = TRUE)
    n<- sum(offspring_distr_gen)
    offspring_distr<- c(offspring_distr, offspring_distr_gen)
    no_nodes <- c(no_nodes,n)
  }
  
  return(list(offspring_distr = offspring_distr, no_nodes = no_nodes))
}



creating_offspringdistr_pois_m1 <- function(lambda, population_generations){
  offspring_distr<- c( 1+rpois(1, lambda) )
  n<- offspring_distr[1]
  no_nodes<- c(1,n)
  
  for (l in 2:population_generations-1) {
    offspring_distr_generation <- rpois(n, lambda)
    n<- sum(offspring_distr_generation)
    offspring_distr<- c(offspring_distr, offspring_distr_generation)
    no_nodes <- c(no_nodes,n)
  }
  
  return(list(offspring_distr = offspring_distr, no_nodes = no_nodes))
}


creating_offspringdistr_binom_m1 <- function(size, prob, population_generations){
  offspring_distr<- c( 1+rbinom(1, size, prob) )
  n<- offspring_distr[1]
  no_nodes<- c(1,n)
  
  for (l in 2:population_generations-1) {
    offspring_distr_generation <- rbinom(n, size, prob)
    n<- sum(offspring_distr_generation)
    offspring_distr<- c(offspring_distr, offspring_distr_generation)
    no_nodes <- c(no_nodes,n)
  }
  
  return(list(offspring_distr = offspring_distr, no_nodes = no_nodes))
}


creating_offspringdistr_nb_m1 <- function(size, prob, population_generations){
  offspring_distr<- c( 1+rnbinom(1, size, prob) )
  n<- offspring_distr[1]
  no_nodes<- c(1,n)
  
  for (l in 2:population_generations-1) {
    offspring_distr_generation <- rnbinom(n, size, prob)
    n<- sum(offspring_distr_generation)
    offspring_distr<- c(offspring_distr, offspring_distr_generation)
    no_nodes <- c(no_nodes,n)
  }
  
  return(list(offspring_distr = offspring_distr, no_nodes = no_nodes))
}


creating_offspringdistr_constant_m1 <- function(expected_deg_trans, population_generations){
  offspring_distr<-c(ceiling(expected_deg_trans))
  n<- offspring_distr[1]
  no_nodes<- c(n)
  
  for (l in 1:population_generations-1) {
    offspring_distr_gen <- sample(c(ceiling(expected_deg_trans) -1, floor(expected_deg_trans) -1), n, prob = c(expected_deg_trans - floor(expected_deg_trans), 1 - (expected_deg_trans - floor(expected_deg_trans))), replace = TRUE)
    n<- sum(offspring_distr_gen)
    offspring_distr<- c(offspring_distr, offspring_distr_gen)
    no_nodes <- c(no_nodes,n)
  }
  
  return(list(offspring_distr = offspring_distr, no_nodes = no_nodes))
}



create_transmission_network_m1 <- function(offspring_distr){
  
  parent <- rep(1:length(offspring_distr), offspring_distr)
  N<- sum(offspring_distr) + 1
  child <- seq(2, N)
  
  contacts_transmission1 <- cbind(parent,child)
  
  return(contacts_transmission1)
}







expected_degreedistr_tracing_m1 <- function(m_tracing, epidemic_generations, noinfected){
  
  degreedistr_tracing <- Matrix::rowSums(matrix_tracing)
  nodes_interested <- sum(noinfected) - noinfected[epidemic_generations] - noinfected[epidemic_generations-1] #laatste generatie heeft lage degree
  degreedistr_tracing <- degreedistr_tracing[1:nodes_interested]
  
  expected_deg_tracing <- mean(degreedistr_tracing )     
  
  return(expected_deg_tracing)
}




get_degreedistr_tracing_m1 <- function(m_tracing, epidemic_generations, noinfected){
  
  degreedistr_tracing <- Matrix::rowSums(m_tracing)
  nodes_interested <- sum(noinfected) - noinfected[epidemic_generations] - noinfected[epidemic_generations-1] #laatste generatie heeft lage degree
  degreedistr_tracing <- degreedistr_tracing[1:nodes_interested]
  
  return(degreedistr_tracing)
}


get_clusteringcoefficient_tracing_m1<- function(matrix_tracing_m1, no_nodes, generation){
  A2<-matrix_tracing_m1%*%matrix_tracing_m1
  A3<- A2%*%matrix_tracing_m1
  
  nodes_interested <- seq(from = (2+ sum(no_nodes[1:(generation-1)]) ), length.out = no_nodes[[generation]])
  m_tracing <- matrix_tracing_m1[,nodes_interested]
  degreedistr_tracing <- Matrix::colSums(m_tracing)
  
  triplets <- choose(degreedistr_tracing, 2)
  
  all_triplets <- sum(triplets, na.rm = TRUE)
  all_triangles <- sum(Matrix::diag(A3[nodes_interested,nodes_interested]))/2
  cluster <- all_triangles / all_triplets
  
  return(cluster)
}




get_clusteringcoefficient_tracing_m1_allgen<- function(matrix_tracing_m1, no_nodes, population_generations){
  A2<-matrix_tracing_m1%*%matrix_tracing_m1
  A3<- A2%*%matrix_tracing_m1
  c_gen <- c()
  
  for(i in 2:population_generations){
    nodes_interested <- seq(from = (2+ sum(no_nodes[1:(i-1)]) ), length.out = no_nodes[[i]]) 
    m_tracing <- matrix_tracing_m1[,nodes_interested]
    degreedistr_tracing <- Matrix::colSums(m_tracing)
    
    triplets <- choose(degreedistr_tracing, 2)
    
   
    all_triplets <- sum(triplets, na.rm = TRUE)
    all_triangles <- sum(Matrix::diag(A3[nodes_interested, nodes_interested]))/2
    cluster <- all_triangles / all_triplets
    c_gen <- c(c_gen, cluster)
  }
  
  return(c_gen)
}




#Tracing

#model 2

sim_epidemic_tracing_m2_v2 <- function(N, m_final, contacts_m2, prob_tracing, prob_detection ){
  
  infected_tracing <- list()                      #list to keep track of the infected per generation
  noinfected_tracing <- c()                       #vector with number of infected per generation
  infected_tracing[[1]]<- c(TRUE, rep(FALSE, N-1))   #initial case
  susc_tracing <- c(FALSE, rep( TRUE, N-1))       #keeping track of the susceptible nodes
  
  infected_notfound_tracing <- list()
  infected_notfound_tracing[[1]]<- c(TRUE, rep(FALSE, N-1))
  
  
  detected_nodes <- sample(c(TRUE, FALSE), N, prob = c(prob_detection, 1 - prob_detection), replace = TRUE)
  
  prob_t <- sample(c(TRUE, FALSE), nrow(contacts_m2), prob = c(prob_tracing, 1 - prob_tracing), replace = TRUE)
  contacts_t <- contacts_m2[prob_t,]
  matrix_tracing <- Matrix::sparseMatrix(i = contacts_t[,1] , j = contacts_t[,2], dims = c(N,N), symmetric = TRUE) 
  m_final_tracing <- m_final
  
  for (gen in 2:epidemic_generations) {
    infected_tracing[[gen]]<- m_final_tracing%*%infected_notfound_tracing[[gen-1]]  & susc_tracing
    susc_tracing <- susc_tracing - infected_tracing[[gen]]
    
    infected_notfound_tracing[[gen]] <- infected_tracing[[gen]]
    
    detected <- which( infected_tracing[[gen-1]] & detected_nodes )    #apply detectie prob op infected_tracing[[gen-1]]
    
    if(length(detected) == 0){
      
    }else if(length(detected) == 1){
      traced_contacts <- unique(which(matrix_tracing[detected,], arr.ind = TRUE))  #vind contacten van deze nodes
      if(length(traced_contacts) != 0){
        infected_notfound_tracing[[gen]][traced_contacts] <- FALSE
      }
    }else {
      traced_contacts <- unique(which(matrix_tracing[detected,], arr.ind = TRUE)[,2])  #vind contacten van deze nodes
      if(length(traced_contacts) != 0){
        infected_notfound_tracing[[gen]][traced_contacts] <- FALSE
        
      }
    }
  }
  
  for (gen in 1:epidemic_generations) {
    noinfected_tracing[gen]<- sum(infected_tracing[[gen]],na.rm = TRUE)
  }
  
  return(noinfected_tracing)
  
}


sim_epidemic_tracing_m2_v1 <- function(N, m_final, contacts_m2, prob_tracing, prob_detection ){
  
  infected_tracing <- list()                      #list to keep track of the infected per generation
  noinfected_tracing <- c()                       #vector with number of infected per generation
  infected_tracing[[1]]<- c(TRUE, rep(FALSE, N-1))   #initial case
  susc_tracing <- c(FALSE, rep( TRUE, N-1))       #keeping track of the susceptible nodes
  
  
  detected_nodes <- sample(c(TRUE, FALSE), N, prob = c(prob_detection, 1 - prob_detection), replace = TRUE)
  
  prob_t <- sample(c(TRUE, FALSE), nrow(contacts_m2), prob = c(prob_tracing, 1 - prob_tracing), replace = TRUE)
  contacts_t <- contacts_m2[prob_t,]
  matrix_tracing <- Matrix::sparseMatrix(i = contacts_t[,1] , j = contacts_t[,2], dims = c(N,N), symmetric = TRUE) 
  m_final_tracing <- m_final
  
  for (gen in 2:epidemic_generations) {
    infected_tracing[[gen]]<- m_final_tracing%*%infected_tracing[[gen-1]]  & susc_tracing
    susc_tracing <- susc_tracing - infected_tracing[[gen]]
    
    
    detected <- which( infected_tracing[[gen-1]] & detected_nodes )    #apply detectie prob op infected_tracing[[gen-1]]
    
    if(length(detected) == 0){
      
    }else if(length(detected) == 1){
      traced_contacts <- unique(which(matrix_tracing[detected,], arr.ind = TRUE))  #vind contacten van deze nodes
      if(length(traced_contacts) != 0){
        m_final_tracing[traced_contacts,] <- FALSE    #verwijder de connecties van die contacten in m_final
        m_final_tracing[,traced_contacts] <- FALSE
        matrix_tracing[traced_contacts,] <- FALSE    
        matrix_tracing[,traced_contacts] <- FALSE
        cat("traced:", traced_contacts, "\n")}
    }else {
      traced_contacts <- unique(which(matrix_tracing[detected,], arr.ind = TRUE)[,2])  #vind contacten van deze nodes
      if(length(traced_contacts) != 0){
        m_final_tracing[traced_contacts,] <- FALSE    #verwijder de connecties van die contacten in m_final
        m_final_tracing[,traced_contacts] <- FALSE
        matrix_tracing[traced_contacts,] <- FALSE    
        matrix_tracing[,traced_contacts] <- FALSE
        cat("traced:",traced_contacts, "\n")}
    }
  }
  
  for (gen in 1:epidemic_generations) {
    noinfected_tracing[gen]<- sum(infected_tracing[[gen]],na.rm = TRUE)
  }
  
  return(noinfected_tracing)
}



#model 1


list_generation_per_node_distr <- function(population_generations, no_nodes){
  inf_gen_node<- c(0)
  nodes <- cumsum(no_nodes)
  for(i in 1:population_generations){
    inf_gen_node[(nodes[i] +1) : nodes[i+1]] <- i
  }
  return(inf_gen_node)
}



list_generation_per_node_const2 <- function(population_generations){
  inf_gen_node<- c(0)
  for(i in 1:population_generations){
    inf_gen_node[(3*2^(i-1) - 2 +1) : (3*2^(i) - 2)] <- i
  }
  return(inf_gen_node)
}



get_all_forward_contacts <- function(quarantaine_nodes, population_generations, generation_node, contacts_og){
  forward_contacts <- c()
  forward_contacts_gen <- quarantaine_nodes
  for(gen in generation_node+1 : (population_generations-1)){
    contacts <- unlist(sapply(forward_contacts_gen, function(node) contacts_og[which(contacts_og[,1]==node),2]))
    
    forward_contacts_gen <- contacts
    forward_contacts <- c(forward_contacts, forward_contacts_gen)
  }
  
  return(forward_contacts)
}






sim_epidemic_tracing_m1 <- function(N, contacts_m1, prob_tracing, prob_detection, generation_per_node, contacts_transmission_m1){
  
  prop_t <-  sample(c(TRUE, FALSE), nrow(contacts_m1), prob = c(prob_tracing, 1 - prob_tracing), replace = TRUE)
  contacts_tracing <- contacts_m1[prop_t,]
  
  nodes <- 1:N
  prob_d <- sample(c(TRUE, FALSE), N, prob = c(prob_detection, 1 - prob_detection), replace = TRUE)
  detected_nodes <- nodes[prob_d]
  detected_nodes_check <- detected_nodes
  
  all_prevented_nodes <- c()
  
  while(length(detected_nodes)!=0){
    node<- detected_nodes[1]
    traced_contacts <- contacts_tracing[which(contacts_tracing[,1]==node),2]
    
    quarantaine <- traced_contacts[which(generation_per_node[traced_contacts] == (generation_per_node[node]+1))]
    prevented_by_tracing <- traced_contacts[which(generation_per_node[traced_contacts] == (generation_per_node[node]+2))]
    
    prevented_by_quarantaine <- get_all_forward_contacts(c(quarantaine,prevented_by_tracing), population_generations, generation_per_node[node], contacts_transmission_m1)
    
    all_prevented_nodes <- unique(c(all_prevented_nodes, prevented_by_quarantaine, prevented_by_tracing))
    
    detected_nodes <- setdiff(detected_nodes, c(all_prevented_nodes, node, quarantaine))
  }
  
  infected_nodes <- setdiff(1:N, all_prevented_nodes)
  
  generation_infected_nodes <- generation_per_node[infected_nodes]
  
  return(generation_infected_nodes)
  
}





