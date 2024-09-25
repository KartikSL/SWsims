
############### Functions for bootstrap simulations (not parallel). Used in the power simulations #############

library(igraph)
library(poweRlaw)
#library(expm)

#Define function to simulate Newman-Watts graphs
simulate_nw <- function(n, d, beta){
  #Generate lattice
  if(beta != 0){
    l <- as_adj(sample_smallworld(1, n, ceiling(d*beta), 0))
  }
  else{
    l <- 0
  }
  #Generate random edges
  e <- as_adj(sample_gnp(n, 2*d*(1 - beta)/(n-1)))
  #Superimpose edges on lattice
  return(graph_from_adjacency_matrix(l + e, "undirected", weighted = TRUE))
}


### NW graph with SBM + lattice

simulate_nw_SBM <- function(n, d, beta){
  #Generate lattice
  l <- as_adj(sample_smallworld(1, n, ceiling(d*beta), 0))
  #Generate random edges
  k=4
  clusters = c(rep(1,n/4),rep(2,n/4),rep(3,n/4),rep(4,n/4))
  phat = matrix(0.5,k,k)
  diag(phat) = rep(2,k)
  ## this makes sure the average degree of the SBM component is 2d(1-beta)
  phat = phat/mean(phat)
  p = 2*(d - ceiling(d*beta))/(n-1)*phat
  P = get_P_SBM(n,k, p,clusters)
  e <- as_adj(sample_from_SBM(n,P))
  #Superimpose edges on lattice
  return(graph_from_adjacency_matrix(l + e, "undirected",weighted=TRUE))
}


### NW graph with DCSBM + lattice
simulate_nw_DCSBM <- function(n, d, beta){
  #Generate lattice
  l <- as_adj(sample_smallworld(1, n, ceiling(d*beta), 0))
  #Generate random edges
  k=4
  clusters = c(rep(1,n/4),rep(2,n/4),rep(3,n/4),rep(4,n/4))
  ## the alpha parameter controls degree heterogeneity...higher values are more homogeneous
  theta = rpldis(n,xmin=1,alpha=3)
  # theta <- 1/runif(n, 1, 2)
  phat = matrix(0.5,k,k)
  diag(phat) = rep(2,k)
    ## this makes sure the average degree of the DCSBM component is 2d(1-beta)
  phat = phat/(mean(phat)*mean(theta))
  p = 2*(d - ceiling(d*beta))/(n-1)*phat
  P = get_P_DCSBM(n,k,theta,p, clusters)
  ### Note some P might be larger than 1 here, but teh function that samples from DCSBM should be able to handle that
  e <- as_adj(sample_from_DCSBM(n,P))
  #Superimpose edges on lattice
  return(graph_from_adjacency_matrix(l + e, "undirected",weighted=TRUE))
}

simulate_nw_DCSBMSDH <- function(n, d, beta){
  #Generate lattice
  l <- as_adj(sample_smallworld(1, n, ceiling(d*beta), 0))
  #Generate random edges
  k=4
  clusters = c(rep(1,n/4),rep(2,n/4),rep(3,n/4),rep(4,n/4))
  ## the alpha parameter controls degree heterogeneity...higher values are more homogeneous
  #theta = rpldis(n,xmin=1,alpha=3)
  theta <- 1/runif(n, 1, 8)
  phat = matrix(0.5,k,k)
  diag(phat) = rep(2,k)
  ## this makes sure the average degree of the DCSBM component is 2d(1-beta)
  phat = phat/(mean(phat)*mean(theta))
  p = 2*(d - ceiling(d*beta))/(n-1)*phat
  P = get_P_DCSBM(n,k,theta,p, clusters)
  ### Note some P might be larger than 1 here, but teh function that samples from DCSBM should be able to handle that
  e <- as_adj(sample_from_DCSBM(n,P))
  #Superimpose edges on lattice
  return(graph_from_adjacency_matrix(l + e, "undirected",weighted=TRUE))
}


simulate_nw_CL <- function(n, d, beta){
  #Generate lattice
  l <- as_adj(sample_smallworld(1, n, ceiling(d*beta), 0))
  #Generate random edges
  ## the alpha parameter controls degree heterogeneity...higher values are more homogeneous
  theta = rpldis(n,xmin=1,alpha=3)
  thetap = (theta/mean(theta))%*%(t(theta)/mean(theta))
  P = 2*(d - ceiling(d*beta))/(n-1)*thetap
  ### Note some P might be larger than 1 here, but teh function that samples from DCSBM should be able to handle that
  e <- as_adj(sample_from_CL(n,P))
  #Superimpose edges on lattice
  return(graph_from_adjacency_matrix(l + e, "undirected",weighted=TRUE))
}

#Adjusted distances for disconnected networks
cl_dist <- function(g){
  if(is_connected(g)){
    return(mean_distance(g))
  }
  dist.mat <- distances(g)
  #Get off diagonal elements
  dist <- dist.mat[row(dist.mat) != col(dist.mat)]
  #Replace inf distances with diamater
  dist[dist == Inf] <- diameter(g)
  return(mean(dist))
}
#Define function to simulate list of NW graphs
simulate_nw_list <- function(m, beta, n, d, null){
  if(null == "ER"){
    return(replicate(m, simulate_nw(n, d, beta), simplify = FALSE))
  }
  else if(null == "SBM"){
    return(replicate(m, simulate_nw_SBM(n, d, beta), simplify = FALSE))
  }
  else if(null == "DCSBM"){
    return(replicate(m, simulate_nw_DCSBM(n, d, beta), simplify = FALSE))
  }
  else if(null == "DCSBMSDH"){
    return(replicate(m, simulate_nw_DCSBMSDH(n, d, beta), simplify = FALSE))
  }
  else if(null == "CL"){
    return(replicate(m, simulate_nw_CL(n, d, beta), simplify = FALSE))
  }
}

#Define function to get C/L
get_t <- function(g){
  L <- mean_distance(g)
  C <- transitivity(g, type = "global")
  return(C/L)
}

#Define function to get MME test statistic
get_t_mme <- function(g){
  g.degree <- degree(g)
  d.bar <- mean(g.degree)
  sd.d <- sd(g.degree)
  n <- gorder(g)
  D <- sqrt((d.bar*(n - 1))^2 - 4*(n - 1)*(d.bar*sd.d)^2)
  y <- (d.bar*(n - 1) - D)/(2*d.bar^2)
  return(1 - y)
}

#Define clustering based statistic 1
get_t_clust_1 <- function(g){
  C <- transitivity(g, type = "global")
  E <- ecount(g)
  n <- gorder(g)
  return(C - 2*E/(n*(n - 1)))
}

#Define clustering based statistic 2
get_t_clust_2 <- function(g){
  C <- transitivity(g)
  E <- ecount(g)
  n <- gorder(g)
  return(0.5*C*n*(n - 1)/E)
}

## Spectral clustering for fitting SBM and DCSBM: using Adjacency matrix, but with projection
spectral<-function(x,n,k)
{
    if(!isSymmetric(x)){
      x <- (x + t(x))/2
    }
    spectra<-eigen(x)
    specmat<-spectra$vectors[,1:k]
    rownorm<-apply(specmat,1,function(a){(sum(a^2))^0.5})
    rownorm<-ifelse(rownorm < 10^(-06),10^(-06),rownorm)
    specnorm<-specmat/rownorm
    speck<-kmeans(specnorm,k,nstart=5)
    return(speck$cluster)
}

#Adjusted distances for disconnected networks
cl_dist <- function(g){
  if(is_connected(g)){
    return(mean_distance(g))
  }
  dist.mat <- distances(g)
  #Get off diagonal elements
  dist <- dist.mat[row(dist.mat) != col(dist.mat)]
  #Replace inf distances with diamater
  dist[dist == Inf] <- diameter(g)
  return(mean(dist))
}

dcspectral<-function(x,n,k)
{ d=rowSums(x)
  d=d+mean(d)
  deg<-diag(1/sqrt(d))
  l=deg%*%x%*%deg
  if(!isSymmetric(l)){
    l <- (l + t(l))/2
  }
  spectra<-eigen(l)
  specmat<-spectra$vectors[,1:k]
  rownorm<-apply(specmat,1,function(a){(sum(a^2))^0.5})
  rownorm<-ifelse(rownorm < 10^(-06),10^(-06),rownorm)
  specnorm<-specmat/rownorm
  speck<-kmeans(specnorm,k,nstart=5)
  return(speck$cluster)
}

get_P_SBM <- function(n, k, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- Z%*%p%*%t(Z)
  return(P)
}

get_P_DCSBM <- function(n, k, theta, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- diag(theta)%*%Z%*%p%*%t(Z)%*%diag(theta)
  return(P)
}

### Sample network from Chung-Lu
sample_from_CL <- function(n,P){
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] <- runif(n*(n - 1)/2)
  A <- A + t(A)
  diag(A) <- runif(n)
  A <- (A < P) + 0
  G <- graph_from_adjacency_matrix(A, mode = "undirected")
  return(G)
}

### Sample network from SBM
sample_from_SBM<-function(n,P)
{
A= matrix(0,n,n)

	for(l in 1:(n-1))
	{
		for(j in (l+1):n)
		{
	    A[l,j]<-rbinom(1,1,P[l,j])
	    A[j,l]<-A[l,j]
	}
	}
	G=graph_from_adjacency_matrix(A,mode="undirected")
return(G)	
}


### Sample network from DCSBM

sample_from_DCSBM<-function(n,P)
{
    A= matrix(0,n,n)
    
    for(l in 1:(n-1))
    {
        for(j in (l+1):n)
        {
            A[l,j]<-ifelse(rpois(1,P[l,j])>=1,1,0)
            A[j,l]<-A[l,j]
        }
    }
    G=graph_from_adjacency_matrix(A)
    return(G)
}




## Note all three functions have number of iterations b set to 10000. It runs on my computer for a small data like Karate club. But for bigger dataset set b to be smaller or use parallel computing

###Compute p-value with Chung-Lu Null model

small_world_CL <- function(G, b){
  n = gorder(G)
  d <- degree(G)
  P <- d %*%t (d)/sum(d)
  
  t <- transitivity(G)
  t_1 <- get_t(G)
  L <- mean_distance(G)
  
  tstar =  rep(NA,b)
  Lstar =  rep(NA,b)
  tLstar =  rep(NA,b)
  
  tick<- proc.time()
cat("Before b iterations ", tick, "\n")
  
  for(iter in 1:b){
    Gstar <- sample_from_CL(n, P)
    Lstar[iter] <- cl_dist(Gstar)
    tstar[iter] <- transitivity(Gstar)
    tLstar[iter] <- transitivity(Gstar, type = "global")/cl_dist(Gstar)
  }
  
  tock <- proc.time() - tick
   cat("After b iterations ", tock, "\n")
  
  pval1 = sum(tstar>t)/b
  pval2 = sum(Lstar<L)/b
  pval3 = sum(tLstar>t_1)/b
  
  return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
              coef1=t, coef2=L, coef3=t_1,
              tstar=tstar,Lstar=Lstar, tLstar=tLstar))
}


### Compute p-value with ER Null model

small_world_ER <-function(G,b=100)
{
	n=gorder(G)
   L = mean_distance(G)
   C= transitivity(G, type = "global")
   #t=sqrt(n)*C/L
   t_1=get_t(G)
   t=transitivity(G)
   phat = gsize(G)/choose(n,2)
   
   tstar =  rep(NA,b)
   Lstar =  rep(NA,b)
   tLstar =  rep(NA,b)

tick<- proc.time()
cat("Before b iterations ", tick, "\n")

   for (iter in 1:b){
   Gstar = sample_gnp(n=n,p=phat)
   Lstar[iter] =mean_distance(Gstar)
   tstar[iter] = transitivity(Gstar)
   tLstar[iter] = transitivity(Gstar, type = "global")/mean_distance(Gstar)
}

tock <- proc.time() - tick
   cat("After b iterations ", tock, "\n") 

pval1 = sum(tstar>t)/b
pval2 = sum(Lstar<L)/b
pval3 = sum(tLstar>t_1)/b

   return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
               coef1=t, coef2=L, coef3=t_1,
               tstar=tstar,Lstar=Lstar, tLstar=tLstar))
}


### Compute p value with SBM null model. For SBM supply number of communities from Louvain method for now. But in future detect number of communities through spectral methods.
small_world_SBM <-function(G,k,b, clusters)
{
		n=gorder(G)
   L = mean_distance(G)
   C= transitivity(G, type = "global")
   #t=sqrt(n)*C/L
   t_1=get_t(G)
   t=transitivity(G)
   A=as_adj(G,type="both",sparse=FALSE)

ones = matrix(1,n,n)
diag(ones) = 0

#tick <- proc.time()
#cat("Computing phat: ", tick, "\n")
   phat = matrix(0,k,k)
   for (i in 1:k)
   {for ( j in i:k)
   {
   phat[i,j] = sum(A[clusters==i,clusters==j])/sum(ones[clusters==i,clusters==j])
   phat[j,i]=phat[i,j]
   }
   }
#tock <- proc.time() - tick
#cat("phat computed. Time taken: ", tock, "\n")
   
   #tick <- proc.time()
#cat("Computing P: ", tick, "\n")
   P <- get_P_SBM(n, k, phat, clusters)
   
   #tock <- proc.time() - tick
#cat("P computed. Time taken: ", tock, "\n")
   
   tstar =  rep(NA,b)
   Lstar =  rep(NA,b)
   tLstar =  rep(NA,b)
   
   #tick <- proc.time()
#cat("Sampling SBM networks: ", tick, "\n")
   
   for (iter in 1:b){
   Gstar = sample_from_SBM(n,P)
   
  #tstar[iter] = sqrt(n)*transitivity(Gstar, type = "global")/mean_distance(Gstar)
   Lstar[iter] =mean_distance(Gstar)
   tstar[iter] = transitivity(Gstar)
   tLstar[iter] = transitivity(Gstar, type = "global")/mean_distance(Gstar)
   }
   #tock <- proc.time() - tick
#cat("Sampling done. Time taken: ", tock, "\n")
  
   #tock <- proc.time() - tick
   #cat("After b iterations ", tock, "\n")
   
   pval1 = sum(tstar>t)/b
   pval2 = sum(Lstar<L)/b
   pval3 = sum(tLstar>t_1)/b
   
   return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
               coef1=t, coef2=L, coef3=t_1,
               tstar=tstar,Lstar=Lstar, tLstar=tLstar))
	
}


### Compute p value with DCSBM null model. 


small_world_DCSBM <-function(G,k,b=100, clusters)
{
    n=gorder(G)
    L = mean_distance(G)
    C= transitivity(G, type = "global")
    t_1=get_t(G)
    t = transitivity(G)
    A=as.matrix(as_adj(G,type="both"))
    
    thetahat=rep(0,n)
    thetau=rowSums(A)
    phat = matrix(0,k,k)
    
    for (i in 1:k){
        sizek = length(thetau[clusters==i])
    thetak=sum(thetau[clusters==i])
    thetahat[clusters==i] = sizek*(thetau[clusters==i])/thetak
    }
    
    for (i in 1:k)
    {for ( j in i:k)
        {
            phat[i,j] = sum(A[clusters==i,clusters==j])/length(A[clusters==i,clusters==j])
            phat[j,i]=phat[i,j]
        }
    }
    
    tstar =  rep(NA,b)
    Lstar =  rep(NA,b)
    tLstar =  rep(NA,b)
    
    P <- get_P_DCSBM(n, k, thetahat, phat, clusters)
    
    tick<- proc.time()
cat("Before b iterations ", tick, "\n")
    
    for (iter in 1:b){
        Gstar = sample_from_DCSBM(n, P)
        Lstar[iter] =mean_distance(Gstar)
        #tstar[iter] = sqrt(n)*transitivity(Gstar, type = "global")/mean_distance(Gstar)
        tstar[iter] = transitivity(Gstar)
        tLstar[iter] = transitivity(Gstar, type = "global")/mean_distance(Gstar)
    }
    
    tock <- proc.time() - tick
   cat("After b iterations ", tock, "\n")
    
    pval1 = sum(tstar>t)/b
    pval2 = sum(Lstar<L)/b
    pval3 = sum(tLstar>t_1)/b
    
    return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
                coef1=t, coef2=L, coef3=t_1,
                tstar=tstar,Lstar=Lstar, tLstar=tLstar))
    
}



