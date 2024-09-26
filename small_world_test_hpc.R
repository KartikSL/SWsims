
###############3 Functions for bootstrap simulations (parallel). Used for real-world simulations ############

library(igraph)
  library(ggplot2)
  library(cowplot)

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

#Define function to simulate list of NW graphs
simulate_nw_list <- function(m, beta, n, d){
  return(replicate(m, simulate_nw(n, d, beta), simplify = FALSE))
}

## Spectral clustering for fitting SBM and DCSBM: using Adjacency matrix, but with projection

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

get_t_mod <- function(g){
  L <- mean_distance(g)
  C <- get_t_clust_2(g)
  return(C/L)
}


plot_real <- function(real_data){

  
  df <- real_data[[1]]
  coef <- real_data[[2]]
  pval <- real_data[[3]]
  
  p.1 <- ggplot(df, aes(x = df$tstar, color = null_model, fill = null_model)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = coef[1], color = "red") + 
    theme(legend.position="none") + 
    xlab("clustering coefficient")
  
  p.2 <- ggplot(df, aes(x = df$lstar, color = null_model, fill = null_model)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = coef[2], color = "red") + 
    xlab("average path length")
  
  p.3 <- ggplot(df, aes(x = df$tlstar, color = null_model, fill = null_model)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = coef[3], color = "red") + 
    xlab("C/L")
  
  p <- list(plot_grid(p.1, p.2, ncol = 2), p.3)
  
  return(p)
}



refactor_data <- function(data.er, data.sbm, data.dcsbm, data.cl){
  
  coef <- c(data.er$coef1, data.er$coef2, data.er$coef3)
  pval <- list(c(data.er$pval1, data.er$pval2, data.er$pval3), 
               c(data.sbm$pval1, data.sbm$pval2, data.sbm$pval3),
               c(data.dcsbm$pval1, data.dcsbm$pval2, data.dcsbm$pval3),
               c(data.cl$pval1, data.cl$pval2, data.cl$pval3))
  
  er <- rep("ER", length(data.er$tstar))
  sbm <- rep("SBM", length(data.sbm$tstar))
  dcsbm <- rep("DCSBM", length(data.dcsbm$tstar))
  cl <- rep("CL", length(data.cl$tstar))
  
  data.er <- data.frame(null_model = er, tstar = data.er$tstar,
                        lstar = data.er$Lstar, tlstar = data.er$tLstar)
  data.sbm <- data.frame(null_model = sbm, tstar = data.sbm$tstar,
                         lstar = data.sbm$Lstar, tlstar = data.sbm$tLstar)
  data.dcsbm <- data.frame(null_model = dcsbm, tstar = data.dcsbm$tstar,
                           lstar = data.dcsbm$Lstar, tlstar = data.dcsbm$tLstar)
  data.cl <- data.frame(null_model = cl, tstar = data.cl$tstar,
                        lstar = data.cl$Lstar, tlstar = data.cl$tLstar)
  
  df <- rbind(data.er, data.sbm, data.dcsbm, data.cl)
  return(list(df, coef, pval))
}

#Get decisions from resuls object
get_decisions <- function(g.list){
  decisions <- lapply(g.list, function(x){
    lapply(x, function(y){
      p.vals <- unlist(y[1:3])
      thresh <- c(0.05, 0.99, 0.05)
      dec <- p.vals < thresh
      return(c(dec[1], dec[2], dec[3]))
    })
  })
  return(decisions)
}

#Get power from decisions
get_power <- function(dec.list){
  power <- lapply(dec.list, function(x){
    power <- colMeans(do.call(rbind, x))
    return(power)
  })
}

#Plot power curves
plot_power <- function(g.power, beta){
  power.df <- as.data.frame(cbind(unlist(beta), g.power))
  colnames(power.df) <- c("Beta", "C", "L", "Small-world coeff")
  power.df <- melt(power.df, id = "Beta")
  ggplot(data = power.df, 
         aes(x = Beta, y = value, color = variable)) + 
    geom_line(size = 1)
}

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

### Sample network from Chung-Lu
sample_from_CL <- function(P){
  n <- nrow(P)
  A <- matrix(0, nrow = n, ncol = n)
  A[col(A) > row(A)] <- runif(n*(n - 1)/2)
  A <- A + t(A)
  diag(A) <- runif(n)
  A <- (A < P) + 0
  G <- graph_from_adjacency_matrix(A, mode = "undirected")
  return(G)
}

get_P_SBM <- function(n, k, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- Z%*%p%*%t(Z)
  return(P)
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


get_P_DCSBM <- function(n, k, theta, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- diag(theta)%*%Z%*%p%*%t(Z)%*%diag(theta)
  return(P)
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
  
  t <- get_t_clust_1(G)
  t_1 <- get_t(G)
  L <- mean_distance(G)
  
mystats<-function(bt,P){
   Gstar = sample_from_CL(P)
   Lstar= cl_dist(Gstar)
   tstar = get_t_clust_1(Gstar)
   tLstar = transitivity(Gstar, type = "global")/cl_dist(Gstar)
    remove(Gstar)
   return(c(Lstar,tstar, tLstar))
}

 bootstats<-foreach(bt=1:b, .packages=c("igraph"),.export=ls(globalenv()),.combine='rbind') %dopar% mystats(bt,P)
  
  Lstar=bootstats[,1]
  tstar = bootstats[,2]
  tLstar = bootstats[,3]


  pval1 = sum(tstar>t)/b
  pval2 = sum(Lstar<L)/b
  pval3 = sum(tLstar>t_1)/b
  
  return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
              coef1=t, coef2=L, coef3=t_1,
              tstar=tstar,Lstar=Lstar, tLstar=tLstar))
}


### Compute p-value with ER Null model

small_world_ER <-function(G,b)
{
	n=gorder(G)
   L = mean_distance(G)
   C= transitivity(G, type = "global")
   #t=sqrt(n)*C/L
   t_1=get_t(G)
   t=get_t_clust_1(G)
   phat = gsize(G)/choose(n,2)


mystats<-function(bt,n,phat){
   Gstar = sample_gnp(n=n,p=phat)
   Lstar= mean_distance(Gstar)
   tstar = get_t_clust_1(Gstar)
   tLstar = transitivity(Gstar, type = "global")/mean_distance(Gstar)
    remove(Gstar)
   return(c(Lstar,tstar, tLstar))
}

 bootstats<-foreach(bt=1:b, .packages=c("igraph"),.export=ls(globalenv()),.combine='rbind') %dopar% mystats(bt,n,phat)
  
  Lstar=bootstats[,1]
  tstar = bootstats[,2]
  tLstar = bootstats[,3]



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
   t=get_t_clust_1(G)
   A=as_adj(G,type="both",sparse=FALSE)

ones = matrix(1,n,n)
diag(ones) = 0
   phat = matrix(0,k,k)
   for (i in 1:k)
   {for ( j in i:k)
   {
   phat[i,j] = sum(A[clusters==i,clusters==j])/sum(ones[clusters==i,clusters==j])
   phat[j,i]=phat[i,j]
   }
   }
   
 P <- get_P_SBM(n, k, phat, clusters)

mystats<-function(bt,n,P){
   Gstar = sample_from_SBM(n,P)
   Lstar= mean_distance(Gstar)
   tstar = get_t_clust_1(Gstar)
   tLstar = transitivity(Gstar, type = "global")/mean_distance(Gstar)
    remove(Gstar)
   return(c(Lstar,tstar, tLstar))
}

 bootstats<-foreach(bt=1:b, .packages=c("igraph"),.export=ls(globalenv()),.combine='rbind') %dopar% mystats(bt,n,P)
  
  Lstar=bootstats[,1]
  tstar = bootstats[,2]
  tLstar = bootstats[,3]

   pval1 = sum(tstar>t)/b
   pval2 = sum(Lstar<L)/b
   pval3 = sum(tLstar>t_1)/b
   
   return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
               coef1=t, coef2=L, coef3=t_1,
               tstar=tstar,Lstar=Lstar, tLstar=tLstar))
	
}





### Compute p value with DCSBM null model. 


small_world_DCSBM <-function(G,k,b, clusters)
{
    n=gorder(G)
    L = mean_distance(G)
    C= transitivity(G, type = "global")
    t_1=get_t(G)
    t = get_t_clust_1(G)
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
     
    P <- get_P_DCSBM(n, k, thetahat, phat, clusters)
    
    mystats<-function(bt,n,P){
   Gstar = sample_from_DCSBM(n, P)
   Lstar= mean_distance(Gstar)
   tstar = get_t_clust_1(Gstar)
   tLstar = transitivity(Gstar, type = "global")/mean_distance(Gstar)
    remove(Gstar)
   return(c(Lstar,tstar, tLstar))
}

 bootstats<-foreach(bt=1:b, .packages=c("igraph"),.export=ls(globalenv()),.combine='rbind') %dopar% mystats(bt,n,P)
  
  Lstar=bootstats[,1]
  tstar = bootstats[,2]
  tLstar = bootstats[,3]

    
    pval1 = sum(tstar>t)/b
    pval2 = sum(Lstar<L)/b
    pval3 = sum(tLstar>t_1)/b
    
    return(list(pval1=pval1,pval2=pval2, pval3=pval3, 
                coef1=t, coef2=L, coef3=t_1,
                tstar=tstar,Lstar=Lstar, tLstar=tLstar))
    
}



