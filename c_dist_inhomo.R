library(igraph)
library(ggplot2)
library(fastRG)


setwd("/Users/paul.963/Library/CloudStorage/OneDrive-TheOhioStateUniversity/small world simulation/clustcoefdist1")
nodes=1000
probs=c(0.004, 0.006, 0.008)
#probs=0.2
b=1000

get_P_SBM <- function(n, k, p, clusters){
  Z <- matrix(0, n, k)
  for(i in 1:n){
    Z[i, clusters[i]] <- 1
  }
  P <- Z%*%p%*%t(Z)
  return(P)
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



for ( j in 1:4){
  tstar =  rep(NA,b)

  n=nodes
  p=probs[j]
  k=4
  # use the sbm function to generate sbm and dcsbm function to generate dcsbm.
  #sbm_obj = fastRG::sbm(n=n,k=k,expected_density=p,poisson_edges = FALSE,allow_self_loops = FALSE)
  sbm_obj = fastRG::dcsbm(theta=rlnorm(n,1,0.5),k=k,expected_density=p,poisson_edges = FALSE,allow_self_loops = FALSE)
  P=2*expectation(sbm_obj)
  diag(P)=0
  expt = sum(diag((P%*%P%*%P)))/(6*choose(n,3))
  ones = matrix(1,n,n)
  exps = sum(diag(t(P%*%P)%*%ones - P%*%P))/(2*choose(n,3))
  P1 = P*(1-P)
  thetaT = sum(diag((P1%*%P1%*%P1)))/(6*choose(n,3))
  
  for (iter in 1:b){


    Gstar = fastRG::sample_igraph(sbm_obj)
    tstar[iter] = transitivity(Gstar, type = "global")
 
  }
  myfile <- file.path(paste0("DCSBM","n=",n, "p=",p*1000,".pdf"))

  
  
  Cstar = (exps/3)*sqrt(choose(n,3))*(tstar-3*expt/exps)/(sqrt(thetaT))
   mydf= data.frame(Cvalue=Cstar)

  myplot = ggplot(data=mydf)+
    geom_histogram(mapping = aes(x=Cvalue, y=after_stat(density)),fill="white", color="blue",alpha=1, position="identity")+
    stat_function(fun = dnorm)+
    theme(legend.position="bottom",legend.title=element_text(size=18,face="bold"),legend.text=element_text(size=18,face="bold")) +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18,face="bold"),
          plot.title=element_text(size=18,face="bold"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) + 
    ggtitle(paste0("n=",n, ", ", "p=",p))+xlab("Scaled clustering coefficient")+ylab("Density")+
    xlim(-5,5)

  ggsave(myfile,plot=myplot)
}




