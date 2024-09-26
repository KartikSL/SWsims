library(igraph)
library(ggplot2)


setwd("/Users/paul.963/Library/CloudStorage/OneDrive-TheOhioStateUniversity/small world simulation/clustcoefdist1")
nodes=1000
probs=c(0.004, 0.006, 0.008, 0.01, 0.02, 0.04, 0.06, 0.10, 0.20)
#probs=0.2
b=10000


  for ( j in 1:9){
    tstar =  rep(NA,b)
    #lstar = rep(NA,b)
    #swstar = rep(NA,b)
    #tcount = rep(NA,b)
    n=nodes
    p=probs[j]
    for (iter in 1:b){
      Gstar = sample_gnp(n=n,p=p)
      tstar[iter] = transitivity(Gstar, type = "global")

    }
    myfile <- file.path(paste0("compare","n=",n, "p=",p*1000,".pdf"))
    sigma11 = 1+((1+p-2*p^2)/(3*p^2*(n-2)))
    sigma22 = 1 + (1-p)/(4*p*(n-2))
    sigma12 = 1 + (1-p)/(2*p*(n-2))
    SigmaC = 9*sigma11 -12*sigma12+4*sigma22
    
    Cstar = sqrt(n*(n-1)/(2*p*(1-p)))*(tstar - p)
    
    Cstar1 = sqrt(n^3*p/(6*(1+2*p)*(1-p)^2))*(tstar-p)
    Cstar2 = sqrt(n*(n-1)/(2*p*(1-p)*SigmaC))*(tstar - p)
    xsq = x=seq(-3.5,3.5,by=0.001)
    sqcdist = dnorm(x, mean=0,sd=1)
    #sqcdist = rlnorm(100000, meanlog = lmean, sdlog = sqrt(lvar))
    
    mydf = data.frame(Cvalue = c(Cstar, Cstar1, Cstar2), scaling = c(rep("scaling 1",b), rep("scaling 2",b), rep("interp", b)))
    
    #pdf(myfile)
    #par(cex=1.4)
    myplot = ggplot(data=mydf)+
      geom_histogram(mapping = aes(x=Cvalue, y=after_stat(density), color=scaling),fill="white", alpha=0.6, position="identity")+
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




