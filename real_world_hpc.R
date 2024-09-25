detachDoParallel <- function() {
  detach("package:doParallel")
  detach("package:foreach")
  detach("package:parallel")
  detach("package:iterators")
}

library(doParallel)
library(Rmpi)

library(igraph)
library(ggplot2)
library(cowplot)

slaves <- detectCores() - 1
#slaves <- as.numeric(Sys.getenv(c("PBS_NP")))-1
{ sink("/dev/null"); cl <- makeCluster(slaves, type="MPI"); sink(); } # number of MPI tasks to use
registerDoParallel(cl)

source("small_world_test.R")

real_world_hpc <- function(g, k, g.name){
  tick<- proc.time()
  cat(tick, ": Started execution\n")
  
  g.er <- small_world_ER(g,b=500)
  tock <- proc.time() - tick
  cat(tock, ": ER done, started evaluating clusters\n")
  
  clusters.1 <- dcspectral(as_adj(g, type="both", sparse=FALSE), gorder(g), k)
  clusters.2 <- SCORE(as_adjacency_matrix(g, sparse = F),K=k)
  
  tock <- proc.time() - tick
  cat(tock, ": Clusters evaluated\n")
  
  g.sbm <- small_world_SBM(g, k,b=500, clusters.1)
  tock <- proc.time() - tick
  cat(tock, ": SBM done\n")
  
  g.dcsbm.1 <- small_world_DCSBM(g, k,b=500, clusters.1)
  g.dcsbm.2 <- small_world_DCSBM(g, k,b=500, clusters.2)
  tock <- proc.time() - tick
  cat(tock, ": DCSBM done\n")
  
  g.cl <- small_world_CL(g,b=500)
  tock <- proc.time() - tick
  cat(tock, ": CL done\n")
  
  cat("\nforeach w/ Rmpi test times using", slaves, "MPI slaves: \n")
  cat(tock, "\n")
  
  g.data <- refactor_data(g.er, g.sbm, g.dcsbm.1, g.cl)
  g.plot <- plot_real(g.data)
  ggsave("Images/wordadj_final.png", g.plot, height = 3.5, width = 11, units = "in")
  
  g.data <- refactor_data_2(g.dcsbm.1, g.dcsbm.2)
  g.plot <- plot_real_2(g.data)
  ggsave("Images/dcbm_wordadj.png", g.plot, height = 3.5, width = 11, units = "in")
  
  save(g.data,file=paste(g.name, ".plot.RData", sep = ""))
  
  #pdf(paste(g.name, ".plot1.pdf", sep = ""), width=8,height=6)
  #g.plot[[1]]
  #dev.off()
  #pdf(paste(g.name, ".plot2.pdf", sep = ""), width=8,height=6)
  #g.plot[[2]]
  #dev.off()
}

#karate test
g <- read_graph("../Data/dolphins.gml", format = "gml")
k <- length(cluster_louvain(g))
g.name <- deparse(substitute(g))
real_world_hpc(karate, k, g.name)

source("SCOREplus.R")

# scorek<-SCORE(as_adjacency_matrix(g, sparse = F),K=k)
# simmonsscore = nmi(scorek,as.factor(simmons$labels))
# 
# scoreplusk<-SCOREplus(simmons$adj,k=k)
# simmonsscoreplus = nmi(scoreplusk$labels,as.factor(simmons$labels))


invisible(stopCluster(cl))
detachDoParallel()
mpi.quit()