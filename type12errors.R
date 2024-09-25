
##### Power curves #####

detachDoParallel <- function() {
  detach("package:doParallel")
  detach("package:foreach")
  detach("package:parallel")
  detach("package:iterators")
}

suppressMessages(library(doParallel))
suppressMessages(library(Rmpi))

suppressMessages(library(igraph))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))


#slaves <- detectCores() - 1
slaves <- 150
#slaves <- as.numeric(Sys.getenv(c("PBS_NP")))-1
{ sink("/dev/null"); cl <- makeCluster(slaves, type="MPI"); sink(); } # number of MPI tasks to use
registerDoParallel(cl)

#source("small_world_test_hpc.R")
source("small_world_test.R")

beta.nw.list <- as.list(seq(0, 1, 0.1))


nw.list.er.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                            d = 20, n = 500, m = 150, null = "ER")
nw.list.sbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                             d = 20, n = 500, m = 150, null = "SBM")
# nw.list.dcsbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
#                                d = 20, n = 500, m = 250, null = "DCSBM")

nw.list.dcsbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                                 d = 20, n = 500, m = 150, null = "DCSBM")

nw.list.dcsbm.sparse.2.sdh <- lapply(beta.nw.list, simulate_nw_list,
                                 d = 20, n = 500, m = 150, null = "DCSBMSDH")

nw.list.cl.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                            d = 20, n = 500, m = 150, null = "CL")

# nw.list.er.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
#                             d =40, n = 1000, m = 250, null = "ER")
# nw.list.sbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
#                              d = 40, n = 1000, m = 250, null = "SBM")
# nw.list.dcsbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
#                                d = 40, n = 1000, m = 250, null = "DCSBM")
# nw.list.cl.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
#                             d = 40, n = 1000, m = 250, null = "CL")

#er.list <- lapply(nw.list.1000, lapply, small_world_ER, b = 200)
#er.dec <- get_decisions(er.list)
#er.power <- do.call(rbind, get_power(er.dec))
#plot_power(er.power, beta.nw.list)
# 
# cl.list <- lapply(nw.list.1000, lapply, small_world_CL)
# cl.dec <- get_decisions(cl.list)
# cl.power <- do.call(rbind, get_power(cl.dec))

#tick <- proc.time()
#cat("\n\nFinding number of clusters: ", tick, "\n")
#k.vec <- lapply(nw.list.1000, lapply, function(x){length(cluster_louvain(x))})
#tock <- proc.time() - tick
#cat("Number of clusters found. Time taken: ", tock, "\n")

#tick <- proc.time()
#cat("Computing clusters: ", tick, "\n")
#bm.clusters <- lapply(nw.list.1000, function(x){
#  lapply(x, function(y){
#    k <- length(cluster_louvain(y))
#    dcspectral(as_adj(y, type="both", sparse=FALSE), gorder(y), k)
#  })
#})
#tock <- proc.time() - tick
#cat("Clusters computed. Time taken: ", tock, "\n")

#sbm.list <- mapply(function(x, y, z){
#  mapply(function(g, k, cluster){
#    small_world_SBM(G = g, k = k, b = 200, clusters = cluster)
#  }, x, y, z)
#}, nw.list.1000, k.vec, bm.clusters)

# result <- foreach(i=1:11, .packages=c("igraph"), 
#                   .export=ls(globalenv()), 
#                   .combine='rbind') %dopar% {
#                     tick <- proc.time()
#                     cat(i, " iteration, start time: " ,tick, "\n")
#                     mapply(function(g, k, cluster){
#                       small_world_SBM(G = g, k = k, b = 200, clusters = cluster)
#                     }, nw.list.1000[[i]], k.vec[[i]], bm.clusters[[i]])
#                     tock <- proc.time() - tick
#                     cat(i, " iteration, total time: " ,tock, "\n")
# }

#Serial for each beta, parallel for network, bootstrap serial

result <- vector(mode = "list", length = 11)

cat("\nStarting execution\n")

######## ER #########

#2.5 hours for 21 iterations

# for(i in 1:11){
#    tick <- proc.time()
#    result[[i]] <- foreach(j = 1:200, .packages = c("igraph"), 
#                          .export = ls(globalenv()), .combine = 'c') %dopar% {
#                            g <- nw.list.1000[[i]][[j]]
#                            small_world_ER(G = g, b = 250)
#                          }
#    tock <- proc.time() - tick
#    cat("Iteration ", i, " done. Time taken: ", tock, "\n")
# }

######### SBM #######

#8 hours for 21 iterations

 for(i in 1:11){
    tick <- proc.time()
   result[[i]] <- foreach(j = 1:200, .packages = c("igraph"), 
                     .export = ls(globalenv()), .combine = 'c') %dopar% {
                       g <- nw.list.1000[[i]][[j]]
                       k <- length(cluster_louvain(g))
                       bm.cluster <- dcspectral(as_adj(g, type="both", sparse=FALSE), gorder(g), k)
                       small_world_SBM(G = g, 
                                       k = k, 
                                       b = 250,
                                       clusters = bm.cluster)
                     }
    tock <- proc.time() - tick
    cat("Iteration ", i, " done. Time taken: ", tock, "\n")
 }

######### DCSBM #######

#8 hours for 11 iterations

#for(i in 1:11){
#    tick <- proc.time()
#    result[[i]] <- foreach(j = 1:300, .packages = c("igraph"), 
#                         .export = ls(globalenv()), .combine = 'c') %dopar% {
#                           g <- nw.list.1000[[i]][[j]]
#                           k <- length(cluster_louvain(g))
#                           bm.cluster <- dcspectral(as_adj(g, type="both", sparse=FALSE), gorder(g), k)
#                           small_world_DCSBM(G = g, 
#                                           k = k, 
#                                           b = 300,
#                                           clusters = bm.cluster)
#                         }
#    tock <- proc.time() - tick
#    cat("Iteration ", i, " done. Time taken: ", tock, "\n")
#}

######## CL #########

#for(i in 1:11){
#    tick <- proc.time()
#    result[[i]] <- foreach(j = 1:300, .packages = c("igraph"), 
#                      .export = ls(globalenv()), .combine = 'c') %dopar% {
#                        g <- nw.list.1000[[i]][[j]]
#                        small_world_CL(G = g, b = 300)
#                          }
#    tock <- proc.time() - tick
#    cat("Iteration ", i, " done. Time taken: ", tock, "\n")
#}

save(result, file = "sbm_1000_200.RData")

#cat("Starting power computation\n")

#for(i in 1:11){
  #tick <- proc.time()
  #cat(i, " iteration, start time: " ,tick, "\n")
#  mapply(function(g, k, cluster){
#    small_world_SBM(G = g, k = k, b = 200, clusters = cluster)
#  }, nw.list.1000[[1]], k.vec[[1]], bm.clusters[[1]])
  #tock <- proc.time() - tick
  #cat(i, " iteration, total time: " ,tock, "\n")
#}

######Theorem 3

# delta <- 10
# n <- 500
# epsilon <- 0.1
# 
# l_df <- as.data.frame(cbind(unlist(beta.nw.list), prob_l))
# colnames(l_df) <- c("Beta", "P")
# 
# ggplot(l_df, aes(x = Beta, y = P)) + geom_point()

invisible(stopCluster(cl))
detachDoParallel()
mpi.quit()


######### Power ##########

nw.er.500.sparse.stat <- compute_cl(nw.list.er.sparse.2, beta.nw.list)
nw.er.500.sparse.stat$null <- rep("ER", dim(nw.er.500.sparse.stat)[1])

nw.sbm.500.sparse.stat <- compute_cl(nw.list.sbm.sparse.2, beta.nw.list)
nw.sbm.500.sparse.stat$null <- rep("SBM", dim(nw.sbm.500.sparse.stat)[1])

nw.dcsbm.500.sparse.stat <- compute_cl(nw.list.dcsbm.sparse.2, beta.nw.list)
nw.dcsbm.500.sparse.stat$null <- rep("DCSBM", dim(nw.dcsbm.500.sparse.stat)[1])

nw.cl.500.sparse.stat <- compute_cl(nw.list.cl.sparse.2, beta.nw.list)
nw.cl.500.sparse.stat$null <- rep("CL", dim(nw.cl.500.sparse.stat)[1])

# nw.dcsbm.400.sparse.stat <- compute_cl(nw.list.dcsbm.sparse.2, beta.nw.list)
# nw.dcsbm.400.sparse.stat$null <- rep("DCSBM", dim(nw.dcsbm.400.sparse.stat)[1])
# 
nw.dcsbmsdh.500.sparse.stat <- compute_cl(nw.list.dcsbm.sparse.2.sdh, beta.nw.list)
nw.dcsbmsdh.500.sparse.stat$null <- rep("DCSBMSDH", dim(nw.dcsbmsdh.500.sparse.stat)[1])

nw.500.sparse.df <- rbind(nw.er.500.sparse.stat, nw.sbm.500.sparse.stat, 
                          nw.dcsbm.500.sparse.stat, nw.cl.500.sparse.stat)

nw.500.dcsbm.df <- rbind(nw.dcsbm.500.sparse.stat, nw.dcsbmsdh.500.sparse.stat)

# generate regular dcsbm using U(0, 1) and compare
# VERIFY if N=1000 in sims

p.c <- ggplot(subset(nw.500.sparse.df, stat == "C"), 
       aes(x = beta, y = mean, colour = null)) + 
  geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
  geom_line() + 
  geom_point() + 
  ggtitle("Observed C") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = "bottom", 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size = 15)) + 
  labs(y = "C", x = expression(beta))

ggsave("Images/CvsBeta.png", p.c, width = 5, height = 5, units = "in")

p.l2 <- ggplot(subset(nw.500.sparse.df, stat == "L" & beta != 1), 
       aes(x = beta, y = mean, colour = null)) + 
  geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
  geom_line() + 
  geom_point() + 
  ggtitle(paste("Observed L, \u03b2", "=1 excluded")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = "bottom", 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size = 15)) + 
  labs(y = "L", x = expression(beta)) + 
  xlim(0, 0.9)

ggsave("Images/LvsBeta2.png", p.l2, width = 5, height = 5, units = "in")

p.l1 <- ggplot(subset(nw.500.sparse.df, stat == "L"), 
       aes(x = beta, y = mean, colour = null)) + 
  geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
  geom_line() + 
  geom_point() + 
  ggtitle("Observed L") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = "bottom", 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size = 15)) + 
  labs(y = "L", x = expression(beta))

ggsave("Images/LvsBeta1.png", p.l1, width = 5, height = 5, units = "in")



############## DCSBM #############

p.c.dcbm <- ggplot(subset(nw.500.dcsbm.df, stat == "C"), 
              aes(x = beta, y = mean, colour = null)) + 
  geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
  geom_line() + 
  geom_point() + 
  ggtitle("Observed C") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = "bottom", 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size = 15)) + 
  labs(y = "C", x = expression(beta))

ggsave("Images/CvsBetaDCBM.png", p.c.dcbm, width = 5, height = 5, units = "in")

p.l2.dcbm <- ggplot(subset(nw.500.dcsbm.df, stat == "L" & beta != 1), 
               aes(x = beta, y = mean, colour = null)) + 
  geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
  geom_line() + 
  geom_point() + 
  ggtitle(paste("Observed L, \u03b2", "=1 excluded")) +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = "bottom", 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size = 15)) + 
  labs(y = "L", x = expression(beta)) + 
  xlim(0, 0.9)

ggsave("Images/LvsBeta2DCBM.png", p.l2.dcbm, width = 5, height = 5, units = "in")

p.l1.dcbm <- ggplot(subset(nw.500.dcsbm.df, stat == "L"), 
               aes(x = beta, y = mean, colour = null)) + 
  geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
  geom_line() + 
  geom_point() + 
  ggtitle("Observed L") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), legend.position = "bottom", 
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.title.y = element_text(size = 20), 
        axis.text.y = element_text(size = 18), 
        legend.title = element_text(size=18), 
        legend.text = element_text(size = 15)) + 
  labs(y = "L", x = expression(beta))

ggsave("Images/LvsBeta1DCBM.png", p.l1.dcbm, width = 5, height = 5, units = "in")





