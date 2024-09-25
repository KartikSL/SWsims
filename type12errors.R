
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





