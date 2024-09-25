
########## Contains all functions used except actual bootstrap simulations ############## 

setwd("~/Dropbox/small world detection SupSBM/Small World Bootstrap")
source("small_world_test.R")

#Define function to simulate WS graphs
simulate_ws <- function(p, d, s, n){
  return(sample_smallworld(d, s, n, p, loops = FALSE, 
                           multiple = FALSE))
}

#Define function to simulate Newman-Watts graphs ER/lattice
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

#Define function to simulate list of WS graphs
simulate_ws_list <- function(p, d, s, n){
  return(replicate(m, simulate_ws(p, d, s, n), simplify = FALSE))
}

#Define function to simulate list of NW graphs
# simulate_nw_list <- function(beta, n, d, m){
#   return(replicate(m, simulate_nw(n, d, beta), simplify = FALSE))
# }

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

#Small world coefficient using normalized C
get_t_mod <- function(g){
  L <- mean_distance(g)
  C <- get_t_clust_2(g)
  return(C/L)
}

#Range of x values for plotting
get_range <- function(df){
  min.t <- min(df)
  max.t <- max(df)
  range.t <- range(df)
  range.t <- range.t[2] - range.t[1]
  return(c(min.t - 0.15*range.t, max.t+ 0.15*range.t))
}


#Function ro simulate NW graphs SBM/lattice
simulate_nw_SBM <- function(n, d, beta){
  #Generate lattice
  l <- as_adj(sample_smallworld(1, n, ceiling(d*beta), 0))
  #Generate random edges
  k=2
  clusters = c(rep(1,n/2),rep(2,n/2))
  phat = matrix(0.5,k,k)
  diag(phat) = rep(1.5,k)
  p = 2*d*(1 - beta)/(n-1)*phat
  e <- as_adj(sample_from_SBM(n,k, p,clusters))
  #Superimpose edges on lattice
  return(graph_from_adjacency_matrix(l + e, "undirected",weighted=TRUE))
}

#Covariance matrix for 
get_sigma <- function(n, p){
  m <- c(4*p^2 + p*(1 - p)/(n - 2), 
         2*p^3 + (p^2)*(1 - p)/(n - 2), 
         2*p^3 + (p^2)*(1 - p)/(n - 2),
         p^4 + (p^2)*(1 + p - 2*p^2)/(3*n - 6))
  return(3*(n - 2)*choose(n, 3)*p*(1 - p)*matrix(m, nrow = 2))
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

#Get data into format for ggplot
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

refactor_data_2 <- function(data.dcsbm.1, data.dcsbm.2){
  
  coef <- c(data.dcsbm.1$coef1, data.dcsbm.1$coef2, data.dcsbm.1$coef3)
  pval <- list(c(data.dcsbm.1$pval1, data.dcsbm.1$pval2, data.dcsbm.1$pval3),
               c(data.dcsbm.2$pval1, data.dcsbm.2$pval2, data.dcsbm.2$pval3))
  
  spectral <- rep("Spectral", length(data.dcsbm.1$tstar))
  score <- rep("SCORE", length(data.dcsbm.2$tstar))
  
  data.spectral <- data.frame(clustering = spectral, tstar = data.dcsbm.1$tstar,
                        lstar = data.dcsbm.1$Lstar, tlstar = data.dcsbm.1$tLstar)
  data.score <- data.frame(clustering = score, tstar = data.dcsbm.2$tstar,
                         lstar = data.dcsbm.2$Lstar, tlstar = data.dcsbm.2$tLstar)
  
  df <- rbind(data.spectral, data.score)
  return(list(df, coef, pval))
}

plot_real <- function(real_data){
  library(ggplot2)
  library(patchwork)
  
  df <- real_data[[1]]
  coef <- real_data[[2]]
  pval <- real_data[[3]]
  colnames(df)[1] <- "null model"
  
  p.1 <- ggplot(df, aes(x = tstar, color = `null model`, fill = `null model`)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = coef[1], color = "red") + 
    xlab("C") + 
    ylab("") + 
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size=16), 
          legend.text = element_text(size = 14))
  
  p.2 <- ggplot(df, aes(x = lstar, color = `null model`, fill = `null model`)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = coef[2], color = "red") + 
    xlab("L") + 
    ylab("") + 
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size=16), 
          legend.text = element_text(size = 14))
  
  p.3 <- ggplot(df, aes(x = tlstar, color = `null model`, fill = `null model`)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = coef[3], color = "red") + 
    xlab("C/L") + 
    ylab("") + 
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size=16), 
          legend.text = element_text(size = 14))
  
  # g_legend<-function(a.gplot){
  #   tmp <- ggplot_gtable(ggplot_build(a.gplot))
  #   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  #   legend <- tmp$grobs[[leg]]
  #   return(legend)}
  # 
  # mylegend<-g_legend(p.3)
  
  # p <- list(plot_grid(p.1, p.2, ncol = 2), p.3)
  # p <- grid.arrange(arrangeGrob(p.3 + theme(legend.position="none"),
  #                               p.1 + theme(legend.position="none"),
  #                                p.2 + theme(legend.position="none"),
  #                                nrow=1),
  #                    mylegend, nrow=1,heights=c(10, 1))
  p <- p.3 + p.1 + p.2 & theme(legend.position = "right")
  p <- p + plot_layout(guides = "collect")
  
  return(p)
}

#Plotting function for real world data
plot_real_2 <- function(real_data){
  library(ggplot2)
  library(gridExtra)
  
  df <- real_data[[1]]
  colnames(df)[1] <- "method"
  df$clustering <- factor(df$method, levels = unique(df$method))
  coef <- real_data[[2]]
  pval <- real_data[[3]]
  
  p.1 <- ggplot(df, aes(x = tstar, color = method, fill = method)) + 
    geom_density(alpha = 0.4) + 
    # geom_vline(xintercept = coef[1], color = "red") + 
    xlab("C") + 
    ylab("") + 
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size=16), 
          legend.text = element_text(size = 14))
  
  p.2 <- ggplot(df, aes(x = lstar, color = method, fill = method)) + 
    geom_density(alpha = 0.4) + 
    # geom_vline(xintercept = coef[2], color = "red") + 
    xlab("L") + 
    ylab("") + 
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size=16), 
          legend.text = element_text(size = 14))
  
  p.3 <- ggplot(df, aes(x = tlstar, color = method, fill = method)) + 
    geom_density(alpha = 0.4) + 
    # geom_vline(xintercept = coef[3], color = "red") + 
    xlab("C/L") + 
    ylab("") + 
    theme(legend.position = "right", 
          axis.title.x = element_text(size = 18),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size=16), 
          legend.text = element_text(size = 14))
  
  p <- p.3 + p.1 + p.2 & theme(legend.position = "right")
  p <- p + plot_layout(guides = "collect")
  
  return(p)
}

get_var_log_c <- function(nw.list){
  sapply(nw.list, function(x){
    nw.c <- sapply(x, transitivity)
    return(var(log(nw.c)))
  })
}

get_mean_log_c <- function(nw.list){
  sapply(nw.list, function(x){
    nw.c <- sapply(x, transitivity)
    return(mean(nw.c))
  })
}

emp_mean_tri_nw <- function(nw.list){
  sapply(nw.list, function(x){
    mean.tri <- mean(sapply(x, function(x){
      # length(cliques(x, min = 3, max = 3))
      sum(count_triangles(x))/3
    }))
  })
}

t_mean_tri_nw <- function(beta, delta, n){
  #triangles in ring lattice
  t.1 <- n*choose(delta*beta, 2) 
  #triangles in random graph
  t.2 <- choose(n, 3)*((2*delta*(1 - beta))/(n - 1))^3
  #triangles through intersection of random edge and two lattice edges
  # t.3 <- ((2*delta*(1 - beta))/(n - 1))*choose(n, 2)*
  #   (delta*beta/n*(n - 1))*(delta*beta*(delta*beta + 1)/2)
  # return(t.1 + t.2 + t.3)
  return(t.1 + t.2)
}

emp_mean_trip_nw <- function(nw.list){
  sapply(nw.list, function(x){
    mean.trip <- mean(sapply(x, function(x){
      A.sq <- as_adjacency_matrix(x)^2
      A.sq[lower.tri(A.sq)]  <- 0
      sum(A.sq)
    }))
  })
}

t_mean_trip_nw <- function(beta, delta, n){
  #triples in lattice
  t.1 <- n*choose(2*delta*beta, 2)
  #triples in random graph
  t.2 <- choose(n, 3)*(((2*delta*(1 - beta))/(n - 1))^2)*(1 - (2*delta*(1 - beta))/(n - 1))
  return(t.1 + t.2)
}

var_log_c <- function(beta, delta, n){
  p <- 2*delta*(1 - beta)/(n - 1)
  sigma_nw <- get_sigma(n, p)
  mean_t <- choose(n, 3)*p^3 #number of triangles
  mean_s <- 3*choose(n, 3)*p^2 - 3*mean_t #number of triples not triangles
  return(sigma_nw[2, 2]/(9*mean_t^2) + 
           sigma_nw[1, 1]/(mean_s + 3*mean_t)^2 - 
           2*sigma_nw[1, 2]/(3*mean_t*(mean_s + 3*mean_t)))
}

mean_log_c <- function(beta, delta){
  c.mean <- 3*(delta - 1)/
    (2*(2*delta - 1) + 
       8*delta*(1 - beta) + 
       4*delta*(1 - beta)^2)
  return(c.mean)
}

#Get decisions from results object
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
    geom_line(position=position_jitter(w=0.002, h=0.001), size = 1)
}

#Theorem 3
# prob_l <- sapply(nw.list.1000, function(x){
#   mean(
#     sapply(x, function(y){
#       phat <- gsize(y)/choose(n,2)
#       K <- (2 + epsilon)*log(n)/log(n*phat)
#       return(mean_distance(y) > K)
#     })
#   )
# })

get_decisions_osc <- function(result){
  seq <- c(1:3)
  dec <- vector(length = 3)
  dec.vec <- vector(mode = "list", length = 300)
  for(i in 1:300){
    p.val <- unlist(result[seq])
    thresh <- c(0.05, 0.99, 0.05)
    dec <- p.val < thresh
    dec.vec[[i]] <- dec
    seq <- seq + 9
  }
  return(do.call(rbind, dec.vec))
}

get_statistics_osc <- function(result){
  seq <- c(7:9)
  tstar <- vector(mode = "list", length = 300)
  lstar <- vector(mode = "list", length = 300)
  tlstar <- vector(mode = "list", length = 300)
  for(i in 1:300){
    tstar[[i]] <- result[seq[1]]
    lstar[[i]] <- result[seq[2]]
    tlstar[[i]] <- result[seq[3]]
    seq <- seq + 9
  }
  g.stat <- list(unlist(tstar), unlist(lstar), unlist(tlstar))
  return(g.stat)
}

summarise_g_stat <- function(g.stat, beta){
  t.stat <- sapply(g.stat, function(x){
    list.t <- x[[1]]
    t.mean <- mean(list.t)
    t.sort <- sort(list.t)
    t.upper <- t.sort[0.95*length(t.sort)]
    t.lower <- t.sort[0.05*length(t.sort)]
    return(c(t.mean, t.lower, t.upper))
  })
  t.stat <- as.data.frame(t(t.stat))
  t.stat$stat <- rep("C", dim(t.stat)[1])
  t.stat$beta <- unlist(beta)
  colnames(t.stat) <- c("mean", "lower", "upper", "stat", "beta")
  
  l.stat <- sapply(g.stat, function(x){
    list.l <- x[[2]]
    l.mean <- mean(list.l)
    l.sort <- sort(list.l)
    l.upper <- l.sort[0.95*length(l.sort)]
    l.lower <- l.sort[0.05*length(l.sort)]
    return(c(l.mean, l.lower, l.upper))
  })
  l.stat <- as.data.frame(t(l.stat))
  l.stat$stat <- rep("L", dim(l.stat)[1])
  l.stat$beta <- unlist(beta)
  colnames(l.stat) <- c("mean", "lower", "upper", "stat", "beta")
  
  return(rbind(t.stat, l.stat))
  
}

plot_stat <- function(g.stat){
  library(cowplot)
  
  g.stat.c <- subset(g.stat, stat == "C")
  g.stat.l <- subset(g.stat, stat == "L")
  pd <- position_dodge(0.1)
  
  plot.c <- ggplot(g.stat.c, aes(x = beta, y = mean)) + 
    geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
    geom_line() + 
    geom_point()
  
  plot.l <- ggplot(g.stat.l, aes(x = beta, y = mean)) + 
    geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
    geom_line() + 
    geom_point()
  
  plot_grid(plot.c, plot.l, ncol = 2)
}

get_power_osc <- function(dec){
  return(colMeans(dec))
}

get_summary_nw <- function(stat.list, beta.list){
  stat.summary <- sapply(stat.list, function(x){
    stat.mean <- median(x)
    stat.upper <- sort(x)[0.99*length(x)]
    stat.lower <- sort(x)[0.01*length(x)]
    return(c(stat.mean, stat.lower, stat.upper))
  })
  stat.summary <- as.data.frame(t(stat.summary))
  colnames(stat.summary) <- c("mean", "lower", "upper")
  stat.summary$beta <- unlist(beta.list)
  return(stat.summary)
}

compute_cl <- function(nw.list, beta.nw.list){
  nw.c <- lapply(nw.list, function(x){
    sapply(x, function(y){
      transitivity(y)
    })
  })
  nw.l <- lapply(nw.list, function(x){
    sapply(x, function(y){
      mean_distance(y)
    })
  })
  nw.c <- get_summary_nw(nw.c, beta.nw.list)
  nw.c$stat <- rep("C", length(beta.nw.list))
  nw.l <- get_summary_nw(nw.l, beta.nw.list)
  nw.l$stat <- rep("L", length(beta.nw.list))
  return(rbind(nw.c, nw.l))
}

plot_stat_nw <- function(nw.list, beta.nw.list){
  g.stat <- compute_cl(nw.list, beta.nw.list)
  ggplot(g.stat, aes(x = beta, y = mean, colour = stat)) + 
    geom_errorbar(aes(min = lower, max = upper, width = 0.1)) + 
    geom_line(size = 2) + 
    geom_point(size = 2)
}
