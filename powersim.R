
beta.nw.list <- as.list(seq(0, 1, 0.25))

nw.list.er.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                              d = 20, n = 400, m = 150, null = "ER")
nw.list.sbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                               d = 20, n = 400, m = 150, null = "SBM")
# nw.list.dcsbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
#                                d = 20, n = 500, m = 250, null = "DCSBM")

nw.list.dcsbm.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                                 d = 20, n = 400, m = 150, null = "DCSBM")

nw.list.dcsbm.sparse.2.sdh <- lapply(beta.nw.list, simulate_nw_list, 
                                     d = 20, n = 400, m = 150, null = "DCSBMSDH")

nw.list.cl.sparse.2 <- lapply(beta.nw.list, simulate_nw_list, 
                              d = 20, n = 500, m = 150, null = "CL")

for(beta in beta.nw.list){
  
}
