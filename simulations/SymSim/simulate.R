library("splatter")
library("MASS")
library("SymSim")

#simulation basal - 2nd best param, club - 2nd best param, ciliated - 1st best param

cell_type<-list("basal", "club", "ciliated")

ngenes<- list(16088, 16056, 13762)
ncells_total<-list(3001, 2001, 300)
depth_mean=list(45000, 45000, 70000)
depth_sd=list(22500, 22500,  35000)
scale_s=list(0.2, 0.3,  0.1)
alpha_mean=list(0.01, 0.01, 0.03)
alpha_sd=list(0.005,0.005, 0.015)
mean_hge=list(4, 5, 4)
prop_hge=list(0.025, 0.02, 0.03)
nPCR1=list(14, 14, 14 )

print('Getting started...')
for (i in 1:length(ngenes)){
  true_counts_res <- SimulateTrueCounts(ncells_total=ncells_total[[i]], ngenes=ngenes[[i]], evf_type="one.population", 
                                        Sigma=0.2, randseed=0,
                                        gene_effect_prob=0.3, gene_effects_sd=2, 
                                        scale_s=scale_s[[i]], mean_hge=mean_hge[[i]],
                                        prop_hge=prop_hge[[i]])
  data(gene_len_pool)
  gene_len <- sample(gene_len_pool, ngenes[[i]], replace = FALSE)
  observed_counts <- True2ObservedCounts(true_counts=true_counts_res[[1]], 
                                         meta_cell=true_counts_res[[3]], protocol="UMI", 
                                         alpha_mean=alpha_mean[[i]], alpha_sd=alpha_sd[[i]], 
                                         gene_len=gene_len, 
                                         depth_mean=depth_mean[[i]], depth_sd=depth_sd[[i]],
                                         nPCR1=nPCR1[[i]])
  write.table(true_counts_res[[1]], paste(paste('symsim', cell_type[[i]], 'true_counts',  sep = '_'), '.txt', sep=''))
  write.table(observed_counts[[1]], paste(paste('symsim', cell_type[[i]], 'obs_counts', sep = '_'), '.txt', sep=''))
  print(as.character(i))
  print('Dataset created')
}

