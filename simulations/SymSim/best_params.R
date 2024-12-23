library("SymSim")
#read table
trachea_droplet<-read.table("~/1tb/data/Montoro_data/Trachea_droplet.tsv", header=TRUE, row.names=1)
print("Trachea read")

#calculate SymSim params for basal type in Montoro data
basal<-trachea_droplet[, unlist(lapply(colnames(trachea_droplet), function(x) strsplit(x, "_")[[1]][1]))=="Basal"]
basal<-basal[, 100:3100]
basal<-basal[!(rowSums(basal)==0),]
best_matches_UMI <- BestMatchParams('UMI',basal,'basal_best_params.umi.qqplot', depth_range = c(20e3,500e3),n_optimal=5)
write.table(best_matches_UMI, "symsim_basal_param.txt")


#calculate SymSim params for club type in Montoro data
club<-trachea_droplet[, unlist(lapply(colnames(trachea_droplet), function(x) strsplit(x, "_")[[1]][1]))=="Club"]
club<-club[, 50:2050]
club<-club[!(rowSums(club)==0),]
best_matches_UMI <- BestMatchParams('UMI',club,'club_best_params.umi.qqplot', depth_range = c(20e3,500e3),n_optimal=5)
write.table(best_matches_UMI, "symsim_club_param.txt")


#calculate SymSim params for ciliated type in Montoro data
ciliated<-trachea_droplet[, unlist(lapply(colnames(trachea_droplet), function(x) strsplit(x, "_")[[1]][1]))=="Ciliated"]
ciliated<-ciliated[, 1:300]
ciliated<-ciliated[!(rowSums(ciliated)==0),]
best_matches_UMI <- BestMatchParams('UMI',ciliated,'ciliated_best_params.umi.qqplot', depth_range = c(20e3,500e3),n_optimal=5)
write.table(best_matches_UMI, "symsim_ciliated_param.txt")



