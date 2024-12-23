library("splatter")
library("MASS")

trachea_droplet<-read.table("~/1tb/data/Montoro_data/Trachea_droplet.tsv", header=TRUE, row.names=1)
print("Trachea read")

#simulate basal type from Montoro data
basal<-trachea_droplet[, unlist(lapply(colnames(trachea_droplet), function(x) strsplit(x, "_")[[1]][1]))=="Basal"]
basal<-basal[, 100:3100]
basal<-basal[!(rowSums(basal)==0),]
basal<-as.matrix(basal)
params_basal <- zinbEstimate(basal) 
print("Estimated")
simulated.basal<- zinbSimulate(params = params_basal, sparsify = TRUE, verbose = TRUE)
print("Simulated")
write.matrix(assay(simulated.basal, "TrueCounts"), "zinb-wave_basal_true_counts.tsv", sep="\t")
write.matrix(assay(simulated.basal, "counts"), "zinb-wave_basal_obs_counts.tsv", sep="\t")
print("Basal done")

#simulate club type from Montoro data
club<-trachea_droplet[, unlist(lapply(colnames(trachea_droplet), function(x) strsplit(x, "_")[[1]][1]))=="Club"]
club<-club[, 50:2050]
club<-club[!(rowSums(club)==0),]
club<-as.matrix(club)
params_club <- zinbEstimate(club) 
print("Estimated")
simulated.club<- zinbSimulate(params = params_club, sparsify = TRUE, verbose = TRUE)
print("Simulated")
write.matrix(assay(simulated.club, "TrueCounts"), "zinb-wave_club_true_counts.tsv", sep="\t")
write.matrix(assay(simulated.club, "counts"), "zinb-wave_club_obs_counts.tsv", sep="\t")
print("Club done")

#simulate ciliated type from Montoro data
ciliated<-trachea_droplet[, unlist(lapply(colnames(trachea_droplet), function(x) strsplit(x, "_")[[1]][1]))=="Ciliated"]
ciliated<-ciliated[, 1:300]
ciliated<-ciliated[!(rowSums(ciliated)==0),]
ciliated<-as.matrix(ciliated)
params_ciliated <- zinbEstimate(ciliated) 
print("Estimated")
simulated.ciliated<- zinbSimulate(params = params_ciliated, sparsify = TRUE, verbose = TRUE)
print("Simulated")
write.matrix(assay(simulated.ciliated, "TrueCounts"), "zinb-wave_ciliated_true_counts.tsv", sep="\t")
write.matrix(assay(simulated.ciliated, "counts"), "zinb-wave_ciliated_obs_counts.tsv", sep="\t")
print("Ciliated done")
