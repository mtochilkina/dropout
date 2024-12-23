library(Rmagic)
library(Seurat)
library(Matrix)

args<-commandArgs(TRUE)
path<-args[1]
path_to_result<-args[2]
param_t<-as.integer(args[3])

dataset<-read.table(path, header=TRUE, row.names=1)
dataset<-CreateSeuratObject(counts=dataset, min.cells=3, min.features=200)
dataset<-NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
A_norm<-t(as.matrix(dataset[["RNA"]]$data)) 
result.completed<-magic(A_norm, genes="all_genes", t=param_t)
A_norm_completed<-result.completed$result
A_norm_completed<-t(A_norm_completed)
write.table(A_norm_completed, paste(path_to_result, '_magic_', as.character(param_t), '.tsv', sep=''), sep='\t')