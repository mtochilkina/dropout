library(rsvd)
library(readr)
library(Matrix)
library(Seurat)
source('alra_script.R')

args<-commandArgs(TRUE)

path<-args[1]    # path to count data 
path_to_result<-args[2]

dataset<-read.table(path, header=TRUE, row.names=1)
dataset<-CreateSeuratObject(counts=dataset, min.cells=3, min.features=200)
dataset<-NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
A_norm<-t(as.matrix(dataset[["RNA"]]$data)) 
result.completed <- alra(A_norm)
A_norm_completed <- result.completed[[3]]
rownames(A_norm_completed)<-rownames(A_norm)
A_norm_completed<-t(A_norm_completed)
write.table(A_norm_completed, paste(path_to_result, '_alra.tsv', sep=''), sep='\t')