library(SAVER)
library(rsvd)
library(readr)
library(Matrix)
library(Seurat)

args<-commandArgs(TRUE)
path<-args[1] 
path_to_result<-args[2]

dataset<-read.table(path, header=TRUE, row.names=1)
dataset<-CreateSeuratObject(counts=dataset, min.cells=3, min.features=200)
A<-dataset[["RNA"]]$counts
A_completed <- saver(A, ncores=28, estimates.only = TRUE, size.factor = 1)
A_completed<- CreateSeuratObject(counts=A_completed)
A_completed<-NormalizeData(A_completed, normalization.method = "LogNormalize", scale.factor = 10000)
write.table(A_completed[["RNA"]]$data, paste(path_to_result,"_saver.tsv", sep=""), sep = "\t")