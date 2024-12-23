library(devtools)
library("scImpute")
library(rsvd)
library(readr)
library(Seurat)
library(Matrix)

args<-commandArgs(TRUE)
path<-args[1] 
outdir<-args[2]
drop_prob<-args[3]
path_to_result<-args[4]

scimpute(# full path to raw count matrix
    count_path = (path), 
    infile = "txt",           # format of input file
    outfile = "txt",          # format of output file
    out_dir = outdir,           # full path to output directory
    labeled = FALSE,          # cell type labels not available
    drop_thre = drop_prob,     # threshold set on dropout probability
    Kcluster = 1,             # 2 cell subpopulations
    ncores = 20)       

A_completed<-read.table(paste(outdir, "scimpute_count.txt", sep='/'),  row.names=1)
A_completed<-as.matrix(A_completed)
A_completed<-CreateSeuratObject(counts=A_completed)
A_completed<-NormalizeData(A_completed, normalization.method = "LogNormalize", scale.factor = 10000)
write.table(as.matrix(A_completed[["RNA"]]$data), paste(path_to_result,"_scimpute_0.5.tsv", sep=""), sep = "\t")
  
