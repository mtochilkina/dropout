library(reticulate)
library(anndata)
library(sceasy)
library(Seurat)

args<-commandArgs(TRUE)
path<-args[1] 
path_to_result<-args[2]

A_counts<-read.table(path, row.names=1)
A_counts<-CreateSeuratObject(counts=A_counts, min.cells=3, min.features=200)
A_counts[["RNA"]] <- as(A_counts[["RNA"]], "Assay")
adata <- convertFormat(A_counts, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
sc <- import('scanpy', convert = FALSE)
scvi<-import("scvi",  convert=FALSE)
scvi$model$SCVI$setup_anndata(adata)
model = scvi$model$SCVI(adata)
model$train()
denoised = model$get_normalized_expression(adata)  
A_completed<-py_to_r(denoised)
A_completed<-t(A_completed)
A_completed<-CreateSeuratObject(counts=A_completed)
A_completed<-NormalizeData(A_completed)
write.table(A_completed[["RNA"]]$data, paste(path_to_result,"_scvi.tsv",  sep=""), sep = "\t")