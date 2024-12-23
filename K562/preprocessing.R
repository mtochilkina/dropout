library(Seurat)
preprocessing<-function(path_to_data, data_name, path_to_resut){
  dataset<-Read10X(data.dir=path_to_data)
  dataset<-CreateSeuratObject(counts=dataset, min.cells=3, min.features=200)
  write.table(as.matrix(dataset[["RNA"]]$data), paste(path_to_result, data_name, "_std.txt", sep=""), sep = "\t" ) # knn, saver, scimpute
}

args = commandArgs(trailingOnly=TRUE)
path_to_data<-args[1]
data_name<-args[2]
path_to_result<-args[3]
preprocessing(path_to_data, data_name, path_to_resut)