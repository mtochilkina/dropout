library(Seurat)
args<-commandArgs(TRUE)
path<-args[1] 
path_to_result<-args[2]
method_name<-args[3]

data<-read.table(path, row.names=1)

#data<-as.matrix(data)

postprocessing<-function(data, path_to_result, method_name){
  print("Starting postprocessing....")
  dataset<-CreateSeuratObject(counts=data)
  dataset<-NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)
  write.table(as.matrix(dataset[["RNA"]]$data), paste(path_to_result,"_", method_name, ".tsv", sep=""), sep = "\t")
  print("Succesfully writen")
}

postprocessing(data, path_to_result, method_name)