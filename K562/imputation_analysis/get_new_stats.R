library(dplyr)
library(tidyr)
library(parallel)
library(stringr)

count_expressing_cells<-function(path, sample, method_name){        #method_name
  data<-read.table(path,row.names=1, header=TRUE)
  rownames(data)<-gsub("-", ".", rownames(data))
  results<-apply(data, 1, function(data_row) sum(data_row>0)/length(data_row)) # calculate proportion of expressing cells for each gene
  results<-as.data.frame(results)
  colnames(results)<-c(paste("sample", sample,"_", method_name ,sep=""))
  return(results)
}
count_expression_level<-function(path, sample){ #aggregate counts and then log them
  data<-read.table(path,row.names=1, header=TRUE)
  results<-rowSums(data)
  results<-log1p(results)
  results<-as.data.frame(results)
  colnames(results)<-c(paste("sample", sub("\\.", "",as.character(sample)),sep=""))
  return(results)
}

param<-c(1,2,3,4,5, 6,7,8, 9, 10)  
magic_names<-lapply(param, function(x) paste("magic_",as.character(x), sep=""))
magic_names<-unlist(magic_names)

param<-c("2", "3", "4", "5", "6", "7", "8", "9", "10")
knn_names<-lapply(param, function(x) paste("knn_",as.character(x), sep=""))
knn_names<-unlist(knn_names)

param<-c("0.3", "0.5", "0.9")
scimpute_names<-lapply(param, function(x) paste("scimpute_",as.character(x), sep=""))
scimpute_names<-unlist(knn_names)

method_name<-c("alra", "scVI", "saver", magic_names, knn_names, scimpute_names)

sample=c("002") #, "004", "006", "008", "01", "012", "014", "016","018","02", "03", "04", "05", "06", "08", "1")  

#!!!
method_name<-c("knn_10","knn_11","knn_12","knn_13", "knn_14","knn_15")  #parApply could overwhelm memory, it is recommended to run script separately for group of files

combination<-tidyr::expand_grid(sample, method_name)

clust <- makeCluster(20)
clusterExport(clust, c("count_expressing_cells", "combination"), 
              envir=environment())

#args<-commandArgs(TRUE)
#path_to_imputed_files<-args[1]


#result_cols_fraction<-parApply(clust, combination, 1,
#function(row) count_expressing_cells(paste(path_to_imputed_files,  "/sample", row[1], '_', row[2] , '.tsv', sep=''),row[1], row[2])) 

result_cols_fraction<-parApply(clust, combination, 1,
function(row) count_expressing_cells(paste("~/1tb/data/SRX15446041/imputation/", strsplit(row[2], "_")[[1]][[1]], "/sample", row[1], '_', row[2] , '.tsv', sep=''),
row[1], row[2]))     
    
result_df_fraction<-result_cols_fraction[[1]]

stopCluster(clust)

for (i in 2:length(result_cols_fraction)){
  result_df_fraction<-merge(result_df_fraction, result_cols_fraction[[i]], by='row.names', all = TRUE)
  rownames(result_df_fraction)<-result_df_fraction$Row.names
  result_df_fraction<-subset(result_df_fraction, select=-c(Row.names))
 }

write.table(result_df_fraction, paste("genes_stat_fraction_imp_knn_sample002.txt", sep=""), sep="\t")