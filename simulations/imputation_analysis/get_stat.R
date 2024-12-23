library(dplyr)
library(tidyverse)
library(Seurat)
library(stringr)
library(readr)

norm_format<-function(path){
  mtrx<-read.table(path, row.names=1, header = TRUE) # table with gene names in the first column and cell names in the
  mtrx<-as.matrix(mtrx)
  return(mtrx)
}

rmse<-function(original_data, imputed_data, method_name){
  common_rownames<-intersect(rownames(original_data), rownames(imputed_data))
  common_colnames<-intersect(colnames(original_data), colnames(imputed_data))
  original_data<-original_data[common_rownames,common_colnames]
  imputed_data<-imputed_data[common_rownames,common_colnames]
  calc_error<-(original_data-imputed_data) #error
  calc_error<-apply(calc_error,c(1,2), function(x) x^2) #squared_error
  calc_error<-sqrt(rowSums(calc_error)*(1/dim(calc_error)[2]))
  print("ok-rmse")
  names(calc_error)<-common_rownames
  calc_error<-as.data.frame(calc_error)
  colnames(calc_error)<-c(method_name)
  print("exit")
  return(calc_error)
}

calc_deviation_help<-function(data_row){
  sd_var<-sd(data_row)
  return(sd_var)
}

calc_deviation<-function(data, method_name){
  rowname<-rownames(data)
  data<-apply(data,1, function(x) calc_deviation_help(x))
  data<-as.data.frame(data)
  colnames(data)<-method_name
  rownames(data)<-rowname
  return(data)
}

calc_mean_help<-function(data_row){
  mean_var<-mean(data_row)
  return(mean_var)
}

calc_mean<-function(data, method_name){
  rowname<-rownames(data)
  data<-apply(data,1, function(x) calc_mean_help(x))
  data<-as.data.frame(data)
  colnames(data)<-method_name
  rownames(data)<-rowname
  return(data)
}

cell_percent_help<-function(row){
  percent<-sum(row>0)/length(row)
  return(percent)
}

cell_percent<-function(data, name){
  rowname<-rownames(data)
  data<-apply(data,1, function(x) cell_percent_help(x))
  print(typeof(data))
  data<-as.data.frame(data)
  colnames(data)<-name
  rownames(data)<-rowname
  return(data)
}

expression_level<-function(observed_data){
  pseudobulk<-rowSums(observed_data)
  pseudobulk<-as.data.frame(pseudobulk)
  rownames(pseudobulk)<-rownames(observed_data)
  colnames(pseudobulk)<-'pseudobulk_expression_obs'
  return(pseudobulk)
}

calc_nrmse<-function(data){
  data$nrmse<-data$rmse/data$true_mean
  return(data)
}

do_analysis<-function(true_data,observed_data, imputed_path_list, method_name_list, name){
  true_data<-read.table(true_data, header=TRUE)
  Seurat_object<-CreateSeuratObject(counts=true_data, min.cells=3, min.features=200)
  Seurat_object<-NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  true_data<-Seurat_object[["RNA"]]$data
  rownames(true_data)<-unlist(lapply(rownames(true_data), function(x) paste('X', x, sep='')))
  observed_data<-read.table(observed_data, header=TRUE)
  Seurat_object<-CreateSeuratObject(counts=observed_data, min.cells=3, min.features=200)
  Seurat_object<-NormalizeData(Seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
  observed_data<-Seurat_object[["RNA"]]$data
  rownames(observed_data)<-unlist(lapply(rownames(observed_data), function(x) paste('X', x, sep='')))
  imputed_data_list<-list()
  for (i in 1:length(imputed_path_list)){
    imputed_data_list[[i]]<-norm_format(imputed_path_list[[i]])
    print(method_name_list[[i]])
    if (substr(rownames(imputed_data_list[[i]])[1], 1, 1) != 'X'){
      rownames(imputed_data_list[[i]])<-unlist(lapply(rownames(imputed_data_list[[i]]), function(x) paste('X', x, sep='')))
      print(rownames(imputed_data_list[[i]])[1:10])
    }
    print(as.character(i))
  }
  imputed_data_list[[i+1]]<-observed_data
  print('Loop done, all data read')
  for (i in 1:length(imputed_data_list)){
    calc_error<-rmse(true_data, imputed_data_list[[i]], method_name_list[[i]])
    if (i==1){
      pict_data<-calc_error
    } else{
      pict_data <- merge(pict_data, calc_error, by=0, all=TRUE)
      rownames(pict_data)<-pict_data[, 'Row.names']
      pict_data<-subset(pict_data, select=-c(Row.names)) 
    }
  }
  pict_data$genes<-rownames(pict_data)
  pict_data<-pivot_longer(pict_data, cols=unlist(method_name_list), names_to = "method", values_to='rmse')
  print(head(pict_data))
  for (i in 1:length(imputed_data_list)){
    calc_perc<-cell_percent(imputed_data_list[[i]], method_name_list[[i]])
    
    if (i==1){
      pict_data2<-calc_perc
    } else{
      pict_data2 <- merge(pict_data2, calc_perc, by=0, all=TRUE)
      rownames(pict_data2)<-pict_data2[, 'Row.names']
      pict_data2<-subset(pict_data2, select=-c(Row.names)) 
    }
  }
  print('Cell percent calculated')
  pict_data2$genes<-rownames(pict_data2)
  pict_data2<-pivot_longer(pict_data2, cols=unlist(method_name_list), names_to = "method", values_to='percent_expressed')
  pict_data<-merge(pict_data, pict_data2, by=c("genes", "method"))
  true_mean<-calc_mean(true_data, 'true_mean')
  true_mean$genes<-rownames(true_mean)
  true_perc<-cell_percent(true_data, "percent_expressed_true")
  true_perc$genes<-rownames(true_perc)
  new_data<-merge(pict_data, true_mean, by.x=c('genes'),by.y=c('genes'))
  new_data<-merge(new_data, true_perc, by.x=c('genes'),by.y=c('genes'))
  obs_mean<-calc_mean(observed_data, "obs_mean")
  obs_mean$genes<-rownames(obs_mean)
  new_data<-merge(new_data, obs_mean, by.x=c('genes'),by.y=c('genes'))
  #nrmse
  new_data<-calc_nrmse(new_data)
  print("nrmse calculated")
  help_data<-new_data %>% filter(method=='observed')
  help_data$percent_expressed_obs<-help_data$percent_expressed
  help_data$nrmse_obs<-help_data$nrmse
  help_data<-subset(help_data, select=c(genes, percent_expressed_obs, nrmse_obs))
  newly_data<-new_data %>% filter(method!='observed')
  new_pict_data<-merge(newly_data, help_data, by.x=c('genes'), by.y=c('genes'))
  new_pict_data$delta_nrmse<-new_pict_data$nrmse_obs-new_pict_data$nrmse
  new_data$nrmse2<-new_data$rmse/new_data$obs_mean # nrmse2 is normalized by observed mean expression
  help_data<-new_data %>% filter(method=='observed')
  help_data$nrmse2_obs<-help_data$nrmse2
  help_data<-subset(help_data, select=c(genes, nrmse2_obs))
  new_data<-new_data %>% filter(method!='observed')
  new_data<-subset(new_data, select=c(genes, nrmse2, method))
  new_pict_data2<-merge(new_data, help_data, by.x=c('genes'), by.y=c('genes'))
  new_pict_data2$delta_nrmse2<-new_pict_data2$nrmse2_obs-new_pict_data2$nrmse2
  full_pict_data<-merge(new_pict_data, new_pict_data2, by=c('genes', 'method'))
  write.table(full_pict_data, paste(name,'_knn_stats.txt', sep=""))
}

args<-commandArgs(TRUE)
dataset_name<-args[1]
sim_type<-"zinb-wave"
name<-paste(sim_type, dataset_name, sep="_")

observed_data<-paste("~/1tb/paper/simulations/ZINB-WaVE/", name, "_obs_counts.tsv", sep="")
true_data<-paste("~/1tb/paper/simulations/ZINB-WaVE/", name, "_true_counts.tsv", sep="")

method_name_list<-c(3, 4, 5 ,6, 7, 8, 9, 10, 11, 12, 13)

imputed_path_list<-c(paste("~/1tb/paper/simulations/imputed_data/knn_params/", name, "_knn_", as.character(2), ".tsv", sep=""))

for (k in method_name_list){
   imputed_path_list<-c(imputed_path_list, (paste("~/1tb/paper/simulations/imputed_data/knn_params/", name, "_knn_", as.character(k), ".tsv", sep="")))
} 
print(imputed_path_list)
method_name_list<-c('2', '3', '4', '5' ,'6', '7', '8', '9', '10', '11', '12', '13', 'observed')

#imputed_path_list<-list(paste("~/1tb/paper/simulations/imputed_data/knn_params/", as.chracter(name), "_knn_",  "/.tsv", sep=""),
#                        "~/1tb/paper/simulations/imputation/knn-smoothing/symsim_ciliated_knn_8.tsv",
#			"~/1tb/paper/simulations/imputation/magic/symsim_ciliated_magic_3.tsv",
#                        "~/1tb/paper/simulations/imputation/saver/symsim_ciliated_saver.tsv",
#                        "~/1tb/paper/simulations/imputation/scimpute/symsim_ciliated_scimpute_0.5.tsv",
#                        "~/1tb/paper/simulations/imputation/scvi/symsim_ciliated_scvi.tsv")

#method_name_list<-list('alra', 'knn', 'magic', 'saver', 'scImpute','scVI', 'observed')

do_analysis(true_data,observed_data, imputed_path_list, method_name_list, name)