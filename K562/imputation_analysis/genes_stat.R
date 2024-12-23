library(dplyr)
library(tidyr)
library(parallel)

count_expressing_cells<-function(path, sample){        
  data<-read.table(path,row.names=1, header=TRUE)
  rownames(data)<-gsub("-", ".", rownames(data))
  results<-apply(data, 1, function(data_row) sum(data_row>0)/length(data_row)) # calculate proportion of expressing cells for each gene
  results<-as.data.frame(results)
  colnames(results)<-c(paste("sample", sample ,sep=""))
  return(results)
}
count_expression_level<-function(path, sample){ #aggregate counts and then log them
  data<-read.table(path,row.names=1, header=TRUE)
  rownames(data)<-gsub("-", ".", rownames(data))
  results<-rowSums(data)
  results<-log1p(results) 
  results<-as.data.frame(results)
  colnames(results)<-c(paste("sample", sample ,sep=""))
  return(results)
}

#
method_name<-c("raw")
sample=c("002","004", "006", "008", "01", "012", "014", "016", "018", "02", "03", "04", "05", "06", "08", "1") 

#for (i in 1:length(sample))  {
#   print(paste("~/1tb/paper/K562/preprocessed_data/sample",sample[i],'_std.txt', sep=''))
#}

result_cols_fraction<-mclapply(sample, function(x) count_expressing_cells(paste("~/1tb/paper/K562/preprocessed_data/sample",x,'_std.txt', sep=''), x), mc.cores=16)
result_df_fraction<-result_cols_fraction[[1]]

result_cols_expression<-mclapply(sample, function(x) count_expression_level(paste("~/1tb/paper/K562/preprocessed_data/sample",x,'_std.txt', sep=''), x), mc.cores=16)
result_df_expression<-result_cols_expression[[1]]

for (i in 2:length(result_cols_fraction)){
  result_df_fraction<-merge(result_df_fraction, result_cols_fraction[[i]], by='row.names', all = TRUE)
  rownames(result_df_fraction)<-result_df_fraction$Row.name
  result_df_fraction<-subset(result_df_fraction, select=-c(Row.names))
  result_df_expression<-merge(result_df_expression, result_cols_expression[[i]], by='row.names', all = TRUE)
  rownames(result_df_expression)<-result_df_expression$Row.names
  result_df_expression<-subset(result_df_expression, select=-c(Row.names))
}

write.table(result_df_fraction, paste("genes_stat_raw_fraction.txt", sep=""), sep="\t")
write.table(result_df_expression, paste("genes_stat_raw_expression.txt", sep=""), sep="\t")