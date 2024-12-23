###merge.R
full.table<-read.table(paste("genes_stat_fraction_imp", as.character(1), ".txt", sep=""), row.names=1, header=TRUE)
print(dim(full.table))
for (i in 2:6){
  table<-read.table(paste("genes_stat_fraction_imp", as.character(i), ".txt", sep=""), row.names=1, header=TRUE)
  print(dim(table))
  full.table<-merge(full.table, table, by='row.names', all = TRUE)
  rownames(full.table)<-full.table$Row.names
  full.table<-subset(full.table, select=-c(Row.names))
} 

                                                       
colnames(table)<-lapply(colnames(table), function(x) paste(strsplit(x, "[.]")[[1]][[1]], "3", sep="."))
full.table<-merge(full.table, table, by='row.names', all = TRUE)
rownames(full.table)<-full.table$Row.names
full.table<-subset(full.table, select=-c(Row.names))

colnames(table)<-lapply(colnames(table), function(x) paste(strsplit(x, "[.]")[[1]][[1]], "9", sep="."))
full.table<-merge(full.table, table, by='row.names', all = TRUE)
rownames(full.table)<-full.table$Row.names
full.table<-subset(full.table, select=-c(Row.names))

print(colnames(full.table))
write.table(full.table, "genes_stats_fraction_imp.txt", sep="\t")
