library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(stringr)
library(readr)


draw_plot<-function(pict_data, dataset_name, simulation_name, path, method_name_list){
  pict_data %>%
   mutate(across(method, factor, levels=unlist(method_name_list))) ->pict_data
  plt<-ggplot(pict_data, aes(x=percent_expressed_true, y=percent_expressed, fill=method)) +
    geom_point(size = 1, alpha = 3/10) +# ylim(0, 0.15)+# xlim(0, 3000)+ylim(0, 0.02)+
    geom_abline(data = pict_data, aes(slope = 1, intercept = 0), colour = "red") +
    theme_bw() + scale_x_continuous(breaks=c(0, 0.5, 1)) + scale_y_continuous(breaks=c(0, 0.5, 1)) +
    theme(plot.title = element_text(size = 24), text=element_text(size=24), axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),legend.position = "none", strip.background=element_rect(colour="black",
                                    fill="white", linetype = "solid")) +
    facet_wrap(~method, ncol = 3)+
    labs(x="'True' fraction of expressing cells", y='Fraction of expressing cells after imputation', title=dataset_name)+
    guides(color = "none", fill = "none")
  svg(file=paste(path,"percent_expressed_true_", dataset_name,"_", simulation_name, '.svg', sep=''))
  print(plt)
  dev.off()
}

change_name<-function(x){
  if (x=="saver"){
    x<-"SAVER"
  } 
  if (x=="alra"){
    x<-"ALRA"
  }
  if (x=="magic"){
    x<-"MAGIC"
  }                       
  if (x=="scImpute"){
    x<-"scImpute"
  }
  if (x=="knn"){
    x<-"kNN-smoothing"
  }
  if (x=="scVI"){
    x<-"scVI"
  }
  return(x)
}

method_name_list<-list('alra', 'knn','magic','saver', 'scImpute', 'scVI')
args<-commandArgs(TRUE)
dataset_name<- args[1]
simulation_name<-args[2]
pict_data<-read.table(paste("~/1tb/paper/simulations/imputation_analysis/",simulation_name,"_", dataset_name, "_stats.txt", sep=""), row.names=1, header=TRUE) 
path<- args[3]
pict_data$method<-sapply(pict_data$method, function(x) change_name(x))
draw_plot(pict_data, dataset_name , simulation_name, path , unique(pict_data$method))