library(dplyr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(stringr)
library(readr)

data_names<-list('alra'="alra", 'knn'="knn", 'magic'="magic",
                 'observed'="observed",'saver'="saver", 'scIV'="scVI", 'scimpute'="scImpute", 'observed_a'="observed")
data_labeller <- function(variable,value){
  return(data_names[value])
}

draw_plot<-function(pict_data, dataset_name, simulation_name, path, method_name_list){
  plt<-ggplot(pict_data, aes(x=percent_expressed_obs, y=delta_nrmse2, fill=method)) +                #delta_nrmse - true_mean and delta_nrmse2 - obs_mean
#percent_expressed_obs,  pseudobulk_expression_obs
    geom_point(size = 1, alpha = 1/10) + ylim(-20, 20)+#xlim(0, 1)+# xlim(0, 3000)+ylim(0, 0.02)+
    stat_summary_bin(fun.y='median', bins=10,
                     color='orange', size=1, geom='point')+
    geom_abline(data = pict_data, aes(slope = 0, intercept = 0), colour = "red")+
    theme_bw()+   
    theme(plot.title = element_text(size = 20), text=element_text(size=20), axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),legend.position = "none", strip.background=element_rect(colour="black",
                                    fill="white", linetype = "solid")) +
    facet_wrap(~method, ncol=6)+ #, labeller=data_labeller)+
    labs(x='Fraction of expressing cells before imputation', y='deltaNRMSE', title=paste(dataset_name, simulation_name, sep=", "))+
# pseudobulk, fraction of expressing cells  
    guides(color = "none", fill = "none")
  svg(file=paste(path,"delta_nrmse2_crop_", dataset_name,"_", "knn_",simulation_name, '.svg', sep=''),height=8, width=14)
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
    x<-"KNN-smoothing"
  }
  if (x=="scIV"){
    x<-"scVI"
  }
  return(x)
}

#method_name_list<-list('alra', 'knn','magic','saver', 'scImpute', 'scVI')
args<-commandArgs(TRUE)
dataset_name<- args[1]
simulation_name<-args[2]
pict_data<-read.table(paste("~/1tb/paper/simulations/imputation_analysis/",simulation_name,"_", dataset_name, "_knn_stats.txt", sep=""), row.names=1, header=TRUE) 
path<-args[3]
pict_data$method<-sapply(pict_data$method, function(x) paste("k=", as.character(x), sep="")) # change_name(x)
pict_data$method<-factor(pict_data$method, levels=c("k=2", "k=3", "k=4", "k=5","k=6", "k=7", "k=8", "k=9", "k=10", "k=11", "k=12","k=13"))
draw_plot(pict_data, dataset_name , simulation_name, path , unique(pict_data$method))
