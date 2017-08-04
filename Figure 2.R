###Generating scale plots for Figure 1a,b with relative importance shaded in the background####
library(lubridate)
library(mgcv)
library(ggplot2)



###calculating importance averages for groups of models###

#load data from location of local machine. Available from https://knb.ecoinformatics.org/#view/doi:10.5063/F17P8WGK
fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/EEE_Cs_melanura_7_27_17.csv", header=T)

#load date and month variable
fulldata$month_num<-month(fulldata$month_year)



#log of mean abundance per night
fulldata$ln_Mean_Abundance_perTrapnight<-log(fulldata$Mean_Abundance_perTrapnight+1)


#Figure 2a
mel_abund_forestsemiclim_graph<-data.frame(buffer=as.numeric(),lag=as.numeric(),aic_score=as.numeric())

buffers<-c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

##testing out different combinations of spatial scales and temporal lags for forested semi-permanent wetland and phdi analysis. Best scales/lags for the other covariates in the model were determined in a separate analysis. Code for this is available from the authors.
for (i in 0:12){
  for (j in 1:length(buffers)){
    
    tryCatch({   
      
      model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_50, k=3)+s(prop_PFO_Ever_3000, k=3)+s(prop_PFO_Decid_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+s(eval(as.name(paste("prop_PFO_F_", buffers[j], sep=""))), eval(as.name(paste("prev_phdi_",i,sep=""))), k=4)+ti(prop_PFO_3000,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5&fulldata$prop_PFO_2000<.15,], select=T)
      
    }, error=function(e){})   
    
    mel_abund_forestsemiclim_graph<-rbind(mel_abund_forestsemiclim_graph, data.frame(buffer=buffers[j],lag=i,aic_score=model$aic))
    
  }}





#Calcuate a delta AIC column
mel_abund_forestsemiclim_graph$delt_AIC<- mel_abund_forestsemiclim_graph$aic_score - min(mel_abund_forestsemiclim_graph$aic_score)


#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_forestsemiclim_graph$weight_AIC<-exp( -0.5 * mel_abund_forestsemiclim_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_forestsemiclim_graph$delt_AIC))



#calcuating importance data for mel forestsemiclim

importance_data_forestforestsemiclim_abund<-data.frame()


for (k in 1:nrow(mel_abund_forestsemiclim_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_forestsemiclim_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_forestsemiclim_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric(), w=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_forestsemiclim_graph[k, "buffer"],x=importance_val, y=species, w=mel_abund_forestsemiclim_graph[k, "lag"])
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"


  colnames(little_data)[2]<-"relative_importance"
  

  
  colnames(little_data)[4]<-"lags"
  
  #binding the data together
  importance_data_forestforestsemiclim_abund<-rbind(importance_data_forestforestsemiclim_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)

 
importance_abund_forestforestsemiclim_bylag<-aggregate(importance_data_forestforestsemiclim_abund$relative_importance~importance_data_forestforestsemiclim_abund$lags, FUN=mean)



###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_forestforestsemiclim_bylag<-data.frame(approx(x=importance_abund_forestforestsemiclim_bylag$'importance_data_forestforestsemiclim_abund$lags', y=log(importance_abund_forestforestsemiclim_bylag$'importance_data_forestforestsemiclim_abund$relative_importance'), n=10000))



#setting relative sizes for plotted points
sizes<-c( 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

#generating the plot
ggplot()+
  geom_rect(data=interp_log_importance_abund_forestforestsemiclim_bylag, aes(xmin=interp_log_importance_abund_forestforestsemiclim_bylag$x-0.0012, xmax=interp_log_importance_abund_forestforestsemiclim_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_abund_forestforestsemiclim_bylag$y))+

  geom_point(data=mel_abund_forestsemiclim_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1
            ,breaks=c(-80, -50, -20), labels=c("1.80e-35 Low","1.93e-22","2.06e-9  High")
   )+
  scale_size_manual("Buffer Size (m)\nfor Forested\nSemi-perm. Wetland", values=sizes/600) +
  geom_line(data=mel_abund_forestsemiclim_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +

  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))






#Figure 2b
mel_abund_forest_graph<-data.frame(buffer=as.numeric(),lag=as.numeric(),aic_score=as.numeric())

buffers<-c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (i in 0:12){
  for (j in 1:length(buffers)){
    
    tryCatch({   
      
      model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_50, k=3)+s(prop_PFO_Ever_3000, k=3)+s(prop_PFO_Decid_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+ti(eval(as.name(paste("prop_PFO_", buffers[j], sep=""))), eval(as.name(paste("prev_phdi_",i,sep=""))), k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
      
    }, error=function(e){})   
    
    mel_abund_forest_graph<-rbind(mel_abund_forest_graph, data.frame(buffer=buffers[j],lag=i,aic_score=model$aic))
    
  }}




mel_abund_forest_graph$delt_AIC<- mel_abund_forest_graph$aic_score - min(mel_abund_forest_graph$aic_score)


mel_abund_forest_graph$weight_AIC<-exp( -0.5 * mel_abund_forest_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_forest_graph$delt_AIC))


importance_data_forestforest_abund<-data.frame()


for (k in 1:nrow(mel_abund_forest_graph)){
  

  importance_val<-sum(mel_abund_forest_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_forest_graph["weight_AIC"])
  
  

  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric(), w=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_forest_graph[k, "buffer"],x=importance_val, y=species, w=mel_abund_forest_graph[k, "lag"])
  

  colnames(little_data)[1]<-"buffers"

  colnames(little_data)[2]<-"relative_importance"

  colnames(little_data)[4]<-"lags"
  

  importance_data_forestforest_abund<-rbind(importance_data_forestforest_abund, little_data)
} 




importance_abund_forestforest_bylag<-aggregate(importance_data_forestforest_abund$relative_importance~importance_data_forestforest_abund$lags, FUN=mean)


interp_log_importance_abund_forestforest_bylag<-data.frame(approx(x=importance_abund_forestforest_bylag$'importance_data_forestforest_abund$lags', y=log(importance_abund_forestforest_bylag$'importance_data_forestforest_abund$relative_importance'), n=10000))


sizes<-c( 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)



ggplot()+
  geom_rect(data=interp_log_importance_abund_forestforest_bylag, aes(xmin=interp_log_importance_abund_forestforest_bylag$x-0.0012, xmax=interp_log_importance_abund_forestforest_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_abund_forestforest_bylag$y))+

  geom_point(data=mel_abund_forest_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1
                      ,breaks=c(-35, -20, -7), labels=c("6.30e-16 Low","2.06e-9","0.0009  High")
  )+
  scale_size_manual("Buffer Size (m)\nfor Forested Wetland", values=sizes/600) +
  geom_line(data=mel_abund_forest_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +

  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))








#Figure 2c
mel_infect_forestsemiclim_graph<-data.frame(buffer=as.numeric(),lag=as.numeric(),aic_score=as.numeric())

buffers<-c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (i in 0:12){
  for (j in 1:length(buffers)){
    

      
      model<-gam(EEE_MIR_PRES~Prop_Gravid_Trap+s(ln_Mean_Abundance_perTrapnight, k=4)+s(prop_PEM_2000, k=4)+s(prop_PSS_1500,k=4)+s(prop_PFO_Ever_5000, k=4)+s(prop_PFO_Decid_5000, k=4)+s(ave_size_PFO_5000_zeros, k=4)+s(eval(as.name(paste("PFO_F_Relative_Temp_", buffers[j], sep=""))), eval(as.name(paste("prev_phdi_",i,sep=""))), k=4)+s(impvMN_5000, k=4)+s(month_num,k=5),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11,],family=binomial, select=T)
      

    mel_infect_forestsemiclim_graph<-rbind(mel_infect_forestsemiclim_graph, data.frame(buffer=buffers[j],lag=i,aic_score=model$aic))
    
  }}






mel_infect_forestsemiclim_graph$delt_AIC<- mel_infect_forestsemiclim_graph$aic_score - min(mel_infect_forestsemiclim_graph$aic_score)



mel_infect_forestsemiclim_graph$weight_AIC<-exp( -0.5 * mel_infect_forestsemiclim_graph$delt_AIC)/sum(exp( -0.5 * mel_infect_forestsemiclim_graph$delt_AIC))



importance_data_forestforestsemiclim_abund<-data.frame()


for (k in 1:nrow(mel_infect_forestsemiclim_graph)){
  

  importance_val<-sum(mel_infect_forestsemiclim_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_infect_forestsemiclim_graph["weight_AIC"])
  
  

  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric(), w=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_infect_forestsemiclim_graph[k, "buffer"],x=importance_val, y=species, w=mel_infect_forestsemiclim_graph[k, "lag"])
  

  colnames(little_data)[1]<-"buffers"

  colnames(little_data)[2]<-"relative_importance"
  

  colnames(little_data)[4]<-"lags"
  

  importance_data_forestforestsemiclim_abund<-rbind(importance_data_forestforestsemiclim_abund, little_data)
} 






importance_infect_forestforestsemiclim_bylag<-aggregate(importance_data_forestforestsemiclim_abund$relative_importance~importance_data_forestforestsemiclim_abund$lags, FUN=mean)


interp_log_importance_infect_forestforestsemiclim_bylag<-data.frame(approx(x=importance_infect_forestforestsemiclim_bylag$'importance_data_forestforestsemiclim_abund$lags', y=log(importance_infect_forestforestsemiclim_bylag$'importance_data_forestforestsemiclim_abund$relative_importance'), n=10000))

sizes<-c( 50,100,200,500, 750, 1000, 1500, 2000, 3000, 4000, 5000)



ggplot()+
  geom_rect(data=interp_log_importance_infect_forestforestsemiclim_bylag, aes(xmin=interp_log_importance_infect_forestforestsemiclim_bylag$x-0.0012, xmax=interp_log_importance_infect_forestforestsemiclim_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_infect_forestforestsemiclim_bylag$y))+
  

  geom_point(data=mel_infect_forestsemiclim_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1,breaks=c(-27, -15, -5), labels=c("1.88e-12 Low","3.06e-7","0.007 High")
  )+
  scale_size_manual("Buffer Size (m) for\nRelative Semi-perm.\nComposition", values=sizes/600) +
  geom_line(data=mel_infect_forestsemiclim_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +

  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))









