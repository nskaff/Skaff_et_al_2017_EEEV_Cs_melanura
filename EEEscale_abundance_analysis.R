#Figures for Chapter 2 manuscript

###creating scale plots with importance shaded in the background####


library(zoo)
library(mgcv)
library(scales)
library(lubridate)
library(ggplot2)

###calculating importance averages for groups of models###


fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/EEE_data_mle_3_10_17.csv", header=T)

fulldata$month_year<-as.Date(fulldata$month_year)
fulldata$month_num<-month(fulldata$month_year)

# fulldata$EEE_MIR_PRES<-fulldata$IR
# fulldata[fulldata$EEE_MIR_PRES>0&!is.na(fulldata$EEE_MIR_PRES),"EEE_MIR_PRES"]<-1
# fulldata[fulldata$EEE_MIR_PRES==0&!is.na(fulldata$EEE_MIR_PRES),"EEE_MIR_PRES"]<-0
# 
# fulldata$EEE_MIR_PRES<-as.factor(fulldata$EEE_MIR_PRES)

buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


#overlaying all WNV vector abundance plots for connectivity


#melanura model 1####

fulldata$ln_Mean_Abundance_perTrapnight<-log(fulldata$Mean_Abundance_perTrapnight+1)




# 
# ####melanura forested semi wetland and clim for abundance####
# for (i in which(colnames(fulldata)=="PFO_PSS_temp_50"):which(colnames(fulldata)=="PEM_temp_750")){
#   
#   fulldata[,paste("prop_", colnames(fulldata)[i], sep="")]<-fulldata[,i]/fulldata[,paste("FIRST_POLY_AREA_", gsub("\\D", "",as.name(as.character(colnames(fulldata)[i]))), sep="")]
#   
# }

##melanura for forested semiclim
mel_abund_forestsemiclim_graph<-data.frame(buffer=as.numeric(),lag=as.numeric(),aic_score=as.numeric())

buffers<-c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

##for forestsemiclim
for (i in 0:12){
  for (j in 1:length(buffers)){
    
    tryCatch({   
      
      model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_50, k=3)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+s(eval(as.name(paste("prop_PFO_F_", buffers[j], sep=""))), eval(as.name(paste("prev_phdi_",i,sep=""))), k=4)+ti(prop_PFO_3000,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5&fulldata$prop_PFO_2000<.15,], select=T)
      
    }, error=function(e){})   
    
    mel_abund_forestsemiclim_graph<-rbind(mel_abund_forestsemiclim_graph, data.frame(buffer=buffers[j],lag=i,aic_score=model$aic))
    
  }}





#melanura forestsemiclim calcuate a delta AIC column
mel_abund_forestsemiclim_graph$delt_AIC<- mel_abund_forestsemiclim_graph$aic_score - min(mel_abund_forestsemiclim_graph$aic_score)


#melanura forestsemiclim aic weights
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
  # colnames(little_data)[1]<-paste(colnames(mel_abund_forestsemiclim_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_forestsemiclim_graph[k,]), "propmodels",values[k], sep="_")
  
  colnames(little_data)[4]<-"lags"
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_forestforestsemiclim_abund<-rbind(importance_data_forestforestsemiclim_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_emergsemiclim_byspecies<-aggregate(importance_data_emergsemiclim_binaryWNV$relative_importance~importance_data_emergsemiclim_binaryWNV$y+importance_data_emergsemiclim_binaryWNV$lags, FUN=max)
# 
importance_abund_forestforestsemiclim_bylag<-aggregate(importance_data_forestforestsemiclim_abund$relative_importance~importance_data_forestforestsemiclim_abund$lags, FUN=mean)

# importance_forestsemiclim_bybuffer<-aggregate(importance_data_forestsemiclim_binaryWNV$relative_importance~importance_data_forestsemiclim_binaryWNV$buffers, FUN=mean)

#importance_window_forestsemiclim<-data.frame(rollapply(importance_forestsemiclim, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_forestforestsemiclim_bylag<-data.frame(approx(x=importance_abund_forestforestsemiclim_bylag$'importance_data_forestforestsemiclim_abund$lags', y=log(importance_abund_forestforestsemiclim_bylag$'importance_data_forestforestsemiclim_abund$relative_importance'), n=10000))
# 
# interp_log_importance_forestsemiclim_bybuffer<-data.frame(approx(x=importance_forestsemiclim_bybuffer$'importance_data_forestsemiclim_binaryWNV$buffers', y=log(importance_forestsemiclim_bybuffer$'importance_data_forestsemiclim_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c( 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/forestsemi_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_forestforestsemiclim_bylag, aes(xmin=interp_log_importance_abund_forestforestsemiclim_bylag$x-0.0012, xmax=interp_log_importance_abund_forestforestsemiclim_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_abund_forestforestsemiclim_bylag$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_point(data=mel_abund_forestsemiclim_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1
            ,breaks=c(-80, -50, -20), labels=c("1.80e-35 Low","1.93e-22","2.06e-9  High")
   )+
  scale_size_manual("Buffer Size (m)\nfor Forested\nSemi-perm. Wetland", values=sizes/600) +
  geom_line(data=mel_abund_forestsemiclim_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()





####

####melanura forested wetland and clim for abund####

##melanura for frstsemiclim
mel_abund_forest_graph<-data.frame(buffer=as.numeric(),lag=as.numeric(),aic_score=as.numeric())

buffers<-c(500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

##for forest
for (i in 0:12){
  for (j in 1:length(buffers)){
    
    tryCatch({   
      
      model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_50, k=3)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+ti(eval(as.name(paste("prop_PFO_", buffers[j], sep=""))), eval(as.name(paste("prev_phdi_",i,sep=""))), k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
      
    }, error=function(e){})   
    
    mel_abund_forest_graph<-rbind(mel_abund_forest_graph, data.frame(buffer=buffers[j],lag=i,aic_score=model$aic))
    
  }}





#melanura forest calcuate a delta AIC column
mel_abund_forest_graph$delt_AIC<- mel_abund_forest_graph$aic_score - min(mel_abund_forest_graph$aic_score)


#melanura forest aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_forest_graph$weight_AIC<-exp( -0.5 * mel_abund_forest_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_forest_graph$delt_AIC))



#calcuating importance data for mel forest

importance_data_forestforest_abund<-data.frame()


for (k in 1:nrow(mel_abund_forest_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_forest_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_forest_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric(), w=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_forest_graph[k, "buffer"],x=importance_val, y=species, w=mel_abund_forest_graph[k, "lag"])
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_forest_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_forest_graph[k,]), "propmodels",values[k], sep="_")
  
  colnames(little_data)[4]<-"lags"
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_forestforest_abund<-rbind(importance_data_forestforest_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_emergsemiclim_byspecies<-aggregate(importance_data_emergsemiclim_binaryWNV$relative_importance~importance_data_emergsemiclim_binaryWNV$y+importance_data_emergsemiclim_binaryWNV$lags, FUN=max)
# 
importance_abund_forestforest_bylag<-aggregate(importance_data_forestforest_abund$relative_importance~importance_data_forestforest_abund$lags, FUN=mean)

# importance_forest_bybuffer<-aggregate(importance_data_forest_binaryWNV$relative_importance~importance_data_forest_binaryWNV$buffers, FUN=mean)

#importance_window_forest<-data.frame(rollapply(importance_forest, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_forestforest_bylag<-data.frame(approx(x=importance_abund_forestforest_bylag$'importance_data_forestforest_abund$lags', y=log(importance_abund_forestforest_bylag$'importance_data_forestforest_abund$relative_importance'), n=10000))
# 
# interp_log_importance_forest_bybuffer<-data.frame(approx(x=importance_forest_bybuffer$'importance_data_forest_binaryWNV$buffers', y=log(importance_forest_bybuffer$'importance_data_forest_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c( 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/forest_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_forestforest_bylag, aes(xmin=interp_log_importance_abund_forestforest_bylag$x-0.0012, xmax=interp_log_importance_abund_forestforest_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_abund_forestforest_bylag$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_point(data=mel_abund_forest_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1
                      ,breaks=c(-35, -20, -7), labels=c("6.30e-16 Low","2.06e-9","0.0009  High")
  )+
  scale_size_manual("Buffer Size (m)\nfor Forested Wetland", values=sizes/600) +
  geom_line(data=mel_abund_forest_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()





################PHDI X wetland panel plot#############################
sizes<-c( 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)



semiforest_phdi<-ggplot()+
  geom_rect(data=interp_log_importance_abund_forestforestsemiclim_bylag, aes(xmin=interp_log_importance_abund_forestforestsemiclim_bylag$x-0.0012, xmax=interp_log_importance_abund_forestforestsemiclim_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_abund_forestforestsemiclim_bylag$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_point(data=mel_abund_forestsemiclim_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1
                       ,breaks=c(-80, -50, -20), labels=c("1.80e-35 Low","1.93e-22","2.06e-9  High")
  )+
  scale_size_manual("Buffer Size (m)\nfor Forested\nSemi-perm. Wetland", values=sizes/600) +
  geom_line(data=mel_abund_forestsemiclim_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))





forest_phdi<-ggplot()+
  geom_rect(data=interp_log_importance_abund_forestforest_bylag, aes(xmin=interp_log_importance_abund_forestforest_bylag$x-0.0012, xmax=interp_log_importance_abund_forestforest_bylag$x+0.0012, ymin=-4,ymax=2,fill=interp_log_importance_abund_forestforest_bylag$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_point(data=mel_abund_forest_graph,aes(lag,scale(aic_score),  size=as.factor(buffer)), shape=19, stroke=1.5)+
  
  scale_x_continuous(name="# Months Lag in PHDI", breaks=seq(0,12,1)) + 
  scale_y_continuous(name="Centered AIC Score", breaks=c(-4,-2,0,2))+ 
  scale_fill_distiller("Importance\nof Lag",palette = "YlOrRd", direction=1
                       ,breaks=c(-35, -20, -7), labels=c("6.30e-16 Low","2.06e-9","0.0009  High")
  )+
  scale_size_manual("Buffer Size (m)\nfor Forested Wetland", values=sizes/600) +
  geom_line(data=mel_abund_forest_graph,aes(lag,scale(aic_score), group=as.factor(buffer)), color="grey") +
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/wetland_PHDI_SCALE_Panel_melabundance.png", width=3400, height=4800, res=300)

multiplot(semiforest_phdi, forest_phdi, cols=1)


dev.off()





####

####melanura emergent wetland NO CLIM for abund####

##melanura for emerg
buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
mel_abund_emerg_graph<-data.frame(buffer=as.numeric(),aic_score=as.numeric())
##for emerg
for (i in 1:length(buffers)){
  
  model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(eval(as.name(paste("prop_PEM_",buffers[i], ".y",sep=""))), k=4)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+ti(prop_PFO_3000, prev_phdi_0, k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
  
  
  mel_abund_emerg_graph<-rbind(mel_abund_emerg_graph, data.frame(buffer=buffers[i],aic_score=model$aic))
  
}





#melanura emerg calcuate a delta AIC column
mel_abund_emerg_graph$delt_AIC<- mel_abund_emerg_graph$aic_score - min(mel_abund_emerg_graph$aic_score)


#melanura emerg aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_emerg_graph$weight_AIC<-exp( -0.5 * mel_abund_emerg_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_emerg_graph$delt_AIC))



#calcuating importance data for mel emerg

importance_data_emergemerg_abund<-data.frame()


for (k in 1:nrow(mel_abund_emerg_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_emerg_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_emerg_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_emerg_graph[k, "buffer"],x=importance_val, y=species)
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_emerg_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_emerg_graph[k,]), "propmodels",values[k], sep="_")
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_emergemerg_abund<-rbind(importance_data_emergemerg_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_emerg_byspecies<-aggregate(importance_data_emerg_binaryWNV$relative_importance~importance_data_emerg_binaryWNV$y+importance_data_emerg_binaryWNV$lags, FUN=max)
# 
# importance_abund_emergemerg_bylag<-aggregate(importance_data_emergemerg_abund$relative_importance~importance_data_emergemerg_abund$buffers, FUN=mean)

# importance_emerg_bybuffer<-aggregate(importance_data_emerg_binaryWNV$relative_importance~importance_data_emerg_binaryWNV$buffers, FUN=mean)

#importance_window_emerg<-data.frame(rollapply(importance_emerg, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_emergemerg_bybuffer<-data.frame(approx(x=importance_data_emergemerg_abund$buffers, y=log(importance_data_emergemerg_abund$relative_importance), n=10000))
# 
# interp_log_importance_emerg_bybuffer<-data.frame(approx(x=importance_emerg_bybuffer$'importance_data_emerg_binaryWNV$buffers', y=log(importance_emerg_bybuffer$'importance_data_emerg_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/emerg_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_emergemerg_bybuffer, aes(xmin=interp_log_importance_abund_emergemerg_bybuffer$x-0.45, xmax=interp_log_importance_abund_emergemerg_bybuffer$x+0.45, ymin=-3,ymax=1.5,fill=interp_log_importance_abund_emergemerg_bybuffer$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_line(data=mel_abund_emerg_graph,aes(buffer,scale(aic_score)), color="grey") +
  geom_point(data=mel_abund_emerg_graph,aes(buffer,scale(aic_score)), shape=19, stroke=1.5, size=3)+
  
  scale_x_continuous(name="Buffer Radius (m) for Emergent Wetland") + 
  scale_y_continuous(name="Centered AIC Score")+ 
  scale_fill_distiller("Importance\nof Scale",palette = "YlOrRd", direction=1
                      ,breaks=c(-175, -90,  -10), labels=c("9.96e-77 Low", 8.19e-40,"4.54e-5 High")
  )+
  # scale_size_manual("Buffer Size (m) for\nForested Wetland Metric", values=sizes/600) +
  
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  #guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()



 


####melanura shrub wetland NO CLIM for abund####

##melanura for shrub
buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
mel_abund_shrub_graph<-data.frame(buffer=as.numeric(),aic_score=as.numeric())
##for shrub
for (i in 1:length(buffers)){
  
  model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(eval(as.name(paste("prop_PSS_",buffers[i],sep=""))), k=4)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PEM_50, k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+ti(prop_PFO_3000, prev_phdi_0, k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
  
  
  mel_abund_shrub_graph<-rbind(mel_abund_shrub_graph, data.frame(buffer=buffers[i],aic_score=model$aic))
  
}





#melanura shrub calcuate a delta AIC column
mel_abund_shrub_graph$delt_AIC<- mel_abund_shrub_graph$aic_score - min(mel_abund_shrub_graph$aic_score)


#melanura shrub aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_shrub_graph$weight_AIC<-exp( -0.5 * mel_abund_shrub_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_shrub_graph$delt_AIC))



#calcuating importance data for mel shrub

importance_data_shrubshrub_abund<-data.frame()


for (k in 1:nrow(mel_abund_shrub_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_shrub_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_shrub_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_shrub_graph[k, "buffer"],x=importance_val, y=species)
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_shrub_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_shrub_graph[k,]), "propmodels",values[k], sep="_")
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_shrubshrub_abund<-rbind(importance_data_shrubshrub_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_shrub_byspecies<-aggregate(importance_data_shrub_binaryWNV$relative_importance~importance_data_shrub_binaryWNV$y+importance_data_shrub_binaryWNV$lags, FUN=max)
# 
# importance_abund_shrubshrub_bylag<-aggregate(importance_data_shrubshrub_abund$relative_importance~importance_data_shrubshrub_abund$buffers, FUN=mean)

# importance_shrub_bybuffer<-aggregate(importance_data_shrub_binaryWNV$relative_importance~importance_data_shrub_binaryWNV$buffers, FUN=mean)

#importance_window_shrub<-data.frame(rollapply(importance_shrub, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_shrubshrub_bybuffer<-data.frame(approx(x=importance_data_shrubshrub_abund$buffers, y=log(importance_data_shrubshrub_abund$relative_importance), n=10000))
# 
# interp_log_importance_shrub_bybuffer<-data.frame(approx(x=importance_shrub_bybuffer$'importance_data_shrub_binaryWNV$buffers', y=log(importance_shrub_bybuffer$'importance_data_shrub_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/shrub_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_shrubshrub_bybuffer, aes(xmin=interp_log_importance_abund_shrubshrub_bybuffer$x-0.45, xmax=interp_log_importance_abund_shrubshrub_bybuffer$x+0.45, ymin=-2,ymax=2.5,fill=interp_log_importance_abund_shrubshrub_bybuffer$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_line(data=mel_abund_shrub_graph,aes(buffer,scale(aic_score)), color="grey") +
  geom_point(data=mel_abund_shrub_graph,aes(buffer,scale(aic_score)), shape=19, stroke=1.5, size=3)+
  
  scale_x_continuous(name="Buffer Radius (m) for Scrub/Shrub Wetland") + 
  scale_y_continuous(name="Centered AIC Score")+ 
  scale_fill_distiller("Importance\nof Scale",palette = "YlOrRd", direction=1
                       ,breaks=c(-55, -32.5,  -10), labels=c("1.30e-24 Low", 7.68e-15,"4.54e-5 High")
  )+
  # scale_size_manual("Buffer Size (m) for\nForested Wetland Metric", values=sizes/600) +
  
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  #guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()



####melanura evergreen forested wetland NO CLIM for abund####

##melanura for everg
buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
mel_abund_everg_graph<-data.frame(buffer=as.numeric(),aic_score=as.numeric())
##for everg
for (i in 1:length(buffers)){
  
  model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(eval(as.name(paste("PFO_Ever.x_",buffers[i], sep=""))), k=4)+s(prop_PEM_50, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+ti(prop_PFO_3000, prev_phdi_0, k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
  
  
  mel_abund_everg_graph<-rbind(mel_abund_everg_graph, data.frame(buffer=buffers[i],aic_score=model$aic))
  
}





#melanura everg calcuate a delta AIC column
mel_abund_everg_graph$delt_AIC<- mel_abund_everg_graph$aic_score - min(mel_abund_everg_graph$aic_score)


#melanura everg aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_everg_graph$weight_AIC<-exp( -0.5 * mel_abund_everg_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_everg_graph$delt_AIC))



#calcuating importance data for mel everg

importance_data_evergeverg_abund<-data.frame()


for (k in 1:nrow(mel_abund_everg_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_everg_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_everg_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_everg_graph[k, "buffer"],x=importance_val, y=species)
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_everg_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_everg_graph[k,]), "propmodels",values[k], sep="_")
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_evergeverg_abund<-rbind(importance_data_evergeverg_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_everg_byspecies<-aggregate(importance_data_everg_binaryWNV$relative_importance~importance_data_everg_binaryWNV$y+importance_data_everg_binaryWNV$lags, FUN=max)
# 
# importance_abund_evergeverg_bylag<-aggregate(importance_data_evergeverg_abund$relative_importance~importance_data_evergeverg_abund$buffers, FUN=mean)

# importance_everg_bybuffer<-aggregate(importance_data_everg_binaryWNV$relative_importance~importance_data_everg_binaryWNV$buffers, FUN=mean)

#importance_window_everg<-data.frame(rollapply(importance_everg, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_evergeverg_bybuffer<-data.frame(approx(x=importance_data_evergeverg_abund$buffers, y=log(importance_data_evergeverg_abund$relative_importance), n=10000))
# 
# interp_log_importance_everg_bybuffer<-data.frame(approx(x=importance_everg_bybuffer$'importance_data_everg_binaryWNV$buffers', y=log(importance_everg_bybuffer$'importance_data_everg_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/everg_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_evergeverg_bybuffer, aes(xmin=interp_log_importance_abund_evergeverg_bybuffer$x-0.45, xmax=interp_log_importance_abund_evergeverg_bybuffer$x+0.45, ymin=-1.5,ymax=2.5,fill=interp_log_importance_abund_evergeverg_bybuffer$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_line(data=mel_abund_everg_graph,aes(buffer,scale(aic_score)), color="grey") +
  geom_point(data=mel_abund_everg_graph,aes(buffer,scale(aic_score)), shape=19, stroke=1.5, size=3)+
  
  scale_x_continuous(name="Buffer Radius (m) for Evergreen Wetland") + 
  scale_y_continuous(name="Centered AIC Score")+ 
  scale_fill_distiller("Importance\nof Scale",palette = "YlOrRd", direction=1
                       ,breaks=c(-110, -60,  -10), labels=c("1.69e-48 Low", 4.54e-5,"8.76e-27 High")
  )+
  # scale_size_manual("Buffer Size (m) for\nForested Wetland Metric", values=sizes/600) +
  
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  #guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()




####melanura deciduous forested wetland NO CLIM for abund####

##melanura for decid
buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
mel_abund_decid_graph<-data.frame(buffer=as.numeric(),aic_score=as.numeric())
##for decid
for (i in 1:length(buffers)){
  
  model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(eval(as.name(paste("prop_PFO_Decid.x_",buffers[i], sep=""))), k=4)+s(prop_PEM_50, k=3)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+ti(prop_PFO_3000, prev_phdi_0, k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
  
  
  mel_abund_decid_graph<-rbind(mel_abund_decid_graph, data.frame(buffer=buffers[i],aic_score=model$aic))
  
}





#melanura decid calcuate a delta AIC column
mel_abund_decid_graph$delt_AIC<- mel_abund_decid_graph$aic_score - min(mel_abund_decid_graph$aic_score)


#melanura decid aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_decid_graph$weight_AIC<-exp( -0.5 * mel_abund_decid_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_decid_graph$delt_AIC))



#calcuating importance data for mel decid

importance_data_deciddecid_abund<-data.frame()


for (k in 1:nrow(mel_abund_decid_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_decid_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_decid_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_decid_graph[k, "buffer"],x=importance_val, y=species)
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_decid_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_decid_graph[k,]), "propmodels",values[k], sep="_")
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_deciddecid_abund<-rbind(importance_data_deciddecid_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_decid_byspecies<-aggregate(importance_data_decid_binaryWNV$relative_importance~importance_data_decid_binaryWNV$y+importance_data_decid_binaryWNV$lags, FUN=max)
# 
# importance_abund_deciddecid_bylag<-aggregate(importance_data_deciddecid_abund$relative_importance~importance_data_deciddecid_abund$buffers, FUN=mean)

# importance_decid_bybuffer<-aggregate(importance_data_decid_binaryWNV$relative_importance~importance_data_decid_binaryWNV$buffers, FUN=mean)

#importance_window_decid<-data.frame(rollapply(importance_decid, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_deciddecid_bybuffer<-data.frame(approx(x=importance_data_deciddecid_abund$buffers, y=log(importance_data_deciddecid_abund$relative_importance), n=10000))
# 
# interp_log_importance_decid_bybuffer<-data.frame(approx(x=importance_decid_bybuffer$'importance_data_decid_binaryWNV$buffers', y=log(importance_decid_bybuffer$'importance_data_decid_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/decid_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_deciddecid_bybuffer, aes(xmin=interp_log_importance_abund_deciddecid_bybuffer$x-0.45, xmax=interp_log_importance_abund_deciddecid_bybuffer$x+0.45, ymin=-2.5,ymax=2,fill=interp_log_importance_abund_deciddecid_bybuffer$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_line(data=mel_abund_decid_graph,aes(buffer,scale(aic_score)), color="grey") +
  geom_point(data=mel_abund_decid_graph,aes(buffer,scale(aic_score)), shape=19, stroke=1.5, size=3)+
  
  scale_x_continuous(name="Buffer Radius (m) for Deciduous Wetland") + 
  scale_y_continuous(name="Centered AIC Score")+ 
  scale_fill_distiller("Importance\nof Scale",palette = "YlOrRd", direction=1
                       ,breaks=c(-250, -130,  -10), labels=c("2.67e-109 Low", 3.48e-57,"4.54e-5  High")
  )+
  # scale_size_manual("Buffer Size (m) for\nForested Wetland Metric", values=sizes/600) +
  
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  #guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()








####melanura stream count for abund####

##melanura for streamCT
buffers<-c( 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
mel_abund_streamCT_graph<-data.frame(buffer=as.numeric(),aic_score=as.numeric())
##for streamCT
for (i in 1:length(buffers)){
  
  model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(eval(as.name(paste("PFO_mean_stream_cnt_",buffers[i],"_zeros", sep=""))), k=4)+s(prop_PEM_50, k=3)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PSS_3000,k=3)+s(prop_PFO_Decid.x_3000, k=4)+ti(prop_PFO_3000, prev_phdi_0, k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
  
  
  mel_abund_streamCT_graph<-rbind(mel_abund_streamCT_graph, data.frame(buffer=buffers[i],aic_score=model$aic))
  
}





#melanura streamCT calcuate a delta AIC column
mel_abund_streamCT_graph$delt_AIC<- mel_abund_streamCT_graph$aic_score - min(mel_abund_streamCT_graph$aic_score)


#melanura streamCT aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_streamCT_graph$weight_AIC<-exp( -0.5 * mel_abund_streamCT_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_streamCT_graph$delt_AIC))



#calcuating importance data for mel streamCT

importance_data_streamCTstreamCT_abund<-data.frame()


for (k in 1:nrow(mel_abund_streamCT_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_streamCT_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_streamCT_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_streamCT_graph[k, "buffer"],x=importance_val, y=species)
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_streamCT_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_streamCT_graph[k,]), "propmodels",values[k], sep="_")
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_streamCTstreamCT_abund<-rbind(importance_data_streamCTstreamCT_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_streamCT_byspecies<-aggregate(importance_data_streamCT_binaryWNV$relative_importance~importance_data_streamCT_binaryWNV$y+importance_data_streamCT_binaryWNV$lags, FUN=max)
# 
# importance_abund_streamCTstreamCT_bylag<-aggregate(importance_data_streamCTstreamCT_abund$relative_importance~importance_data_streamCTstreamCT_abund$buffers, FUN=mean)

# importance_streamCT_bybuffer<-aggregate(importance_data_streamCT_binaryWNV$relative_importance~importance_data_streamCT_binaryWNV$buffers, FUN=mean)

#importance_window_streamCT<-data.frame(rollapply(importance_streamCT, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_streamCTstreamCT_bybuffer<-data.frame(approx(x=importance_data_streamCTstreamCT_abund$buffers, y=log(importance_data_streamCTstreamCT_abund$relative_importance), n=10000))
# 
# interp_log_importance_streamCT_bybuffer<-data.frame(approx(x=importance_streamCT_bybuffer$'importance_data_streamCT_binaryWNV$buffers', y=log(importance_streamCT_bybuffer$'importance_data_streamCT_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c(200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/streamCT_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_streamCTstreamCT_bybuffer, aes(xmin=interp_log_importance_abund_streamCTstreamCT_bybuffer$x-0.45, xmax=interp_log_importance_abund_streamCTstreamCT_bybuffer$x+0.45, ymin=-2.5,ymax=2,fill=interp_log_importance_abund_streamCTstreamCT_bybuffer$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_line(data=mel_abund_streamCT_graph,aes(buffer,scale(aic_score)), color="grey") +
  geom_point(data=mel_abund_streamCT_graph,aes(buffer,scale(aic_score)), shape=19, stroke=1.5, size=3)+
  
  scale_x_continuous(name="Buffer Radius (m) for Stream Connectivity") + 
  scale_y_continuous(name="Centered AIC Score")+ 
  scale_fill_distiller("Importance\nof Scale",palette = "YlOrRd", direction=1
,breaks=c(-30, -17,  -5), labels=c("9.36e-14 Low", 4.14e-8,"0.007    High")
  )+
  # scale_size_manual("Buffer Size (m) for\nForested Wetland Metric", values=sizes/600) +
  
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  #guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()





####melanura impervious for abund####

##melanura for impv
buffers<-c( 50,100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
mel_abund_impv_graph<-data.frame(buffer=as.numeric(),aic_score=as.numeric())
##for impv
for (i in 1:length(buffers)){
  
  model<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(eval(as.name(paste("impvMN_",buffers[i], sep=""))), k=4)+s(prop_PEM_50, k=3)+s(prop_PFO_Ever.x_3000, k=3)+s(prop_PSS_3000,k=3)+s(prop_PFO_Decid.x_3000, k=4)+ti(prop_PFO_3000, prev_phdi_0, k=4)+s(prop_PFO_F_500,prev_phdi_8,k=4)+s(PFO_mean_stream_cnt_1500_zeros, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$PFO_mean_stream_cnt_1500_zeros<5,], select=T)
  
  
  mel_abund_impv_graph<-rbind(mel_abund_impv_graph, data.frame(buffer=buffers[i],aic_score=model$aic))
  
}





#melanura impv calcuate a delta AIC column
mel_abund_impv_graph$delt_AIC<- mel_abund_impv_graph$aic_score - min(mel_abund_impv_graph$aic_score)


#melanura impv aic weights
#caclulate AIC weights in order to get importance measurement http://brianomeara.info/tutorials/aic/ 
mel_abund_impv_graph$weight_AIC<-exp( -0.5 * mel_abund_impv_graph$delt_AIC)/sum(exp( -0.5 * mel_abund_impv_graph$delt_AIC))



#calcuating importance data for mel impv

importance_data_impvimpv_abund<-data.frame()


for (k in 1:nrow(mel_abund_impv_graph)){
  
  #calculating importance for each buffer/lag combo
  importance_val<-sum(mel_abund_impv_graph[k,"weight_AIC" ], na.rm = T)/sum(mel_abund_impv_graph["weight_AIC"])
  
  
  #creating a blank dataframe for the data from just buffer/buffer then filling it with the data
  little_data<-data.frame(z=as.numeric(),x=as.numeric(),y=as.numeric())
  species<-"mel"
  little_data<-data.frame(z=mel_abund_impv_graph[k, "buffer"],x=importance_val, y=species)
  
  #changing column names to match each covariate
  #the buffer/buffer label
  colnames(little_data)[1]<-"buffers"
  # colnames(little_data)[1]<-paste(colnames(mel_abund_impv_graph[k,]), "import",values[k], sep="_")
  #the importance measurement
  colnames(little_data)[2]<-"relative_importance"
  
  # colnames(little_data)[2]<-paste(colnames(mel_abund_impv_graph[k,]), "propmodels",values[k], sep="_")
  
  #binding the data for all the buffers/buffers for 1 coviarate together
  importance_data_impvimpv_abund<-rbind(importance_data_impvimpv_abund, little_data)
} 






#only considering the best importance value for any given climate lag (this removes all but the best buffer size)
# importance_binaryWNV_impv_byspecies<-aggregate(importance_data_impv_binaryWNV$relative_importance~importance_data_impv_binaryWNV$y+importance_data_impv_binaryWNV$lags, FUN=max)
# 
# importance_abund_impvimpv_bylag<-aggregate(importance_data_impvimpv_abund$relative_importance~importance_data_impvimpv_abund$buffers, FUN=mean)

# importance_impv_bybuffer<-aggregate(importance_data_impv_binaryWNV$relative_importance~importance_data_impv_binaryWNV$buffers, FUN=mean)

#importance_window_impv<-data.frame(rollapply(importance_impv, width=3, by=1,FUN=mean))

###plotting scale and importance##

#interpolation of ln importance
interp_log_importance_abund_impvimpv_bybuffer<-data.frame(approx(x=importance_data_impvimpv_abund$buffers, y=log(importance_data_impvimpv_abund$relative_importance), n=10000))
# 
# interp_log_importance_impv_bybuffer<-data.frame(approx(x=importance_impv_bybuffer$'importance_data_impv_binaryWNV$buffers', y=log(importance_impv_bybuffer$'importance_data_impv_binaryWNV$relative_importance'), n=10000))

#setting colors for lines and points
sizes<-c(50,100,200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_Figures_3_10_17/Melanura_abund/impv_importance_melabundance.png", width=3200, height=2400, res=300)

ggplot()+
  geom_rect(data=interp_log_importance_abund_impvimpv_bybuffer, aes(xmin=interp_log_importance_abund_impvimpv_bybuffer$x-0.45, xmax=interp_log_importance_abund_impvimpv_bybuffer$x+0.45, ymin=-2.5,ymax=2,fill=interp_log_importance_abund_impvimpv_bybuffer$y))+
  
  #scale_fill_distiller("ln(Relative\nImportance)",palette = "YlOrRd", direction=1 )+
  geom_line(data=mel_abund_impv_graph,aes(buffer,scale(aic_score)), color="grey") +
  geom_point(data=mel_abund_impv_graph,aes(buffer,scale(aic_score)), shape=19, stroke=1.5, size=3)+
  
  scale_x_continuous(name="Buffer Radius (m) for Impervious Surfaces") + 
  scale_y_continuous(name="Centered AIC Score")+ 
  scale_fill_distiller("Importance\nof Scale",palette = "YlOrRd", direction=1
                       ,breaks=c(-225, -122.5,  -20), labels=c("1.92e-98 Low", 6.29e-54,"2.06e-9  High")
  )+
  # scale_size_manual("Buffer Size (m) for\nForested Wetland Metric", values=sizes/600) +
  
  # scale_colour_manual(name  = "",values=c(sal_purple="#BA0CE8", rest_blue="#00B3FF", pip_green="#52E80C"), breaks=c("sal_purple", "rest_blue", "pip_green"),labels=c("Culex salinarius", "Culex restuans", "Culex pipiens"))+
  #guides( fill=guide_colorbar(order = 1),  size = guide_legend(order = 2))+
  theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),axis.text=element_text(size=22),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title.x=element_text(size=22), axis.title.y=element_text(size=22), legend.text=element_text(face="italic",size=19), legend.title=element_text(size=22))
dev.off()



