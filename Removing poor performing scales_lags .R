####Main Regression Analsysis
library(mgcv)
library(gamm4)
library(bbmle)
library(gam)
library(itsadug)
library(lubridate)
library(car)
library(ggplot2)

#fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/fulldataset_mle_3_11_16.csv", header=T)

#fulldata<-read.csv("C:/Users/FWL/Documents/Nick's Work/R Code/aggregated_restpip_4_26_16.csv", header=T)

#fulldata<-read.csv("C:/Users/FWL/Documents/Nick's Work/R Code/aggregated_restpip_5_27_16.csv", header=T)
# 
# fulldata<-read.csv("C:/Users/FWL/Documents/Nick's Work/R Code/fulldataset_mle_7_19_16.csv", header=T)

fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/EEE_fulldataset_mle_7_26_16.csv", header=T, stringsAsFactors = F)

fulldata$month_year<-as.Date(fulldata$month_year)

fulldata$EEE_MIR_Pres<-fulldata$IR
fulldata[fulldata$EEE_MIR_Pres>0&!is.na(fulldata$EEE_MIR_Pres),"EEE_MIR_Pres"]<-1
fulldata[fulldata$EEE_MIR_Pres==0&!is.na(fulldata$EEE_MIR_Pres),"EEE_MIR_Pres"]<-0

fulldata$EEE_MIR_Pres<-as.factor(fulldata$EEE_MIR_Pres)


base_model_aic<-gam(EEE_MIR_Pres~s(log(Mean_Abundance_perTrapnight+1), k=4)+s(month_num, k=5)+s(Prop_Gravid_Trap, k=5), select=T, family=binomial, data=fulldata)$aic


####for loops to do the analyses quickly####
####creating graphs for each combination of phdi month (lags and just month see hashtags) and semi buffer distance####
#monthphdi_semi<-data.frame(gcv_50=as.numeric(),gcv_100=as.numeric(),gcv_200=as.numeric() ,gcv_500=as.numeric(), gcv_750=as.numeric(),gcv_1000=as.numeric(), gcv_1500=as.numeric(),gcv_2000=as.numeric(),gcv_3000=as.numeric(),gcv_4000=as.numeric(),gcv_5000=as.numeric(),phdi=as.character(), stringsAsFactors = F)

#monthphdi_semi<-data.frame(R2_50=as.numeric(),R2_100=as.numeric(),R2_200=as.numeric() ,R2_500=as.numeric(), R2_750=as.numeric(),R2_1000=as.numeric(), R2_1500=as.numeric(),R2_2000=as.numeric(),month=as.character(), stringsAsFactors = F)


monthphdi_semi<-data.frame(buffer=as.numeric(),phdi_lag=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (j in 1:11){
  for (i in 0:12){
    
    #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
    
    data<-fulldata[,c(paste("prev_phdi_",i,sep=""),paste("prop_semi_",buffers[j],"_ring", sep=""),"EEE_MIR_Pres","Mean_Abundance_perTrapnight", "month_num", "Prop_Gravid_Trap")]
    
    IR<-data[,3]
    semi<-data[,2]
    phdi<-data[,1]
    month_num<-data[,5]
    abundance<-data[,4]
    gravid<-data[,6]
    
    
    
    semi_month<-gam(IR~s(log(abundance+1), k=4)+ti(phdi,semi, k=4)+s(month_num, k=5)+s(gravid, k=5), select=T, family=binomial)
    
    
    monthphdi_semi<-rbind(monthphdi_semi, data.frame(buffer=buffers[j],phdi_lag=i,gcv_score=semi_month$aic))
    
    
  }
}



ggplot(monthphdi_semi, aes(x = phdi_lag, y = gcv_score)) + geom_point(aes(size = as.factor(buffer), color=as.factor(buffer)), shape=1, stroke=1)+ theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_continuous(name="PHDI Lag", breaks=seq(1,24,1))+ scale_y_continuous(name="GCV Score")+ scale_size_manual("Buffer Size", values=buffers/750) + scale_color_discrete(name="Buffer Size")+guides(color= guide_legend(), size=guide_legend())


#identifying minima for phdi and semi buffer. The first part figures out if the gcv score is a minima for a lag within all 50m buffers, all 100m buffers etc., whereas the second part figures out whether it's a minima for a buffer size among all the 1 lags, 2 lags etc. It has to be a minima for both its buffer and its lag to be included.

lags<-c(0,1,2,3,4,5,6,7,8,9,10,11,12)
semi_selected<-data.frame(buffers=as.character(), lags=as.character())
for (i in 1:13){
  for (j in 1:length(buffers)){
    
    
    if((monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i,"gcv_score"]<base_model_aic & 
        (monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i,"gcv_score"]<monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i+1,"gcv_score"] | is.na(monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i+1,"gcv_score"])) & 
        ((length(monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i-1,"gcv_score"])==0) | monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i,"gcv_score"]<sum(monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i-1,"gcv_score"])))
       
       &
       
       (monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j,"gcv_score"]<base_model_aic & 
        (monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j,"gcv_score"]<monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j+1,"gcv_score"] | is.na(monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j+1,"gcv_score"])) & 
        ((length(monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j-1,"gcv_score"])==0) | monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j,"gcv_score"]<sum(monthphdi_semi[monthphdi_semi$phdi_lag==lags[i],][j-1,"gcv_score"])))
       
    ) {
      
      buffer_size<-paste("prop_semi_",monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i,"buffer"], "_ring", sep="")
      lag_num<-paste("prev_phdi_", monthphdi_semi[monthphdi_semi$buffer==buffers[j],][i,"phdi_lag"], sep="")
      
      semi_selected<-rbind(semi_selected,data.frame(buffer=as.character(buffer_size), lags=as.character(lag_num)))
      
    }
  }
}




####using months rather than lags

# months<-c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
# 
# 
# #renaming columsn for current year so they can be included more easily in analysis
# for (i in 440:451){
#   colnames(fulldata)[i]<-paste("phdi","current", months[i-439], sep="_")
# }

##renaming columsn for previous year so they can be included more easily in analysis
# for (i in 452:463){
#   colnames(fulldata)[i]<-paste("phdi","previous", months[i-451], sep="_")
# }
# 
# 
# yearphdi_semi<-data.frame(buffer=as.numeric(),phdi_month=as.character(),gcv_score=as.numeric(), stringsAsFactors = F)
# buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
# 
# month_year_combos<-c("current_Jan","current_Feb", "current_Mar", "current_Apr","current_May"  , "current_Jun","current_Jul", "current_Aug","current_Sep","current_Oct" ,"current_Nov" ,"current_Dec", "previous_Jan", "previous_Feb" , "previous_Mar" , "previous_Apr","previous_May", "previous_Jun" ,"previous_Jul" , "previous_Aug", "previous_Sep", "previous_Oct","previous_Nov" , "previous_Dec")
# 
# for (j in 1:11){
#   for (i in 1:24){
#     
#     #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
#     
#     data<-fulldata[,c(paste("phdi_",month_year_combos[i],sep=""),paste("prop_semi_",buffers[j],"_ring", sep=""),"EEE_MIR_Pres","Mean_Abundance_perPool", "month_num")]
#     
#     IR<-data[,3]
#     semi<-data[,2]
#     phdi<-data[,1]
#     month_num<-data[,5]
#     abundance<-data[,4]
#     
#     
#     
#     semi_month<-gam(IR~s(abundance, k=4)+ti(phdi,semi, k=4)+s(month_num, k=5), select=T, family=binomial)
#     
#     
#     yearphdi_semi<-rbind(yearphdi_semi, data.frame(buffer=buffers[j],phdi_month=month_year_combos[i],gcv_score=semi_month$aic))
#     
#     
#   }
# }
# 
# 
# ggplot(yearphdi_semi, aes(x = as.numeric(phdi_month), y = gcv_score)) + geom_point(aes(size = as.factor(buffer), color=as.factor(buffer)), shape=1, stroke=1)+ theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_continuous(name="PHDI Lag", breaks=seq(1,24,1))+ scale_y_continuous(name="GCV Score")+ scale_size_manual("Buffer Size", values=sizes/750) + scale_color_discrete(name="Buffer Size")+guides(color= guide_legend(), size=guide_legend())
# 




####identifying best buffer for average size####
size_comparison<-data.frame(buffer=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
#buffer_size<-data.frame(R2_50=as.numeric(),R2_100=as.numeric(),R2_200=as.numeric() ,R2_500=as.numeric(), R2_750=as.numeric(),R2_1000=as.numeric(), R2_1500=as.numeric(),R2_2000=as.numeric(),month=as.character(), stringsAsFactors = F)

buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (j in 1:length(buffers)){
  
  #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
  
  data<-fulldata[,c(paste("ave_size_",buffers[j],"_zeros", sep=""),"EEE_MIR_Pres","Mean_Abundance_perTrapnight", "month_num", "Prop_Gravid_Trap")]
  
  IR<-data[,2]
  size<-data[,1]
  abundance<-data[,3]
  month<-data[,4]
  gravid<-data[,5]
  
  
  size_model<-gam(IR~s(size,k=4)+s(log(abundance+1),k=4)+s(month, k=5)+s(gravid,k=5), select=T, family=binomial)
  
  size_comparison<-rbind(size_comparison, data.frame(buffer=buffers[j],gcv_score=size_model$aic))
  
  #pvisgam(semi_month, view=c("semi", "month"), plot.type="contour", color="heat", ylab=paste("current_phdi_",month.abb[i],sep=""), xlab=paste("prop_semi_",buffers[j], sep=""), main=summary(semi_month)$s.table[1,4])
  
}


plot(size_comparison[,2]~size_comparison[,1], type="o")

###identifying local minima in plot
for (i in 1:length(size_comparison$gcv_score)){
  if(size_comparison[i,"gcv_score"]<base_model_aic & 
     (size_comparison[i,"gcv_score"]<size_comparison[i+1,"gcv_score"]| is.na(size_comparison[i+1,"gcv_score"])) & 
     ((length(size_comparison[i-1,"gcv_score"])==0) | size_comparison[i,"gcv_score"]<sum(size_comparison[i-1,"gcv_score"]))) {
    
    print(size_comparison$buffer[i])
  }
  
}


impv_comparison<-data.frame(buffer=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
#buffer_impv<-data.frame(R2_50=as.numeric(),R2_100=as.numeric(),R2_200=as.numeric() ,R2_500=as.numeric(), R2_750=as.numeric(),R2_1000=as.numeric(), R2_1500=as.numeric(),R2_2000=as.numeric(),month=as.character(), stringsAsFactors = F)



buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (j in 1:length(buffers)){
  
  #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
  
  data<-fulldata[,c(paste("impvMN_",buffers[j], sep=""),"EEE_MIR_Pres","Mean_Abundance_perTrapnight", "month_num", "Prop_Gravid_Trap")]
  
  IR<-data[,2]
  impv<-data[,1]
  abundance<-data[,3]
  month<-data[,4]
  gravid<-data[,5]
  
  impv_model<-gam(IR~s(impv,k=4)+s(log(abundance+1),k=4)+s(month, k=5)+s(gravid,k=5), select=T, family=binomial)
  
  impv_comparison<-rbind(impv_comparison, data.frame(buffer=buffers[j],gcv_score=impv_model$aic))
  
  #pvisgam(semi_month, view=c("semi", "month"), plot.type="contour", color="heat", ylab=paste("current_phdi_",month.abb[i],sep=""), xlab=paste("prop_semi_",buffers[j], sep=""), main=summary(semi_month)$s.table[1,4])
  
}



plot(impv_comparison[,2]~impv_comparison[,1], type="o")

###identifying local minima in plot
for (i in 1:length(impv_comparison$gcv_score)){
  if(impv_comparison[i,"gcv_score"]<base_model_aic & 
     (impv_comparison[i,"gcv_score"]<impv_comparison[i+1,"gcv_score"]| is.na(impv_comparison[i+1,"gcv_score"])) & 
     ((length(impv_comparison[i-1,"gcv_score"])==0) | impv_comparison[i,"gcv_score"]<sum(impv_comparison[i-1,"gcv_score"]))) {
    
    print(impv_comparison$buffer[i])
  }
  
}



####climate by month and year on it's own
# clim_comparison<-data.frame(lag=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
# #buffer_clim<-data.frame(R2_50=as.numeric(),R2_100=as.numeric(),R2_200=as.numeric() ,R2_500=as.numeric(), R2_750=as.numeric(),R2_1000=as.numeric(), R2_1500=as.numeric(),R2_2000=as.numeric(),month=as.character(), stringsAsFactors = F)
# 
# month_year_combos<-c("current_Jan","current_Feb", "current_Mar", "current_Apr","current_May"  , "current_Jun","current_Jul", "current_Aug","current_Sep","current_Oct" ,"current_Nov" ,"current_Dec", "previous_Jan", "previous_Feb" , "previous_Mar" , "previous_Apr","previous_May", "previous_Jun" ,"previous_Jul" , "previous_Aug", "previous_Sep", "previous_Oct","previous_Nov" , "previous_Dec")
# 
# 
# for (j in 1:24){
#   
#   #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
#   
#   data<-fulldata[,c(paste("phdi_",month_year_combos[j], sep=""),"EEE_MIR_Pres","Mean_Abundance_perPool", "month_num")]
#   
#   IR<-data[,2]
#   clim<-data[,1]
#   abundance<-data[,3]
#   month<-data[,4]
#   
#   
#   clim_model<-gam(IR~s(clim,k=4)+s(abundance,k=4)+s(month, k=5), select=T, family=binomial)
#   
#   clim_comparison<-rbind(clim_comparison, data.frame(lag=month_year_combos[j],gcv_score=clim_model$aic))
#   
#   #pvisgam(semi_month, view=c("semi", "month"), plot.type="contour", color="heat", ylab=paste("current_phdi_",month.abb[i],sep=""), xlab=paste("prop_semi_",buffers[j], sep=""), main=summary(semi_month)$s.table[1,4])
#   
# }
# 
# 
# 
# plot(clim_comparison[,2]~as.numeric(clim_comparison[,1]), type="o")
# 
# 
# 
# 
# ###identifying local minima in plot
# for (i in 1:length(clim_comparison$gcv_score)){
#   if(clim_comparison[i,"gcv_score"]<base_model_aic & 
#      (clim_comparison[i,"gcv_score"]<clim_comparison[i+1,"gcv_score"]| is.na(clim_comparison[i+1,"gcv_score"])) & 
#      ((length(clim_comparison[i-1,"gcv_score"])==0) | clim_comparison[i,"gcv_score"]<sum(clim_comparison[i-1,"gcv_score"]))) {
#     
#     print(clim_comparison$lag[i])
#   }
#   
# }
# 


####


####identifying best buffer for average clim1####
clim1_comparison<-data.frame(lag=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
#buffer_clim1<-data.frame(R2_50=as.numeric(),R2_100=as.numeric(),R2_200=as.numeric() ,R2_500=as.numeric(), R2_750=as.numeric(),R2_1000=as.numeric(), R2_1500=as.numeric(),R2_2000=as.numeric(),month=as.character(), stringsAsFactors = F)

buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (j in 0:12){
  
  #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
  
  data<-fulldata[,c(paste("prev_phdi_",j, sep=""),"EEE_MIR_Pres","Mean_Abundance_perTrapnight", "month_num", "Prop_Gravid_Trap")]
  
  IR<-data[,2]
  clim1<-data[,1]
  abundance<-data[,3]
  month<-data[,4]
  gravid<-data[,5]
  
  
  clim1_model<-gam(IR~s(clim1,k=4)+s(log(abundance+1),k=4)+s(month, k=5)+s(gravid, k=5), select=T, family=binomial)
  
  clim1_comparison<-rbind(clim1_comparison, data.frame(lag=j,gcv_score=clim1_model$aic))
  
  #pvisgam(semi_month, view=c("semi", "month"), plot.type="contour", color="heat", ylab=paste("current_phdi_",month.abb[i],sep=""), xlab=paste("prop_semi_",buffers[j], sep=""), main=summary(semi_month)$s.table[1,4])
  
}


plot(clim1_comparison[,2]~clim1_comparison[,1], type="o")

###identifying local minima in plot
for (i in 1:length(clim1_comparison$gcv_score)){
  if(clim1_comparison[i,"gcv_score"]<base_model_aic & 
     (clim1_comparison[i,"gcv_score"]<clim1_comparison[i+1,"gcv_score"]| is.na(clim1_comparison[i+1,"gcv_score"])) & 
     ((length(clim1_comparison[i-1,"gcv_score"])==0) | clim1_comparison[i,"gcv_score"]<sum(clim1_comparison[i-1,"gcv_score"]))) {
    
    print(clim1_comparison$lag[i])
  }
  
}







# ###temp analysis on its own
# temp_comparison<-data.frame(buffer=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
# #buffer_temp<-data.frame(R2_50=as.numeric(),R2_100=as.numeric(),R2_200=as.numeric() ,R2_500=as.numeric(), R2_750=as.numeric(),R2_1000=as.numeric(), R2_1500=as.numeric(),R2_2000=as.numeric(),month=as.character(), stringsAsFactors = F)
# 
# 
# 
# buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
# 
# 
# for (j in 1:length(buffers)){
#   
#   #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_semi_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
#   
#   data<-fulldata[,c(paste("prop_temp_",buffers[j], sep=""),"EEE_MIR_Pres","Mean_Abundance_perPool", "month_num", "Prop_Gravid_Trap")]
#   
#   IR<-data[,2]
#   temp<-data[,1]
#   abundance<-data[,3]
#   month<-data[,4]
#   gravid<-data[,5]
#   
#   
#   temp_model<-gam(IR~s(temp,k=4)+s(abundance,k=4)+s(month, k=5)+s(gravid, k=5), select=T, family=binomial)
#   
#   temp_comparison<-rbind(temp_comparison, data.frame(buffer=buffers[j],gcv_score=temp_model$aic))
#   
#   #pvisgam(semi_month, view=c("semi", "month"), plot.type="contour", color="heat", ylab=paste("current_phdi_",month.abb[i],sep=""), xlab=paste("prop_semi_",buffers[j], sep=""), main=summary(semi_month)$s.table[1,4])
#   
# }
# 
# 
# 
# plot(temp_comparison[,2]~temp_comparison[,1], type="o")
# 
# ###identifying local minima in plot
# for (i in 1:length(temp_comparison$gcv_score)){
#   if(temp_comparison[i,"gcv_score"]<base_model_aic & 
#      (temp_comparison[i,"gcv_score"]<temp_comparison[i+1,"gcv_score"]| is.na(temp_comparison[i+1,"gcv_score"])) & 
#      ((length(temp_comparison[i-1,"gcv_score"])==0) | temp_comparison[i,"gcv_score"]<sum(temp_comparison[i-1,"gcv_score"]))) {
#     
#     print(temp_comparison$buffer[i])
#   }
#   
# }



# ####month by temp interaction plot
# monthphdi_temp<-data.frame(buffer=as.numeric(),phdi_lag=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
# buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
# 
# 
# for (j in 1:11){
#   for (i in 0:12){
#     
#     #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_temp_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
#     
#     data<-fulldata[,c(paste("prev_phdi_",i,sep=""),paste("prop_temp_",buffers[j], sep=""),"EEE_MIR_Pres","Mean_Abundance_perPool", "month_num", "Prop_Gravid_Trap")]
#     
#     IR<-data[,3]
#     temp<-data[,2]
#     phdi<-data[,1]
#     month_num<-data[,5]
#     abundance<-data[,4]
#     gravid<-data[,6]
#     
#     
#     
#     temp_month<-gam(IR~s(abundance, k=4)+ti(phdi,temp, k=4)+s(month_num, k=5)+s(gravid, k=5), select=T, family=binomial)
#     
#     
#     monthphdi_temp<-rbind(monthphdi_temp, data.frame(buffer=buffers[j],phdi_lag=i,gcv_score=temp_month$aic))
#     
#     
#   }
# }
# 
# 
# ggplot(monthphdi_temp, aes(x = phdi_lag, y = gcv_score)) + geom_point(aes(size = as.factor(buffer), color=as.factor(buffer)), shape=1, stroke=1)+ theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_continuous(name="PHDI Lag", breaks=seq(1,24,1))+ scale_y_continuous(name="GCV Score")+ scale_size_manual("Buffer Size", values=buffers/750) + scale_color_discrete(name="Buffer Size")+guides(color= guide_legend(), size=guide_legend())
# 
# 
# 
# #identifying minima for phdi and temp buffer. The first part figures out if the gcv score is a minima for a lag within all 50m buffers, all 100m buffers etc., whereas the second part figures out whether it's a minima for a buffer size among all the 1 lags, 2 lags etc. It has to be a minima for both its buffer and its lag to be included.
# 
# lags<-c(0,1,2,3,4,5,6,7,8,9,10,11,12)
# temp_selected<-data.frame(buffers=as.character(), lags=as.character())
# for (i in 1:13){
#   for (j in 1:length(buffers)){
#     
#     
#     if((monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"gcv_score"]<base_model_aic & 
#         (monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"gcv_score"]<monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i+1,"gcv_score"] | is.na(monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i+1,"gcv_score"])) & 
#         ((length(monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i-1,"gcv_score"])==0) | monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"gcv_score"]<sum(monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i-1,"gcv_score"])))
#        
#        &
#        
#        (monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j,"gcv_score"]<base_model_aic & 
#         (monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j,"gcv_score"]<monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j+1,"gcv_score"] | is.na(monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j+1,"gcv_score"])) & 
#         ((length(monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j-1,"gcv_score"])==0) | monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j,"gcv_score"]<sum(monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j-1,"gcv_score"])))
#        
#     ) {
#       
#       buffer_size<-paste("prop_temp_",monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"buffer"], sep="")
#       lag_num<-paste("prev_phdi_", monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"phdi_lag"], sep="")
#       
#       temp_selected<-rbind(temp_selected,data.frame(buffer=as.character(buffer_size), lags=as.character(lag_num)))
#       
#     }
#   }
# }
# 






####month by temp WITH RINGS interaction plot
monthphdi_temp<-data.frame(buffer=as.numeric(),phdi_lag=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)


for (j in 1:11){
  for (i in 0:12){
    
    #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_temp_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
    
    data<-fulldata[,c(paste("prev_phdi_",i,sep=""),paste("prop_temp_",buffers[j],"_ring", sep=""),"EEE_MIR_Pres","Mean_Abundance_perTrapnight", "month_num", "Prop_Gravid_Trap")]
    
    IR<-data[,3]
    temp<-data[,2]
    phdi<-data[,1]
    month_num<-data[,5]
    abundance<-data[,4]
    gravid<-data[,6]
    
    
    
    temp_month<-gam(IR~s(log(abundance+1), k=4)+ti(phdi,temp, k=4)+s(month_num, k=5)+s(gravid, k=5), select=T, family=binomial)
    
    
    monthphdi_temp<-rbind(monthphdi_temp, data.frame(buffer=buffers[j],phdi_lag=i,gcv_score=temp_month$aic))
    
    
  }
}


ggplot(monthphdi_temp, aes(x = phdi_lag, y = gcv_score)) + geom_point(aes(size = as.factor(buffer), color=as.factor(buffer)), shape=1, stroke=1)+ theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_continuous(name="PHDI Lag", breaks=seq(1,24,1))+ scale_y_continuous(name="GCV Score")+ scale_size_manual("Buffer Size", values=buffers/750) + scale_color_discrete(name="Buffer Size")+guides(color= guide_legend(), size=guide_legend())



#identifying minima for phdi and temp buffer. The first part figures out if the gcv score is a minima for a lag within all 50m buffers, all 100m buffers etc., whereas the second part figures out whether it's a minima for a buffer size among all the 1 lags, 2 lags etc. It has to be a minima for both its buffer and its lag to be included.

lags<-c(0,1,2,3,4,5,6,7,8,9,10,11,12)
temp_selected<-data.frame(buffers=as.character(), lags=as.character())
for (i in 1:13){
  for (j in 1:length(buffers)){
    
    
    if((monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"gcv_score"]<base_model_aic & 
        (monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"gcv_score"]<monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i+1,"gcv_score"] | is.na(monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i+1,"gcv_score"])) & 
        ((length(monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i-1,"gcv_score"])==0) | monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"gcv_score"]<sum(monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i-1,"gcv_score"])))
       
       &
       
       (monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j,"gcv_score"]<base_model_aic & 
        (monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j,"gcv_score"]<monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j+1,"gcv_score"] | is.na(monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j+1,"gcv_score"])) & 
        ((length(monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j-1,"gcv_score"])==0) | monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j,"gcv_score"]<sum(monthphdi_temp[monthphdi_temp$phdi_lag==lags[i],][j-1,"gcv_score"])))
       
    ) {
      
      buffer_size<-paste("prop_temp_",monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"buffer"], "_ring",sep="")
      lag_num<-paste("prev_phdi_", monthphdi_temp[monthphdi_temp$buffer==buffers[j],][i,"phdi_lag"], sep="")
      
      temp_selected<-rbind(temp_selected,data.frame(buffer=as.character(buffer_size), lags=as.character(lag_num)))
      
    }
  }
}



# ####month by total interaction plot
# monthphdi_total<-data.frame(buffer=as.numeric(),phdi_lag=as.numeric(),gcv_score=as.numeric(), stringsAsFactors = F)
# buffers<-c(50, 100, 200, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000)
# 
# 
# for (j in 1:11){
#   for (i in 0:12){
#     
#     #data<-fulldata[,c(paste("current_phdi_",month.abb[i],sep=""),paste("prop_totalwetland_",buffers[j], sep=""),"EEE_MIR_Pres", "LONG_DD", "LAT_DD")]
#     
#     data<-fulldata[,c(paste("prev_phdi_",i,sep=""),paste("prop_totalwetland_",buffers[j], sep=""),"EEE_MIR_Pres","Mean_Abundance_perPool", "month_num")]
#     
#     IR<-data[,3]
#     total<-data[,2]
#     phdi<-data[,1]
#     month_num<-data[,5]
#     abundance<-data[,4]
#     
#     
#     
#     total_month<-gam(IR~s(abundance, k=4)+ti(phdi,total, k=4)+s(month_num, k=5), select=T, family=binomial)
#     
#     
#     monthphdi_total<-rbind(monthphdi_total, data.frame(buffer=buffers[j],phdi_lag=i,gcv_score=total_month$aic))
#     
#     
#   }
# }
# 
# 
# ggplot(monthphdi_total, aes(x = phdi_lag, y = gcv_score)) + geom_point(aes(size = as.factor(buffer), color=as.factor(buffer)), shape=1, stroke=1)+ theme(axis.line.x = element_line(color = "black"),axis.line.y = element_line(color = "black"),panel.background=element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+ scale_x_continuous(name="PHDI Lag", breaks=seq(1,24,1))+ scale_y_continuous(name="GCV Score")+ scale_size_manual("Buffer Size", values=buffers/750) + scale_color_discrete(name="Buffer Size")+guides(color= guide_legend(), size=guide_legend())
# 
# 
# 
# #identifying minima for phdi and total buffer. The first part figures out if the gcv score is a minima for a lag within all 50m buffers, all 100m buffers etc., whereas the second part figures out whether it's a minima for a buffer size among all the 1 lags, 2 lags etc. It has to be a minima for both its buffer and its lag to be included.
# 
# lags<-c(0,1,2,3,4,5,6,7,8,9,10,11,12)
# total_selected<-data.frame(buffers=as.character(), lags=as.character())
# for (i in 1:13){
#   for (j in 1:length(buffers)){
#     
#     
#     if((monthphdi_total[monthphdi_total$buffer==buffers[j],][i,"gcv_score"]<base_model_aic & 
#         (monthphdi_total[monthphdi_total$buffer==buffers[j],][i,"gcv_score"]<monthphdi_total[monthphdi_total$buffer==buffers[j],][i+1,"gcv_score"] | is.na(monthphdi_total[monthphdi_total$buffer==buffers[j],][i+1,"gcv_score"])) & 
#         ((length(monthphdi_total[monthphdi_total$buffer==buffers[j],][i-1,"gcv_score"])==0) | monthphdi_total[monthphdi_total$buffer==buffers[j],][i,"gcv_score"]<sum(monthphdi_total[monthphdi_total$buffer==buffers[j],][i-1,"gcv_score"])))
#        
#        &
#        
#        (monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j,"gcv_score"]<base_model_aic & 
#         (monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j,"gcv_score"]<monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j+1,"gcv_score"] | is.na(monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j+1,"gcv_score"])) & 
#         ((length(monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j-1,"gcv_score"])==0) | monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j,"gcv_score"]<sum(monthphdi_total[monthphdi_total$phdi_lag==lags[i],][j-1,"gcv_score"])))
#        
#     ) {
#       
#       buffer_size<-paste("prop_totalwetland_",monthphdi_total[monthphdi_total$buffer==buffers[j],][i,"buffer"], sep="")
#       lag_num<-paste("prev_phdi_", monthphdi_total[monthphdi_total$buffer==buffers[j],][i,"phdi_lag"], sep="")
#       
#       total_selected<-rbind(total_selected,data.frame(buffer=as.character(buffer_size), lags=as.character(lag_num)))
#       
#     }
#   }
# }


##final model based on individual selection of covariates

fullmodel_IRcontinuous<-gam(EEE_MIR_Pres~s(log(Mean_Abundance_perTrapnight+1), k=4)+s(ave_size_4000_zeros, k=4)+ti(prop_semi_4000_ring, prev_phdi_3)+ti(prop_temp_5000, prev_phdi_9)+s(impvMN_3000, k=4)+s(prev_phdi_1, k=4)+s(month_num,k=5)+s(Prop_Gravid_Trap, k=5),data=fulldata,family=binomial, select=T)

summary(fullmodel_IRcontinuous)

plot(fullmodel_IRcontinuous, select=1, scale=0)
plot(fullmodel_IRcontinuous, select=2, scale=0)

plot(fullmodel_IRcontinuous, select=3, scale=0)
plot(fullmodel_IRcontinuous, select=4, scale=0)
plot(fullmodel_IRcontinuous, select=5, scale=0)
plot(fullmodel_IRcontinuous, select=6, scale=0)
plot(fullmodel_IRcontinuous, select=7, scale=0)
plot(fullmodel_IRcontinuous, select=8, scale=0)
vis.concurvity(fullmodel_IRcontinuous)
pvisgam(fullmodel_IRcontinuous, view=c("prop_semi_4000_ring", "prev_phdi_3"), plot.type="contour", color="heat", type="response", too.far=0.2)

pvisgam(fullmodel_IRcontinuous, view=c("prop_temp_5000", "prev_phdi_9"), plot.type="contour", color="heat", type="response", too.far=0.1)
