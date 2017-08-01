###doing gam analysis for IR positive but after doing a repeated regression analysis to limit the number of potential combinations of covariates


library(parallel)

library(mgcv)

library(data.table)


##using expand grid to generate all combinations of models
# 
# fulldata<-read.csv("C:/Users/FWL/Documents/Nick's Work/R Code/aggregated_restpip_4_26_16.csv", header=T)

fulldata<-read.csv("C:/Users/FWL/Documents/Nick's Work/R Code/aggregated_restpip_5_27_16.csv", header=T)

buffers<-c(50, 100, 200, 500, 750,1000, 1500, 2000, 3000, 4000, 5000)

static<-data.frame(fulldata[,c("MYSiteCode","month_year","logit_IR","month_num","Mean_Abundance_perTrapnight",  "Prop_Gravid_Trap")])

# clim1<-data.frame(fulldata[,c("current_phdi_Jan","current_phdi_Apr","current_phdi_Jun","current_phdi_Aug", "current_phdi_Oct", "previous_phdi_Feb", "previous_phdi_Jun", "previous_phdi_Nov")])

clim1<-data.frame(fulldata[,c("prev_phdi_1","prev_phdi_4","prev_phdi_7","prev_phdi_11")])

size<-data.frame(fulldata[,c("ave_size_50_zeros","ave_size_200_zeros","ave_size_4000_zeros")])



impv<-data.frame(fulldata[,c("impvMN_1000","impvMN_3000","impvMN_5000")])

#total<-data.frame(fulldata[,c("prop_totalwetland_100","prop_totalwetland_4000")])

#temp<-data.frame(fulldata[,c("prop_temp_50","prop_temp_4000")])



######################
##expand grid creates all combinations of variables from each of the subsets above. This is the basis for each model combinations. I'm skipping the first 6 columns becasue they are stuff like abundance and IR
variables<-expand.grid(size=colnames(size[1:length(size)]), impv=colnames(impv[,1:length(impv)]),  paste(semi_selected[,1],semi_selected[,2], sep=" ")
                       ,clim1=colnames(clim1[,1:length(clim1)])
                       #,temp=colnames(temp[,1:length(temp)])
                       #,total=colnames(total[,1:length(total)])
                       #, paste(total_selected[,1],total_selected[,2], sep=" ")
                       , paste(temp_selected[,1],temp_selected[,2], sep=" ")
)

variables2 <- data.frame(do.call(rbind, strsplit(as.character(variables$Var3), ' ')),size=variables[,1],impv=variables[,2]
                         , clim1=variables[,4]
                        # ,temp=variables[,5]
                         #,total=variables[,5]
                         ,do.call(rbind, strsplit(as.character(variables$Var5), ' '))
                         
)

colnames(variables2)[1:2]<-c("semi", "phdi_semi")

#colnames(variables2)[6:7]<-c("total", "phdi_total")
colnames(variables2)[6:7]<-c("temp", "phdi_temp")

View(variables2)
##addign a column with the separator indicated so that it can be used in the paste below
#vars<-c(variables, sep="+")

##have to use do.call to do the paste on all rows and columns
#form1<-do.call(paste, vars)
##adding the response variable to the model
#form2<-paste("logit_IR~",form1, sep="")

###running the gam analysis on all the combinations of variables2
gam_analysis<-function(index){
  #simple the gam model with the formula put in. Index is just a vector of numbers from 1 to the last model number so that we do one model at a time. Need to convert the model to a formula with as.formula
  #model<-gam(as.formula(form2[index]),data=fulldata, select=T)
  
  model<-gam(WNV_MIR_Pres~s(log(Mean_Abundance_perTrapnight+1), k=4)+
               s(month_num, k=5)+s(Prop_Gravid_Trap, k=5)+
               #ti(eval(as.name(as.character(variables2[index, "semi"]))), k=4)+
               #ti(eval(as.name(as.character(variables2[index, "temp"]))), k=4)+
               s(eval(as.name(as.character(variables2[index, "clim1"]))), k=4)+
               s(eval(as.name(as.character(variables2[index, "size"]))), k=4)+
               s(eval(as.name(as.character(variables2[index, "impv"]))), k=4)+
               ti(eval(as.name(as.character(variables2[index, "semi"]))),
                  eval(as.name(as.character(variables2[index, "phdi_semi"]))))+
             #ti(eval(as.name(as.character(variables2[index, "total"]))),
                # eval(as.name(as.character(variables2[index, "phdi_total"]))), k=4)+
              ti(eval(as.name(as.character(variables2[index, "temp"]))),
               eval(as.name(as.character(variables2[index, "phdi_temp"])))) 
             ,data=fulldata[fulldata$Species_fix=="CULEX RESTUANS",], select=T, family=binomial)
  
  ##the below are extract the buffer size from each component of the model. The beginning part after gsub is just saying remove all the text from the name. The stuff after as.formula is taking the formula that was put into the model above, first extracting the part with the covariates, then extracting each individual covariate. Also, since there are multple numbers in each name, the last [[2]] is just taking the second one which is the buffer size
  semi_buff<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index,"semi"]))))
  #total_buff<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "total"]))))
  temp_buff<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "temp"]))))
  impv_buff<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "impv"]))))
  phdi_semi<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "phdi_semi"]))))
 # phdi_total<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "phdi_total"]))))
  phdi_temp<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "phdi_temp"]))))
  size_buff<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "size"]))))
  clim1_buff<-as.numeric(gsub("\\D", "",as.name(as.character(variables2[index, "clim1"]))))
  #create a row with all the important components as columns. This will later be bound together to create the final result
  return(data.frame(gcv.score=as.numeric(summary(model)$sp.criterion),aic.score=as.numeric(model$aic), semi_buffer=as.numeric(semi_buff), phdi_semi=as.numeric(phdi_semi), 
                     size_buff=as.numeric(size_buff),
                   # total_buffer=as.numeric(total_buff),
                   # phdi_total=as.numeric(phdi_total), 
                     temp_buffer=as.numeric(temp_buff),
                    phdi_temp=as.numeric(phdi_temp), 
                    impv_buff=as.numeric(impv_buff),
                clim1_buff=as.character(clim1_buff)
))
  
}




#how many cores are in the computer
detectCores()
#tell the program how many clusters to use
#c1<-makeCluster(3)
c1<-makeCluster(11)
#you have to export all the dataframes called by the parLapply/gam_analysis function to the clusters. not necessary with FORK cluster on mac

clusterExport(c1, c("variables2", "fulldata"))

#have to export the libraries that are used by the gam_analuysis function. Not necessary with FORK cluster on mac

clusterEvalQ(c1, library(mgcv))

#starting the timer
ptm <- proc.time()
#creating an index so that each component of the list of formulas will be run in the function
index<-1:nrow(variables2)
#using a parallel version of lapply, saying send index to the gam_analysis function
results<-parLapply(c1, index, gam_analysis)

#end timer
proc.time() - ptm
#turn off the clusters so it doesn't keep using tons of memory/processor
stopCluster(c1);gc()
closeAllConnections() 

#bind all the results columns together
gcv_table<-rbindlist(results)

write.csv(gcv_table, "C:/Users/fwl/Documents/Nick's Work/GCV_AUC_output_tables/restuans_WNV_binary_7_18_16/restuans_binary_tempSemiinteract_zeros_7_18_16.csv")
head(gcv_table[(order(gcv_table$aic.score)),] ,n=30)
