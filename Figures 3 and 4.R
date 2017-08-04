library(mgcv)
library(itsadug)

#Generating GAM response curves for Figure 3: relationships between wetland characteristics and abundance.

#load data from location of local machine. Available from https://knb.ecoinformatics.org/#view/doi:10.5063/F17P8WGK
fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/EEE_Cs_melanura_8_4_17.csv", header=T)



#Abundance model #1 with 0 month hydrological wetness lag
fullmodel_melanura_abund<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_100.y, k=3)+s(prop_PFO_Ever.x_1000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_2000_zeros, k=4)+s(prop_PFO_F_3000,prev_phdi_1,k=4)+ti(prop_PFO_3000,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$prop_PFO_Ever.x_500<.3&fulldata$PFO_mean_stream_cnt_1500_zeros<5, ], select=T)



#Abundance model #2 with 8 month hydrological wetness lag
fullmodel_melanura_abund_1<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_100.y, k=3)+s(prop_PFO_Ever.x_1000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_2000_zeros, k=4)+s(prop_PFO_F_3000,prev_phdi_1,k=4)+ti(prop_PFO_2000,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$prop_PFO_Ever.x_500<.3&fulldata$PFO_mean_stream_cnt_1500_zeros<5&fulldata$prop_PFO_2000<.15, ], select=T)







# #Fig. 3a
plot_smooth(fullmodel_melanura_abund, view="prop_PFO_Decid_3000",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of deciduous\nforested wetland (3000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(0,15))
axis(2,las=1)
axis(1)
mtext("(a)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 13%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)


# #Fig. 3b
plot_smooth(fullmodel_melanura_abund, view="prop_PFO_Ever.x_1000",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of evergreen\nforested wetland (1000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(1,6))
axis(2,las=1)
axis(1)
mtext("(b)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 7%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)


# #Fig. 3c
plot_smooth(fullmodel_melanura_abund, view="prop_PEM_100.y",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of emergent\nwetland (100m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim = c(0,2.5), shade = T)
axis(2,las=1)
axis(1)
mtext("(c)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 1%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)


#Fig. 3d
plot_smooth(fullmodel_melanura_abund, view="prop_PSS_3000",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of scrub/shurb\nwetland (3000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(0,2.5), xlim=c(0,.05))
axis(2,las=1)
axis(1)
mtext("(d)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 2%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)

#Fig. 3e
plot_smooth(fullmodel_melanura_abund, view="PFO_mean_stream_cnt_2000_zeros",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Mean # stream connections\n(2000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(.5,2), xlim=c(.5,3))
axis(2,las=1)
axis(1)
mtext("(e)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 1%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)

#Fig. 3f
plot_smooth(fullmodel_melanura_abund, view="impvMN_200",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Mean land-cover imperviousness\n(200m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")),ylim=c(0,3), xlim=c(0,50))
axis(2,las=1)
axis(1)
mtext("(f)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 11%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)





#Fig 4a
fvisgam(fullmodel_melanura_abund, view=c("prop_PFO_F_3000", "prev_phdi_1"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (1 month lag)" ,xlab="Proportional area semi-permanent\n forested wetland (3000m buffer)", xlim=c(0,.0012))
mtext("(a)", side = 3, adj = 0, line = .5, cex=2)

#Fig 4b
fvisgam(fullmodel_melanura_abund, view=c("prop_PFO_3000","prev_phdi_0"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (0 month lag)" ,xlab="Proportional area of forested\n wetlands (3000m buffer)")
mtext("(b)", side = 3, adj = 0, line = .5, cex=2)

#Fig 4c
fvisgam(fullmodel_melanura_abund_1, view=c("prop_PFO_2000","prev_phdi_8"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (8 month lag)" ,xlab="Proportional area of forested\n wetlands (2000m buffer)")
mtext("(c)", side = 3, adj = 0, line = .5, cex=2)
