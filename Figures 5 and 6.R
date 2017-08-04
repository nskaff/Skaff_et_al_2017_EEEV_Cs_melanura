library(mgcv)
library(itsadug)

#Generating GAM response curves for Figure 5: relationships between wetland characteristics and EEEV infection.

#load data from location of local machine. Available from https://knb.ecoinformatics.org/#view/doi:10.5063/F17P8WGK
fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/EEE_Cs_melanura_8_4_17.csv", header=T)


fullmodel_mel_binaryEEE<-gam(EEE_MIR_PRES~Prop_Gravid_Trap+s(ln_Mean_Abundance_perTrapnight, k=3)+s(prop_PEM_2000.y, k=4)+s(prop_PSS_5000,k=4)+s(prop_PFO_Ever.x_1000, k=4)+s(prop_PFO_Decid.x_5000, k=4)+s(ave_size_PFO_1500_zeros, k=4)+s(PFO_F_Relative_Temp_500,prev_phdi_10,k=4)+s(impvMN_200, k=4)+s(month_num, k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11,],family=binomial, select=T)

fullmodel_mel_binaryEEE_1<-gam(EEE_MIR_PRES~Prop_Gravid_Trap+s(ln_Mean_Abundance_perTrapnight, k=3)+s(prop_PEM_2000.y, k=4)+s(prop_PSS_5000,k=4)+s(prop_PFO_Ever.x_1000, k=4)+s(prop_PFO_Decid.x_5000, k=4)+s(ave_size_PFO_1500_zeros, k=4)+s(PFO_F_Relative_Temp_500,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num, k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11,],family=binomial, select=T)

#figure 5a
plot_smooth(fullmodel_mel_binaryEEE, view="prop_PFO_Decid.x_5000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of deciduous-forested\nwetland (5000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-10,-4), xlim=c(0,.1))
#something weird is happening and number isn't printing so have to add manually
axis(2, at=c(-10,-9,-8,-7,-6,-5,-4),label=c(-10,-9,-8,-7,-6,-5,-4),col="black",lwd=1, line=0, las=1)
axis(1)

mtext("(a)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 3%\np=0.005", side = 3,adj = 0.9, line = -1, cex=1.3)


#Figure 5b
plot_smooth(fullmodel_mel_binaryEEE, view="prop_PFO_Ever.x_1000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of evergreen-forested\nwetland (1000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,-2))

#something weird is happening and number isn't printing so have to add manually
axis(2, at=c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3,-2),label=c(-12,NA,-10,NA,-8,NA,-6,NA,-4, NA,-2),col="black",lwd=1, line=0, las=1 )
axis(1)

mtext("(b)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 0%\np=0.26", side = 3,adj = 0.9, line = -1, cex=1.3)


#figure 5c
plot_smooth(fullmodel_mel_binaryEEE, view="prop_PEM_2000.y",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of emergent\nwetland (2000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-25,-5))

axis(2, at=c(-25,-20,-15,-10,-5),label=c(-25,-20,-15,-10,-5),col="black",lwd=1, line=0, las=1 )
axis(1)


mtext("(c)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 1%\np=0.01", side = 3,adj = 0.9, line = -1, cex=1.3)



#figure 5d
plot_smooth(fullmodel_mel_binaryEEE, view="prop_PSS_5000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of scrub/shurb\nwetland (5000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,-2), xlim=c(0,.03))

#something weird is happening and number isn't printing so have to add manually
axis(2, at=c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3,-2),label=c(-12,NA,-10,NA,-8,NA,-6,NA,-4, NA,-2),col="black",lwd=1, line=0, las=1 )
axis(1)


mtext("(d)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 0%\np=0.74", side = 3,adj = 0.9, line = -1, cex=1.3)



#figure 5e
plot_smooth(fullmodel_mel_binaryEEE, view="ave_size_PFO_1500_zeros",col='red',lwd=3, rug=F, main="", hide.label = T,h0=NA ,axes=F, xlab="Average forested wetland size\n(km2; 1500m Buffer)", ylab="Log odds of EEEV positive pool", ylim=c(-10,-5))


axis(2, at=c(-10,-9,-8,-7,-6,-5, -4),label=c(-10,-9,-8,-7,-6,-5,-4),col="black",lwd=1, line=0, las=1 )
#axis(1)
axis(1, at=(c(0,20,40,60,80,100)*1000),label=c(0,20,40,60,80,100) )

mtext("(e)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 1%\np=0.078", side = 3,adj = 0.9, line = -1, cex=1.3)



#figure 5f
plot_smooth(fullmodel_mel_binaryEEE, view="impvMN_200",col='red',lwd=3, rug=F, main="", hide.label = T,h0=NA ,axes=F, xlab="Mean land-cover imperviousness\n(200m Buffer)", ylab="Log odds of EEEV positive pool", ylim=c(-14, -4))

axis(2, at=c(-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4),label=c(-14,NA,-12,NA,-10,NA,-8,NA,-6,NA,-4),col="black",lwd=1, line=0, las=1 )
axis(1)


mtext("(f)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 1%\np=0.1", side = 3,adj = 0.9, line = -1, cex=1.3)



#figure 5g
plot_smooth(fullmodel_mel_binaryEEE, view="ln_Mean_Abundance_perTrapnight",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab=substitute(paste("Mean ",italic("Cs. melanura "),"capture per night", sep=" ")), ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,0))

axis(2, at=c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3,-2,-1,0),label=c(-12,NA,-10,NA,-8,NA,-6,NA,-4, NA,-2,NA,0),col="black",lwd=1, line=0, las=1 )
#axis(1)

axis(1, at=log(c(0,1,5,15,50,150,500, 1000)+1), label=c(0,1,5,15,50,150,500,NA))


mtext("(g)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 15%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)



#figure 6a
fvisgam(fullmodel_mel_binaryEEE, view=c("PFO_F_Relative_Temp_500","prev_phdi_10"), color="topo",  hide.label = T, add.color.legend = T,too.far=.4,main="",ylab="PHDI (10 month lag)" ,xlab="Relative area of semi-permanent forested\nwetland (500m buffer)")
mtext("(a)", side = 3, adj = 0, line = .5, cex=2)


#figure 6b
fvisgam(fullmodel_mel_binaryEEE_1, view=c("PFO_F_Relative_Temp_500","prev_phdi_0"), color="topo",  hide.label = T, add.color.legend = T,too.far=.4,main="",ylab="PHDI (0 month lag)" ,xlab="Relative area of semi-permanent forested\nwetland (500m buffer)")
mtext("(b)", side = 3, adj = 0, line = .5, cex=2)


