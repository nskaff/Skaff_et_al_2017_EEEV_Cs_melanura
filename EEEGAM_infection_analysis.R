
fulldata<-read.csv("~/Documents/Grad School/West Nile Virus/Data/Chapter 2 Datafiles/CT Regression Datasets/EEE_data_mle_3_10_17.csv", header=T)

library(mgcv)
library(lubridate)
library(itsadug)
load("~/Documents/Grad School/Edited R Packages/plot_error_grey.Rdata")
load("~/Documents/Grad School/Edited R Packages/plot_smooth_grey.Rdata")

fulldata$month_year<-as.Date(fulldata$month_year)
fulldata$month_num<-month(fulldata$month_year)
##maybe only use scales of figure out a way to remove outliers



fullmodel_mel_binaryEEE<-gam(EEE_MIR_PRES~Prop_Gravid_Trap+s(ln_Mean_Abundance_perTrapnight, k=3)+s(prop_PEM_2000.y, k=4)+s(prop_PSS_5000,k=4)+s(prop_PFO_Ever.x_1000, k=4)+s(prop_PFO_Decid.x_5000, k=4)+s(ave_size_PFO_1500_zeros, k=4)+s(PFO_F_Relative_Temp_500,prev_phdi_10,k=4)+s(impvMN_200, k=4)+s(month_num, k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11,],family=binomial, select=T)


summary(fullmodel_mel_binaryEEE)
gam.check(fullmodel_mel_binaryEEE)

plot(fullmodel_mel_binaryEEE)

fullmodel_mel_binaryEEE_1<-gam(EEE_MIR_PRES~Prop_Gravid_Trap+s(ln_Mean_Abundance_perTrapnight, k=3)+s(prop_PEM_2000.y, k=4)+s(prop_PSS_5000,k=4)+s(prop_PFO_Ever.x_1000, k=4)+s(prop_PFO_Decid.x_5000, k=4)+s(ave_size_PFO_1500_zeros, k=4)+s(PFO_F_Relative_Temp_500,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num, k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11,],family=binomial, select=T)

summary(fullmodel_mel_binaryEEE_1)



png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/abund_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)
plot_smooth_grey(fullmodel_mel_binaryEEE, view="ln_Mean_Abundance_perTrapnight",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=T, xlab=substitute(paste("ln(mean ",italic("Cs. melanura "),"abundance)", sep=" ")), ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,0))

# axis(2, at=c(-13.8, -11.5,-9.2,-6.9, -4.6, -2.3,0),labels = c(.000001,.00001,.0001,.001,.01,.1, 1),col="black",lwd=2, line=0)



dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/PFO_ever_1000_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PFO_Ever.x_1000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=T, xlab="Proportional area of evergreen-forested\nwetland (1000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,0))

dev.off()


png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/PFO_decid_5000_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PFO_Decid.x_5000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=T, xlab="Proportional area of deciduous-forested\nwetland (5000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-10,-4), xlim=c(0,.1))


dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/shrub_5000_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PSS_5000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=T, xlab="Proportional area of scrub/shurb\nwetland (5000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,0), xlim=c(0,.03))

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/emerg_2000_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PEM_2000.y",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=T, xlab="Proportional area of emergent\nwetland (2000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-25,-5))

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/ave_size_1500_zeros_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)


plot_smooth_grey(fullmodel_mel_binaryEEE, view="ave_size_PFO_1500_zeros",col='red',lwd=3, rug=F, main="", hide.label = T,h0=NA ,axes=T, xlab="Average forested wetland size\n(m2; 1500m Buffer)", ylab="Log odds of EEEV positive pool", ylim=c(-10,-5))

dev.off()


png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/impv_200_infect.png", width=2400, height=2400, res=300)

par(cex=1.8)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="impvMN_200",col='red',lwd=3, rug=F, main="", hide.label = T,h0=NA ,axes=T, xlab="Mean land-cover imperviousness (200m Buffer)", ylab="Log odds of EEEV positive pool", ylim=c(-14, -4))

dev.off()



# plot_smooth_grey(fullmodel_mel_binaryEEE, view="month_num",col='red',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Proportional area of Evergreen-forested wetland (500m Buffer)", ylab="Mean Cs. melanura Abundance per Night")


png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/abundance_relativesemi500_10.png", width=2400, height=2400, res=300)

par(cex=1.8)
fvisgam(fullmodel_mel_binaryEEE, view=c("PFO_F_Relative_Temp_500","prev_phdi_10"), color="topo",  hide.label = T, add.color.legend = T,too.far=.3,main="",ylab="PHDI (10 month lag)" ,xlab="Relative area semi-permanent forested\nwetland (500m buffer)")

dev.off()




##creating a panel of plots###

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/infection_panel_plot_3.png", width=4800, height=4800, res=300)

#layout(matrix(c(1,2,3,4,5,6,7,7,7), 3, 3, byrow = TRUE))
#par(cex=1.5)
par(mfrow=c(3,3),cex=1.5)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PFO_Decid.x_5000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of deciduous-forested\nwetland (5000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-10,-4), xlim=c(0,.1))
#something weird is happening and number isn't printing so have to add manually
axis(2, at=c(-10,-9,-8,-7,-6,-5,-4),label=c(-10,-9,-8,-7,-6,-5,-4),col="black",lwd=1, line=0, las=1)
axis(1)

mtext("(a)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 3%\np=0.005", side = 3,adj = 0.9, line = -1, cex=1.3)


plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PFO_Ever.x_1000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of evergreen-forested\nwetland (1000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,-2))

#something weird is happening and number isn't printing so have to add manually
axis(2, at=c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3,-2),label=c(-12,NA,-10,NA,-8,NA,-6,NA,-4, NA,-2),col="black",lwd=1, line=0, las=1 )
axis(1)

mtext("(b)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 0%\np=0.26", side = 3,adj = 0.9, line = -1, cex=1.3)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PEM_2000.y",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of emergent\nwetland (2000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-25,-5))

axis(2, at=c(-25,-20,-15,-10,-5),label=c(-25,-20,-15,-10,-5),col="black",lwd=1, line=0, las=1 )
axis(1)


mtext("(c)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 1%\np=0.01", side = 3,adj = 0.9, line = -1, cex=1.3)



plot_smooth_grey(fullmodel_mel_binaryEEE, view="prop_PSS_5000",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab="Proportional area of scrub/shurb\nwetland (5000m Buffer)", ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,-2), xlim=c(0,.03))

#something weird is happening and number isn't printing so have to add manually
axis(2, at=c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3,-2),label=c(-12,NA,-10,NA,-8,NA,-6,NA,-4, NA,-2),col="black",lwd=1, line=0, las=1 )
axis(1)


mtext("(d)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 0%\np=0.74", side = 3,adj = 0.9, line = -1, cex=1.3)


plot_smooth_grey(fullmodel_mel_binaryEEE, view="ave_size_PFO_1500_zeros",col='red',lwd=3, rug=F, main="", hide.label = T,h0=NA ,axes=F, xlab="Average forested wetland size\n(km2; 1500m Buffer)", ylab="Log odds of EEEV positive pool", ylim=c(-10,-5))


axis(2, at=c(-10,-9,-8,-7,-6,-5, -4),label=c(-10,-9,-8,-7,-6,-5,-4),col="black",lwd=1, line=0, las=1 )
#axis(1)
axis(1, at=(c(0,20,40,60,80,100)*1000),label=c(0,20,40,60,80,100) )


mtext("(e)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 1%\np=0.078", side = 3,adj = 0.9, line = -1, cex=1.3)

plot_smooth_grey(fullmodel_mel_binaryEEE, view="impvMN_200",col='red',lwd=3, rug=F, main="", hide.label = T,h0=NA ,axes=F, xlab="Mean land-cover imperviousness\n(200m Buffer)", ylab="Log odds of EEEV positive pool", ylim=c(-14, -4))

axis(2, at=c(-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4),label=c(-14,NA,-12,NA,-10,NA,-8,NA,-6,NA,-4),col="black",lwd=1, line=0, las=1 )
axis(1)


mtext("(f)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 1%\np=0.1", side = 3,adj = 0.9, line = -1, cex=1.3)



#par(pin=par()$pin, cex=1.5)
# plot_smooth_grey(fullmodel_mel_binaryEEE, view="ln_Mean_Abundance_perTrapnight",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab=substitute(paste("ln(mean ",italic("Cs. melanura "),"abundance + 1)", sep=" ")), ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,0))

plot_smooth_grey(fullmodel_mel_binaryEEE, view="ln_Mean_Abundance_perTrapnight",col='red',lwd=3, rug=FALSE, main="", hide.label = T, h0=NA ,axes=F, xlab=substitute(paste("Mean ",italic("Cs. melanura "),"capture per night", sep=" ")), ylab="Log odds of EEEV positive pool", shade = T, ylim=c(-12,0))

axis(2, at=c(-12,-11,-10,-9,-8,-7,-6,-5,-4, -3,-2,-1,0),label=c(-12,NA,-10,NA,-8,NA,-6,NA,-4, NA,-2,NA,0),col="black",lwd=1, line=0, las=1 )
#axis(1)

axis(1, at=log(c(0,1,5,15,50,150,500, 1000)+1), label=c(0,1,5,15,50,150,500,NA))



mtext("(g)", side = 3, adj = 0.05, line=0, las=1, cex=1.7)
mtext("%dev.= 15%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)

dev.off()




png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_infect/wetland_phdi_interaction_plots_infect.png", width=2400, height=4800, res=300)
par(mar=c(5,4,1.5,2.3),mfrow=c(2,1),cex=2.2)


 fvisgam(fullmodel_mel_binaryEEE, view=c("PFO_F_Relative_Temp_500","prev_phdi_10"), color="topo",  hide.label = T, add.color.legend = T,too.far=.4,main="",ylab="PHDI (10 month lag)" ,xlab="Relative area of semi-permanent forested\nwetland (500m buffer)")
 mtext("(a)", side = 3, adj = 0, line = .5, cex=2)
 
 
 fvisgam(fullmodel_mel_binaryEEE_1, view=c("PFO_F_Relative_Temp_500","prev_phdi_0"), color="topo",  hide.label = T, add.color.legend = T,too.far=.4,main="",ylab="PHDI (0 month lag)" ,xlab="Relative area of semi-permanent forested\nwetland (500m buffer)")
 mtext("(b)", side = 3, adj = 0, line = .5, cex=2)
 
 
dev.off()



