
load("~/Documents/Grad School/Edited R Packages/plot_error_grey.Rdata")
load("~/Documents/Grad School/Edited R Packages/plot_smooth_grey.Rdata")


##maybe only use scales of figure out a way to remove outliers

fullmodel_melanura_abund<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_100.y, k=3)+s(prop_PFO_Ever.x_1000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_2000_zeros, k=4)+s(prop_PFO_F_3000,prev_phdi_1,k=4)+ti(prop_PFO_3000,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$prop_PFO_Ever.x_500<.3&fulldata$PFO_mean_stream_cnt_1500_zeros<5, ], select=T)


summary(fullmodel_melanura_abund)
gam.check(fullmodel_melanura_abund)

plot(fullmodel_melanura_abund)


fullmodel_melanura_abund_1<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_100.y, k=3)+s(prop_PFO_Ever.x_1000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_2000_zeros, k=4)+s(prop_PFO_F_3000,prev_phdi_1,k=4)+ti(prop_PFO_2000,prev_phdi_8,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$prop_PFO_Ever.x_500<.3&fulldata$PFO_mean_stream_cnt_1500_zeros<5&fulldata$prop_PFO_2000<.15, ], select=T)


summary(fullmodel_melanura_abund_1)

plot(fullmodel_melanura_abund_1)

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/PEM_100_abund.png", width=2400, height=2400, res=300)

par(cex=1.9)
plot_smooth_grey(fullmodel_melanura_abund, view="prop_PEM_100.y",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Proportional area of\n emergent wetland (100m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim = c(0,2.5), shade = T)

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/PFO_ever_1000_abund.png", width=2400, height=2400, res=300)

par(cex=1.9)
plot_smooth_grey(fullmodel_melanura_abund, view="prop_PFO_Ever.x_1000",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Proportional area of\n evergreen-forested wetland (1000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(1,6))

dev.off()


png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/PFO_decid_3000_abund.png", width=2400, height=2400, res=300)

par(cex=1.9)
plot_smooth_grey(fullmodel_melanura_abund, view="prop_PFO_Decid.x_3000",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Proportional area of\n deciduous-forested wetland (3000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(0,15))

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/shrub_3000_abund.png", width=2400, height=2400, res=300)

par(cex=1.9)
plot_smooth_grey(fullmodel_melanura_abund, view="prop_PSS_3000",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Proportional area of\n scrub/shurb wetland (3000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(0,2.5), xlim=c(0,.05))

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/stream_connections_2000_zeros_abund.png", width=2400, height=2400, res=300)

par(cex=1.9)
plot_smooth_grey(fullmodel_melanura_abund, view="PFO_mean_stream_cnt_2000_zeros",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Mean # stream connections (2000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(.5,2), xlim=c(.5,3))

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/impv_200_abund.png", width=2400, height=2400, res=300)

par(cex=1.9)
plot_smooth_grey(fullmodel_melanura_abund, view="impvMN_200",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Mean land-cover imperviousness (200m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")),ylim=c(0,3), xlim=c(0,50))

dev.off()


# plot_smooth_grey(fullmodel_melanura_abund, view="month_num",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=T, xlab="Proportional area of Evergreen-forested wetland (500m Buffer)", ylab="Mean Cs. melanura Abundance per Night")


png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/abundance_semiforest3000_1.png", width=2400, height=2400, res=300)

par(cex=1.9)
fvisgam(fullmodel_melanura_abund, view=c("prop_PFO_F_3000", "prev_phdi_1"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (1 month lag)" ,xlab="Proportional area semi-permanent\n forested wetland (3000m buffer)")

dev.off()

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/abundance_forest3000_0.png", width=2400, height=2400, res=300)

par(cex=1.9)
fvisgam(fullmodel_melanura_abund, view=c("prop_PFO_3000","prev_phdi_0"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (0 month lag)" ,xlab="Proportional area of forested\n wetlands (3000m buffer)")

dev.off()


##creating a panel of plots###

png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/abundance_panel_plot.png", width=4800, height=3600, res=300)

par(mfrow=c(2,3),cex=1.5)

plot_smooth_grey(fullmodel_melanura_abund, view="prop_PFO_Decid.x_3000",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of deciduous\nforested wetland (3000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(0,15))
axis(2,las=1)
axis(1)
mtext("(a)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 13%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)




plot_smooth_grey(fullmodel_melanura_abund, view="prop_PFO_Ever.x_1000",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of evergreen\nforested wetland (1000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(1,6))
axis(2,las=1)
axis(1)
mtext("(b)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 7%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)


plot_smooth_grey(fullmodel_melanura_abund, view="prop_PEM_100.y",col='blue',lwd=3, rug=FALSE, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of emergent\nwetland (100m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim = c(0,2.5), shade = T)
axis(2,las=1)
axis(1)
mtext("(c)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 1%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)



plot_smooth_grey(fullmodel_melanura_abund, view="prop_PSS_3000",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Proportional area of scrub/shurb\nwetland (3000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(0,2.5), xlim=c(0,.05))
axis(2,las=1)
axis(1)
mtext("(d)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 2%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)


plot_smooth_grey(fullmodel_melanura_abund, view="PFO_mean_stream_cnt_2000_zeros",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Mean # stream connections\n(2000m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")), ylim=c(.5,2), xlim=c(.5,3))
axis(2,las=1)
axis(1)
mtext("(e)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 1%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)


plot_smooth_grey(fullmodel_melanura_abund, view="impvMN_200",col='blue',lwd=3, rug=F, main="", hide.label = T, transform=expm1,h0=NA ,axes=F, xlab="Mean land-cover imperviousness\n(200m Buffer)", ylab=substitute(paste("Mean ",italic("Cs. melanura "),"abundance", sep=" ")),ylim=c(0,3), xlim=c(0,50))
axis(2,las=1)
axis(1)
mtext("(f)", side = 3, adj = 0.05, line = 0, cex=1.7)
mtext("%dev.= 11%\np<0.0001", side = 3,adj = 0.9, line = -1, cex=1.3)

dev.off()


png("~/Documents/Grad School/Dissertation Chapt. 2/Manuscript Prep/Manuscript_figures_3_10_17/Melanura_abund/wetland_phdi_interaction_plots_abund_new.png", width=2400, height=7000, res=300)
par(mar=c(5,4,2,2),mfrow=c(3,1),cex=2.2)

fvisgam(fullmodel_melanura_abund, view=c("prop_PFO_F_3000", "prev_phdi_1"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (1 month lag)" ,xlab="Proportional area semi-permanent\n forested wetland (3000m buffer)", xlim=c(0,.0012))
mtext("(a)", side = 3, adj = 0, line = .5, cex=2)

fvisgam(fullmodel_melanura_abund, view=c("prop_PFO_3000","prev_phdi_0"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (0 month lag)" ,xlab="Proportional area of forested\n wetlands (3000m buffer)")
mtext("(b)", side = 3, adj = 0, line = .5, cex=2)

fvisgam(fullmodel_melanura_abund_1, view=c("prop_PFO_2000","prev_phdi_8"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (8 month lag)" ,xlab="Proportional area of forested\n wetlands (2000m buffer)")
mtext("(c)", side = 3, adj = 0, line = .5, cex=2)


dev.off()



########
fullmodel_melanura_abundXXX<-gam(ln_Mean_Abundance_perTrapnight~Prop_Gravid_Trap+s(prop_PEM_100.y, k=3)+s(prop_PFO_Ever.x_1000, k=3)+s(prop_PFO_Decid.x_3000, k=4)+s(prop_PSS_3000,k=3)+s(PFO_mean_stream_cnt_2000_zeros, k=4)+s(prop_PFO_F_3000,prev_phdi_1,k=4)+ti(prop_PFO_3000,prev_phdi_0,k=4)+s(impvMN_200, k=4)+s(month_num,k=4),data=fulldata[month(fulldata$month_year)!=5&month(fulldata$month_year)!=11 &fulldata$prop_PFO_Ever.x_500<.3&fulldata$PFO_mean_stream_cnt_1500_zeros<5&fulldata$prop_PFO_F_3000<.002, ], select=T)

fvisgam(fullmodel_melanura_abund_1, view=c("prop_PFO_F_3000", "prev_phdi_1"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.3,main="",ylab="PHDI (1 month lag)" ,xlab="Proportional area semi-permanent\n forested wetland (3000m buffer)", xlim=c(0,.0015))
mtext("(a)", side = 3, adj = 0, line = .5, cex=2)

fvisgam(fullmodel_melanura_abund_1, view=c("prop_PFO_2000","prev_phdi_8"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (8 month lag)" ,xlab="Proportional area of forested\n wetlands (2000m buffer)",xlim=c(0,0.15))
mtext("(c)", side = 3, adj = 0, line = .5, cex=2)

fvisgam(fullmodel_melanura_abundXXX, view=c("prop_PFO_3000","prev_phdi_0"), color="topo", transform=expm1, hide.label = T, add.color.legend = T,too.far=.2,main="",ylab="PHDI (8 month lag)" ,xlab="Proportional area of forested\n wetlands (2000m buffer)")
mtext("(c)", side = 3, adj = 0, line = .5, cex=2)



