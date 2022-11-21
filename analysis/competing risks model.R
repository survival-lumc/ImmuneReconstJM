library(targets)
library(dplyr)
library(survival)
library(prodlim)

tar_load(NMA_preDLI_datasets)
dat_wide <- NMA_preDLI_datasets$wide
table(dat_wide$endpoint_specify7)

fit_compEvents_ITT <- prodlim(
  Hist(endpoint7, endpoint7_s, cens.code = "cens") ~ hirisk,
  data = NMA_preDLI_datasets$wide%>%
    #filter(TCD%in%c("NMA RD: ALT","NMA UD: ALT + 1mg ATG","NMA UD: ALT + 2mg ATG"))%>%
    mutate(endpoint7_s=ifelse(endpoint7_s=="7 days after cellular intervention",
                              "censored",endpoint7_s))
)

summary(fit_compEvents_ITT, times=c(0,3,6))

#png("compEvents_ITT_all.png", width=18, height = 9, units="cm", pointsize = 11, bg = "white", res=300)
par(mar=c(6, 4, 2, 0) + 0.1, lend=1, cex=0.5,ljoin=1, mfrow=c(1,3))
plot(fit_compEvents_ITT, cause="non-relapse failure: GvHD", atrisk.at=c(0:6), xlim=c(0,6), col=c("blue","red"), marktime = T,
     atrisk.line=c(3.5,4.5), atrisk.cex=0.5,
     axes=T,percent=F, ylab="cumulative incidence", xlab="months since alloSCT\n ",
     axis2.las=1, confint = T, lwd=2, automar=F, ylim=c(0,1.01), axis2.at=seq(0,1,0.1), background.border="transparent",
     background = F, #background.horizontal=(0:4)*0.25,
     plot.main="clinically significant GvHD", axis1.pos=0, axis2.pos=0, legend=F,
     atrisk.title="No. at risk  ", atrisk.labels=c("non-high risk  ", "high risk  "))
lines(x=c(0,6),y=c(1,1))
lines(x=c(6,6),y=c(0,1))
legend("topright", bty="n", legend=c("non-high risk", "high risk"),
       col=c("blue","red"), title="", lty=1, lwd=2, inset=c(0.05,0))
plot(fit_compEvents_ITT, cause="relapse", atrisk.at=c(0:6), xlim=c(0,6), col=c("blue","red"), marktime = T,
     atrisk.line=c(3.5,4.5), atrisk.cex=0.5,
     axes=T,percent=F, ylab="", xlab="months since alloSCT\n ",
     axis2.las=1, confint = T, lwd=2, automar=F, ylim=c(0,1.01), axis2.at=seq(0,1,0.1), background.border="transparent",
     background = F, #background.horizontal=(0:4)*0.25,
     plot.main="relapse", axis1.pos=0, axis2.pos=0, legend=F,
     atrisk.title="", atrisk.labels=c("", ""))
lines(x=c(0,6),y=c(1,1))
lines(x=c(6,6),y=c(0,1))
plot(fit_compEvents_ITT, cause="non-relapse failure: other", atrisk.at=c(0:6), xlim=c(0,6), col=c("blue","red"), marktime = T,
     atrisk.line=c(3.5,4.5), atrisk.cex=0.5,
     axes=T,percent=F, ylab="", xlab="months since alloSCT\n ",
     axis2.las=1, confint = T, lwd=2, automar=F, ylim=c(0,1.01), axis2.at=seq(0,1,0.1), background.border="transparent",
     background = F, #background.horizontal=(0:4)*0.25,
     plot.main="other failure", axis1.pos=0, axis2.pos=0, legend=F,
     atrisk.title="", atrisk.labels=c("", ""))
lines(x=c(0,6),y=c(1,1))
lines(x=c(6,6),y=c(0,1))
#dev.off()
