### Supplementary figure of sigma, p0, and abundance across projects
### The purpose of this code is to summarize results from the separate analyses of real datasets
### All analyses are from the Closed Population (CP)

library(jagsUI);library(ggplot2);library(RColorBrewer);library(tidyverse);library(ggpubr);library(HDInterval)

#### VISUAL SURVEYS ####
## Dataset analyses below are NOT archived on ScienceBase as we used combined analysis for rest of manuscript
## However, we include the code for how we created this plot
mdlsVIS <- c("Real Data Analysis/Visual surveys/Results/Unused/NWFNVIS2_SCRpstarvisCATsizeCAT.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISHL1_SCRpstarvisCATsizeCAT.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISHL2_SCRpstarvisCATsizeCAT.Rdata",
          'Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPREBT2_SCRpstarvisCATsizeCAT.Rdata',
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTBT2_SCRpstarvisCATsizeCAT.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTKB1_SCRpstarvisCATsizeCAT.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTKB2_SCRpstarvisCATsizeCAT.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTKB3_SCRpstarvisCATsizeCAT.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISTRAPVIS_SCRpstarvisCATsizeCAT.Rdata")

yearsVIS <- c(2006,2007,2007,2009,2010,2011,2012,2012,2015)
areas <- c(rep(50000,9))

## read in list of model results
datVIS <- as.data.frame(matrix(NA, nrow = (14*length(mdlsVIS)), ncol = 9))
key <- 1
for(i in 1:length(mdlsVIS)){
  load(mdlsVIS[i])
  datVIS[key,1] <- out$mean$sigma
  datVIS[key,2] <- out$q2.5$sigma
  datVIS[key,3] <- out$q97.5$sigma
  datVIS[key,4] <- c("sigma")
  datVIS[key,5] <- c("sigma")
  datVIS[key,6] <- paste("Model",i,sep="")
  datVIS[key,7] <- i
  datVIS[key,8] <- yearsVIS[i]
  datVIS[key,9] <- out$mean$N/areas[i]
  datVIS[key,10] <- hdi(out$sims.list$sigma)[[1]]
  datVIS[key,11] <- hdi(out$sims.list$sigma)[[2]]
  key <- key + 1
  for(j in 1:length(out$mean$pstar)){
    datVIS[key,1] <- out$mean$pstar[j]
    datVIS[key,2] <- out$q2.5$pstar[j]
    datVIS[key,3] <- out$q97.5$pstar[j]
    datVIS[key,4] <- c(paste("pstar",j,sep=""))
    datVIS[key,5] <- c("pstar")
    datVIS[key,6] <- paste("Model",i,sep="")
    datVIS[key,7] <- paste(i,j,sep=".")
    datVIS[key,8] <- yearsVIS[i]
    datVIS[key,9] <- out$mean$N/areas[i]
    datVIS[key,10] <- hdi(out$sims.list$pstar)[[1]]
    datVIS[key,11] <- hdi(out$sims.list$pstar)[[2]]
    key <- key + 1
  }
  for(j in 1:length(out$mean$p0)){
    datVIS[key,1] <- out$mean$p0[j]
    datVIS[key,2] <- out$q2.5$p0[j]
    datVIS[key,3] <- out$q97.5$p0[j]
    datVIS[key,4] <- c(paste("p0",j,sep=""))
    datVIS[key,5] <- c("p0")
    datVIS[key,6] <- paste("Model",i,sep="")
    datVIS[key,7] <- paste(i,j,sep=".")
    datVIS[key,8] <- yearsVIS[i]
    datVIS[key,9] <- out$mean$N/areas[i]
    datVIS[key,10] <- hdi(out$sims.list$p0)[[1]]
    datVIS[key,11] <- hdi(out$sims.list$p0)[[2]]
    key <- key + 1
  }
  datVIS[key,1] <- out$mean$N
  datVIS[key,2] <- out$q2.5$N
  datVIS[key,3] <- out$q97.5$N
  datVIS[key,4] <- c("N")
  datVIS[key,5] <- c("N")
  datVIS[key,6] <- paste("Model",i,sep="")
  datVIS[key,7] <- i
  datVIS[key,8] <- yearsVIS[i]
  datVIS[key,9] <- out$mean$N/areas[i]
  datVIS[key,10] <- hdi(out$sims.list$N)[[1]]
  datVIS[key,11] <- hdi(out$sims.list$N)[[2]]
  key <- key + 1
  for(j in 1:length(out$mean$Ngroup)){
    datVIS[key,1] <- out$mean$Ngroup[j]
    datVIS[key,2] <- out$q2.5$Ngroup[j]
    datVIS[key,3] <- out$q97.5$Ngroup[j]
    datVIS[key,4] <- c(paste("Ngroup",j,sep=""))
    datVIS[key,5] <- c("Ngroup")
    datVIS[key,6] <- paste("Model",i,sep="")
    datVIS[key,7] <- paste(i,j,sep=".")
    datVIS[key,8] <- yearsVIS[i]
    datVIS[key,9] <- out$mean$N/areas[i]
    datVIS[key,10] <- hdi(out$sims.list$Ngroup)[[1]]
    datVIS[key,11] <- hdi(out$sims.list$Ngroup)[[2]]
    key <- key + 1
  }
}
colnames(datVIS) <- c("Mean","Q2.5","Q97.5","Parameter","Type","Model","ParModel","Year","Density","HDPI2.5","HDPI97.5")
## Remove NA rows when a project dataset was missing a snake size category
datVIS <- na.omit(datVIS)

cols <- brewer.pal(length(mdlsVIS), "Greens") 


### Code for appendix - compare N and sigma across projects

dV2Mean <- datVIS %>%
  group_by(Type) %>%
  filter(Type == "sigma" | Type == "N"  | Type == "p0") %>%
  summarise(Mean = mean(Mean))

### Figure option 1 - separate N, sigma, and p0 plots

pVISapp1 <- ggplot(data = subset(datVIS, Type == c("sigma") | Type == c("N") | Type == c("p0")), aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = 3, labeller = label_parsed) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=subset(datVIS, Type == c("sigma") | Type == c("N") | Type == c("p0")), aes(ymin=HDPI2.5, ymax=HDPI97.5)) +
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 12)) + 
  ggtitle("Visual Surveys") +
  geom_hline(data = dV2Mean, aes(yintercept = Mean), alpha=0.4, colour="#1d92f1")

#combine data into different format for different kind of plot (add NA row to Model 6 where no data were available on size category 4 snakes - so no p04)
altdatVIS <- cbind(subset(datVIS, Parameter == "N")[,c(1,6,8,10:11)],subset(datVIS, Parameter == "sigma")[,c(1,10:11)],subset(datVIS, Parameter == "p01")[,c(1,10:11)],subset(datVIS, Parameter == "p02")[,c(1,10:11)],subset(datVIS, Parameter == "p03")[,c(1,10:11)],add_row(subset(datVIS, Parameter == "p04")[,c(1,10:11)], Mean = NA, HDPI2.5 = NA, HDPI97.5 = NA, .before = 6))
colnames(altdatVIS) <- c("NMean","Model","Year","NHDPI2.5","NHDPI97.5","sigMean","sigHDPI2.5","sigHDPI97.5","p01Mean","p01HDPI2.5","p01HDPI97.5","p02Mean","p02HDPI2.5","p02HDPI97.5","p03Mean","p03HDPI2.5","p03HDPI97.5","p04Mean","p04HDPI2.5","p04HDPI97.5")

pVISapp1alt <- ggplot(data = altdatVIS, aes(x=NMean, y=sigMean)) + 
  geom_point(aes(shape=Model), size=2) +
  scale_shape_manual(values=c(0,1,2,5,6,15,16,17,18)) +
  geom_linerange(data=altdatVIS, aes(xmin=NHDPI2.5, xmax=NHDPI97.5)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 12)) +
  ggtitle("Visual Surveys")

pVISapp1altP0 <- ggplot() + 
  geom_point(data = altdatVIS, aes(x=NMean, y=p01Mean), size=2, shape=21, fill=cols[1]) +
  geom_point(data = altdatVIS, aes(x=NMean, y=p02Mean), size=2, shape=22, fill=cols[3]) +
  geom_point(data = altdatVIS, aes(x=NMean, y=p03Mean), size=2, shape=23, fill=cols[5]) +
  geom_point(data = altdatVIS, aes(x=NMean, y=p04Mean), size=2, shape=24, fill=cols[7]) +
  ggtitle("Visual Surveys") +
  ylab("Mean p0")
#ignore warning about missing value - there should be as there is an NA row in the dataset. One dataset lacked info on size category 4 and did not have a p0 for this category.
  


#### TRAPPING ####
## Dataset analyses below not archived on ScienceBase as used combined analysis for rest of manuscript
mdlsTRAP <- c("Real Data Analysis/Trapping/Results/Unused/NWFNTRAP1_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNTRAP2LINVIS_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNTRAP3_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNTRAP4LCM_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNVISTRAPTRAP_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPOSTBT2TRAP_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPOSTKBTRAP1_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPOSTKBTRAP2_SCRpstartrapCATsizeCAT.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPREBT1TRAP_SCRpstartrapCATsizeCAT.Rdata")

yearsTRAP <- c(2004,2005,2007,2007,2007,2009,2012,2013,2008)
areas <- c(rep(50000,9))

## read in list of model results
datTRAP <- as.data.frame(matrix(NA, nrow = (14*length(mdlsTRAP)), ncol = 9))
key <- 1
for(i in 1:length(mdlsTRAP)){
  load(mdlsTRAP[i])
  datTRAP[key,1] <- out$mean$sigma
  datTRAP[key,2] <- out$q2.5$sigma
  datTRAP[key,3] <- out$q97.5$sigma
  datTRAP[key,4] <- c("sigma")
  datTRAP[key,5] <- c("sigma")
  datTRAP[key,6] <- paste("Model",i,sep="")
  datTRAP[key,7] <- i
  datTRAP[key,8] <- yearsTRAP[i]
  datTRAP[key,9] <- out$mean$N/areas[i]
  datTRAP[key,10] <- hdi(out$sims.list$sigma)[[1]]
  datTRAP[key,11] <- hdi(out$sims.list$sigma)[[2]]
  key <- key + 1
  for(j in 1:length(out$mean$pstar)){
    datTRAP[key,1] <- out$mean$pstar[j]
    datTRAP[key,2] <- out$q2.5$pstar[j]
    datTRAP[key,3] <- out$q97.5$pstar[j]
    datTRAP[key,4] <- c(paste("pstar",j,sep=""))
    datTRAP[key,5] <- c("pstar")
    datTRAP[key,6] <- paste("Model",i,sep="")
    datTRAP[key,7] <- paste(i,j,sep=".")
    datTRAP[key,8] <- yearsTRAP[i]
    datTRAP[key,9] <- out$mean$N/areas[i]
    datTRAP[key,10] <- hdi(out$sims.list$pstar)[[1]]
    datTRAP[key,11] <- hdi(out$sims.list$pstar)[[2]]
    key <- key + 1
  }
  for(j in 1:length(out$mean$p0)){
    datTRAP[key,1] <- out$mean$p0[j]
    datTRAP[key,2] <- out$q2.5$p0[j]
    datTRAP[key,3] <- out$q97.5$p0[j]
    datTRAP[key,4] <- c(paste("p0",j,sep=""))
    datTRAP[key,5] <- c("p0")
    datTRAP[key,6] <- paste("Model",i,sep="")
    datTRAP[key,7] <- paste(i,j,sep=".")
    datTRAP[key,8] <- yearsTRAP[i]
    datTRAP[key,9] <- out$mean$N/areas[i]
    datTRAP[key,10] <- hdi(out$sims.list$p0)[[1]]
    datTRAP[key,11] <- hdi(out$sims.list$p0)[[2]]
    key <- key + 1
  }
  datTRAP[key,1] <- out$mean$N
  datTRAP[key,2] <- out$q2.5$N
  datTRAP[key,3] <- out$q97.5$N
  datTRAP[key,4] <- c("N")
  datTRAP[key,5] <- c("N")
  datTRAP[key,6] <- paste("Model",i,sep="")
  datTRAP[key,7] <- i
  datTRAP[key,8] <- yearsTRAP[i]
  datTRAP[key,9] <- out$mean$N/areas[i]
  datTRAP[key,10] <- hdi(out$sims.list$N)[[1]]
  datTRAP[key,11] <- hdi(out$sims.list$N)[[2]]
  key <- key + 1
  for(j in 1:length(out$mean$Ngroup)){
    datTRAP[key,1] <- out$mean$Ngroup[j]
    datTRAP[key,2] <- out$q2.5$Ngroup[j]
    datTRAP[key,3] <- out$q97.5$Ngroup[j]
    datTRAP[key,4] <- c(paste("Ngroup",j,sep=""))
    datTRAP[key,5] <- c("Ngroup")
    datTRAP[key,6] <- paste("Model",i,sep="")
    datTRAP[key,7] <- paste(i,j,sep=".")
    datTRAP[key,8] <- yearsTRAP[i]
    datTRAP[key,9] <- out$mean$N/areas[i]
    datTRAP[key,10] <- hdi(out$sims.list$Ngroup)[[1]]
    datTRAP[key,11] <- hdi(out$sims.list$Ngroup)[[2]]
    key <- key + 1
  }
}
colnames(datTRAP) <- c("Mean","Q2.5","Q97.5","Parameter","Type","Model","ParModel","Year","Density","HDPI2.5","HDPI97.5")
## Remove NA rows when a project dataset was missing a snake size category
datTRAP <- na.omit(datTRAP)

cols <- brewer.pal(length(mdlsTRAP), "Greens") 


### Code for appendix - compare N and sigma across projects

dT2Mean <- datTRAP %>%
  group_by(Type) %>%
  filter(Type == "sigma" | Type == "N" | Type == "p0") %>%
  summarise(Mean = mean(Mean))

pTRAPapp1 <- ggplot(data = subset(datTRAP, Type == c("sigma") | Type == c("N") | Type == c("p0")), aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = 3, labeller = label_parsed) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=subset(datTRAP, Type == c("sigma") | Type == c("N") | Type == c("p0")), aes(ymin=HDPI2.5, ymax=HDPI97.5)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 12)) + 
  ggtitle("Trapping Surveys") +
  xlab("Project") +
  geom_hline(data = dT2Mean, aes(yintercept = Mean), alpha=0.4, colour="#1d92f1")


#combine data into different format for different kind of plot
altdatTRAP <- cbind(subset(datTRAP, Parameter == "N")[,c(1,6,8,10:11)],subset(datTRAP, Parameter == "sigma")[,c(1,10:11)],subset(datTRAP, Parameter == "p01")[,c(1,10:11)],subset(datTRAP, Parameter == "p02")[,c(1,10:11)],subset(datTRAP, Parameter == "p03")[,c(1,10:11)],add_row(subset(datTRAP, Parameter == "p04")[,c(1,10:11)], Mean = NA, HDPI2.5 = NA, HDPI97.5 = NA, .before = 8))
colnames(altdatTRAP) <- c("NMean","Model","Year","NHDPI2.5","NHDPI97.5","sigMean","sigHDPI2.5","sigHDPI97.5","p01Mean","p01HDPI2.5","p01HDPI97.5","p02Mean","p02HDPI2.5","p02HDPI97.5","p03Mean","p03HDPI2.5","p03HDPI97.5","p04Mean","p04HDPI2.5","p04HDPI97.5")

pTRAPapp1alt <- ggplot(data = altdatTRAP, aes(x=NMean, y=sigMean)) + 
  geom_point(aes(shape=Model), size=2) +
  scale_shape_manual(values=c(0,1,2,5,6,15,16,17,18)) +
  geom_linerange(data=altdatTRAP, aes(xmin=NHDPI2.5, xmax=NHDPI97.5)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 12)) +
  ggtitle("Trapping Surveys")

### Create combined plot for appendix

png(file="Real Data Analysis/Figures/VISTRAPsigma&N&p0comparison.png",width=8,height=8.5,units="in",res=600)
ggarrange(pVISapp1, pTRAPapp1, nrow = 2, ncol = 1, labels = "AUTO")
dev.off()

png(file="Real Data Analysis/Figures/VISTRAPsigmaVSNcomparison.png",width=8,height=8.5,units="in",res=600)
ggarrange(pVISapp1alt, pTRAPapp1alt, nrow = 2, ncol = 1, labels = "AUTO")
dev.off()

