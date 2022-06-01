### Summarize results from the analysis of real datasets
### Analyses from Closed Population (CP), maybe Habitat Management Unit (HMU), and others

library(jagsUI);library(ggplot2);library(RColorBrewer);library(tidyverse);library(ggpubr);library(HDInterval)

#### VISUAL SURVEYS ####
## Dataset analyses below not archived on ScienceBase as used combined analysis for rest of manuscript
mdlsVIS <- c("Real Data Analysis/Visual surveys/Results/Unused/NWFNVIS2_SCRpstar.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISHL1_SCRpstar.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISHL2_SCRpstar.Rdata",
          'Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPREBT2_SCRpstar.Rdata',
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTBT2_SCRpstar.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTKB1_SCRpstar.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTKB2_SCRpstar.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISPOSTKB3_SCRpstar.Rdata",
          "Real Data Analysis/Visual surveys/Results/Unused/NWFNVISTRAPVIS_SCRpstar.Rdata")

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
  filter(Type == "sigma" | Type == "N") %>%
  summarise(Mean = mean(Mean))


pVISapp1 <- ggplot(data = subset(datVIS, Type == c("sigma") | Type == c("N")), aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = 2, labeller = label_parsed) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=subset(datVIS, Type == c("sigma") | Type == c("N")), aes(ymin=HDPI2.5, ymax=HDPI97.5)) +
  theme(legend.position = "none", axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 12)) + 
  ggtitle("Visual Surveys") +
  geom_hline(data = dV2Mean, aes(yintercept = Mean), alpha=0.4, colour="#1d92f1")
  


#### TRAPPING ####
## Dataset analyses below not archived on ScienceBase as used combined analysis for rest of manuscript
mdlsTRAP <- c("Real Data Analysis/Trapping/Results/Unused/NWFNTRAP1_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNTRAP2LINVIS_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNTRAP3_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNTRAP4LCM_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNVISTRAPTRAP_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPOSTBT2TRAP_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPOSTKBTRAP1_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPOSTKBTRAP2_SCRpstar.Rdata",
          "Real Data Analysis/Trapping/Results/Unused/NWFNPREBT1TRAP_SCRpstar.Rdata")

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
  filter(Type == "sigma" | Type == "N") %>%
  summarise(Mean = mean(Mean))

pTRAPapp1 <- ggplot(data = subset(datTRAP, Type == c("sigma") | Type == c("N")), aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = 2, labeller = label_parsed) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=subset(datTRAP, Type == c("sigma") | Type == c("N")), aes(ymin=HDPI2.5, ymax=HDPI97.5)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 14), strip.text.x = element_text(size = 12)) + 
  ggtitle("Trapping Surveys") +
  xlab("Project") +
  geom_hline(data = dT2Mean, aes(yintercept = Mean), alpha=0.4, colour="#1d92f1")


### Create combined plot for appendix

png(file="Real Data Analysis/Figures/VISTRAPsigma&Ncomparison.png",width=8,height=8.5,units="in",res=600)
ggarrange(pVISapp1, pTRAPapp1, nrow = 2, ncol = 1, labels = "AUTO")
dev.off()

