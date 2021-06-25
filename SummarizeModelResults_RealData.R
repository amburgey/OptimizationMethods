### Summarize results from the analysis of real datasets
### Analyses from Closed Population (CP), maybe Habitat Management Unit (HMU), and others

library(jagsUI);library(ggplot2);library(RColorBrewer);library(tidyverse)

#### VISUAL SURVEYS ####
## Only showing NWFN so far
mdlsVIS <- c("Visual surveys/Results/NWFNVIS2_SCRpstarvisCATsizeCAT.Rdata",
          "Visual surveys/Results/NWFNVISHL1_SCRpstarvisCATsizeCAT.Rdata",
          "Visual surveys/Results/NWFNVISHL2_SCRpstarvisCATsizeCAT.Rdata",
          'Visual surveys/Results/NWFNVISPREBT2_SCRpstarvisCATsizeCAT.Rdata',
          "Visual surveys/Results/NWFNVISPOSTBT2_SCRpstarvisCATsizeCAT.Rdata",
          "Visual surveys/Results/NWFNVISPOSTKB1_SCRpstarvisCATsizeCAT2sizes.Rdata",
          "Visual surveys/Results/NWFNVISPOSTKB2_SCRpstarvisCATsizeCAT.Rdata",
          "Visual surveys/Results/NWFNVISPOSTKB3_SCRpstarvisCATsizeCAT3sizes.Rdata",
          "Visual surveys/Results/NWFNVISTRAPVIS_SCRpstarvisCATsizeCAT.Rdata")

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
    key <- key + 1
  }
}
colnames(datVIS) <- c("Mean","Q2.5","Q97.5","Parameter","Type","Model","ParModel","Year","Density")
## Remove NA rows when a reduced number of size categories were used
datVIS <- na.omit(datVIS)

cols <- brewer.pal(length(mdlsVIS), "Greens") 

pVISall <- ggplot(dat = datVIS, aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = length(unique(datVIS$Type))) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=datVIS, aes(ymin=Q2.5, ymax=Q97.5)) +
  # scale_fill_manual(values=c("#00AFBB","#eda306","#dce035","#0e07a9","#84eca4","#be713d","#2aa012","#fe5ea4")) +
  theme(legend.position = "none")

dVMean <- datVIS %>%
  group_by(Type) %>%
  filter(Type == "sigma" | Type == "p0" | Type == "N") %>%
  summarise(Mean = mean(Mean))

pVISsub <- ggplot(data = subset(datVIS, Type == c("sigma") | Type == c("p0") | Type == c("N")), aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = length(unique(datVIS$Type))) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=subset(datVIS, Type == c("sigma") | Type == c("p0") | Type == c("N")), aes(ymin=Q2.5, ymax=Q97.5)) +
  # scale_fill_manual(values=c("#00AFBB","#eda306","#dce035","#0e07a9","#84eca4","#be713d","#2aa012","#fe5ea4")) +
  theme(legend.position = "none") +
  geom_hline(data = dVMean, aes(yintercept = Mean), alpha=0.4, colour="#1d92f1")
  


#### TRAPPING ####
## Only showing NWFN so far
mdlsTRAP <- c("Trapping/Results/NWFNTRAP1_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNTRAP2LINVIS_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNTRAP3_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNTRAP4LCM_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNVISTRAPTRAP_SCRpstartrapCATsizeCAT3500.Rdata",
          "Trapping/Results/NWFNPOSTBT2TRAP_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNPOSTKBTRAP1_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNPOSTKBTRAP2_SCRpstartrapCATsizeCAT.Rdata",
          "Trapping/Results/NWFNPREBT1TRAP_SCRpstartrapCATsizeCAT.Rdata")

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
    key <- key + 1
  }
}
colnames(datTRAP) <- c("Mean","Q2.5","Q97.5","Parameter","Type","Model","ParModel","Year","Density")
## Remove NA rows when a reduced number of size categories were used
datTRAP <- na.omit(datTRAP)

cols <- brewer.pal(length(mdlsTRAP), "Greens") 

pTRAPall <- ggplot(data = datTRAP, aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = length(unique(datTRAP$Type))) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=datTRAP, aes(ymin=Q2.5, ymax=Q97.5)) +
  # scale_fill_manual(values=c("#00AFBB","#eda306","#dce035","#0e07a9","#84eca4","#be713d","#2aa012","#fe5ea4")) +
  theme(legend.position = "none")

dTMean <- datTRAP %>%
  group_by(Type) %>%
  filter(Type == "sigma" | Type == "p0" | Type == "N") %>%
  summarise(Mean = mean(Mean))

pTRAPsub <- ggplot(data = subset(datTRAP, Type == c("sigma") | Type == c("p0") | Type == c("N")), aes(x=ParModel, y=Mean)) + 
  geom_point(aes(shape=Model, fill=Model), size=2) +
  facet_wrap(vars(Type), scales = "free", nrow = length(unique(datTRAP$Type))) +
  scale_shape_manual(values=c(21,17,23,19,21,17,23,19,21)) +
  scale_fill_manual(values=cols) +
  geom_linerange(data=subset(datTRAP, Type == c("sigma") | Type == c("p0") | Type == c("N")), aes(ymin=Q2.5, ymax=Q97.5)) +
  # scale_fill_manual(values=c("#00AFBB","#eda306","#dce035","#0e07a9","#84eca4","#be713d","#2aa012","#fe5ea4")) +
  theme(legend.position = "none") +
  geom_hline(data = dTMean, aes(yintercept = Mean), alpha=0.4, colour="#1d92f1")



#### COMPARE BETWEEN VIS AND TRAP and ACROSS TIME ####

### Sigma
sig <- rbind(subset(datVIS,Type == "sigma"),subset(datTRAP,Type == "sigma"))
sig$Method <- c(rep("VIS",nrow(sig)/2),rep("TRAP",nrow(sig)/2))

dSMean <- sig %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

lm_fit <- lm(Mean ~ Method, data=sig)

pVISTRAPsigma <- ggplot(data=sig, aes(x=Year,y=Mean, shape=Method, fill=Method)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23)) +
  scale_fill_manual(values=c("#00FA9A","#FFFACD")) +
  ylab("Mean - sigma") +
  geom_linerange(data=sig, aes(ymin=Q2.5, ymax=Q97.5)) +
  geom_smooth(method="lm")
  # geom_hline(data = dSMean, aes(yintercept = Mean), alpha=0.4, colour=c("#004d2f","#998a00"))

### Abundance
abund <- rbind(subset(datVIS,Type == "N"),subset(datTRAP,Type == "N"))
abund$Method <- c(rep("VIS",nrow(abund)/2),rep("TRAP",nrow(abund)/2))

dNMean <- abund %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

lm_fit <- lm(Mean ~ Method, data=abund)

pVISTRAPabund <- ggplot(data=abund, aes(x=Year,y=Mean, shape=Method, fill=Method)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23)) +
  scale_fill_manual(values=c("#00FA9A","#FFFACD")) +
  ylab("Mean - N") +
  geom_linerange(data=abund, aes(ymin=Q2.5, ymax=Q97.5))
  # geom_smooth(method="lm")

### Encounter probability
enc <- rbind(subset(datVIS,Type == "p0"),subset(datTRAP,Type == "p0"))
enc$Method <- c(rep("VIS",nrow(enc)/2),rep("TRAP",nrow(enc)/2))
enc$combo <- paste(enc$Type,enc$Parameter,enc$Method,sep=".")

dEMean <- enc %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

aggregate(enc[, c(1)], list(enc$Parameter), mean)

lm_fit <- lm(Mean ~ combo, data=enc)

pVISTRAPenc <- ggplot(data=enc, aes(x=Year,y=Mean, shape=combo, fill=combo)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23,21,23,21,23,21,23)) +
  scale_fill_manual(values=c("#00a868","#c7b300","#00c77b","#ffe605","#00FA9A","#fff38a","#a3ffdc","#FFFACD")) +
  # geom_linerange(data=enc, aes(ymin=Q2.5, ymax=Q97.5)) +
  ylab("Mean - p0") #+
  # geom_smooth(method="lm", formula= y~x, se=FALSE)

pVISTRAPenc2 <- ggplot(data=enc, aes(x=Year,y=Mean, shape=Method, fill=Method)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23)) +
  scale_fill_manual(values=c("#00FA9A","#FFFACD")) +
  ylab("Mean - p0") +
  geom_smooth(method="lm")

### Detection probability
det <- rbind(subset(datVIS,Type == "pstar"),subset(datTRAP,Type == "pstar"))
det$Method <- c(rep("VIS",nrow(det)/2),rep("TRAP",nrow(det)/2))
det$combo <- paste(det$Type,det$Parameter,det$Method,sep=".")

dDMean <- det %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

lm_fit <- lm(Mean ~ combo, data=det)

pVISTRAPdet <- ggplot(data=det, aes(x=Year,y=Mean, shape=combo, fill=combo)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23,21,23,21,23,21,23)) +
  scale_fill_manual(values=c("#00a868","#c7b300","#00c77b","#ffe605","#00FA9A","#fff38a","#a3ffdc","#FFFACD")) +
  # geom_linerange(data=det, aes(ymin=Q2.5, ymax=Q97.5)) +
  ylab("Mean - pstar") #+
# geom_smooth(method="lm", formula= y~x, se=FALSE)

### Sigma by Abundance
sig <- rbind(subset(datVIS,Type == "sigma"),subset(datTRAP,Type == "sigma"))
sig$Method <- c(rep("VIS",nrow(sig)/2),rep("TRAP",nrow(sig)/2))

dSMean <- sig %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

lm_fit <- lm(Mean ~ Method, data=sig)

pVISTRAPsigmaabund <- ggplot(data=sig, aes(x=Density,y=Mean, shape=Method, fill=Method)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23)) +
  scale_fill_manual(values=c("#00FA9A","#FFFACD")) +
  ylab("Mean - sigma") + xlab("Mean - Density") +
  geom_linerange(data=sig, aes(ymin=Q2.5, ymax=Q97.5)) +
  geom_smooth(method="lm")
# geom_hline(data = dSMean, aes(yintercept = Mean), alpha=0.4, colour=c("#004d2f","#998a00"))

### Encounter probability by Abundance
enc <- rbind(subset(datVIS,Type == "p0"),subset(datTRAP,Type == "p0"))
enc$Method <- c(rep("VIS",nrow(enc)/2),rep("TRAP",nrow(enc)/2))
enc$combo <- paste(enc$Type,enc$Parameter,enc$Method,sep=".")

dEMean <- enc %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

lm_fit <- lm(Mean ~ combo, data=enc)

pVISTRAPenc <- ggplot(data=enc, aes(x=Density,y=Mean, shape=combo, fill=combo)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23,21,23,21,23,21,23)) +
  scale_fill_manual(values=c("#00a868","#c7b300","#00c77b","#ffe605","#00FA9A","#fff38a","#a3ffdc","#FFFACD")) +
  # geom_linerange(data=enc, aes(ymin=Q2.5, ymax=Q97.5)) +
  ylab("Mean - p0") + xlab("Mean - Density") #+
# geom_smooth(method="lm", formula= y~x, se=FALSE)

pVISTRAPenc2 <- ggplot(data=enc, aes(x=Density,y=Mean, shape=Method, fill=Method)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23)) +
  scale_fill_manual(values=c("#00FA9A","#FFFACD")) +
  ylab("Mean - p0") +
  geom_smooth(method="lm")

### Detection probability by Abundance
det <- rbind(subset(datVIS,Type == "pstar"),subset(datTRAP,Type == "pstar"))
det$Method <- c(rep("VIS",nrow(det)/2),rep("TRAP",nrow(det)/2))
det$combo <- paste(det$Type,det$Parameter,det$Method,sep=".")

dDMean <- det %>%
  group_by(Method) %>%
  summarise(Mean = mean(Mean))

lm_fit <- lm(Mean ~ combo, data=det)

pVISTRAPdet <- ggplot(data=det, aes(x=Density,y=Mean, shape=combo, fill=combo)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(21,23,21,23,21,23,21,23)) +
  scale_fill_manual(values=c("#00a868","#c7b300","#00c77b","#ffe605","#00FA9A","#fff38a","#a3ffdc","#FFFACD")) +
  # geom_linerange(data=det, aes(ymin=Q2.5, ymax=Q97.5)) +
  ylab("Mean - pstar") #+
# geom_smooth(method="lm", formula= y~x, se=FALSE)


