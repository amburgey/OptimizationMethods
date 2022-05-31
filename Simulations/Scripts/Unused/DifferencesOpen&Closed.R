### Compare closed to open SCR analysis
### Redo when have 100 simulations

library(ggplot2); library(HDInterval)

nsims <- 50

## read in "closed" pop analysis
resC <- matrix(NA, nrow = nsims, ncol = 6, dimnames = list(seq(1:nsims),c("MnN","Mnp01","Mnp02","Mnp03","Mnp04","sigma")))

for(i in 1:nsims){
  load(paste("Simulations/Results/RESULTS_closedVIS120small60half",i,".Rdata",sep=""))
  print(i)
  resC[i,1] <- out$mean$N
  resC[i,2:5] <- out$mean$p0
  resC[i,6] <- out$mean$sigma
}
resC <- t(resC)

## read in "open" pop analysis
resO <- matrix(NA, nrow = nsims, ncol = 6, dimnames = list(seq(1:nsims),c("MnN","Mnp01","Mnp02","Mnp03","Mnp04","sigma")))

for(i in 1:nsims){
  load(paste("Simulations/Results/RESULTS_openVIS120small60half",i,".Rdata",sep=""))
  print(i)
  resO[i,1] <- out$mean$N
  resO[i,2:5] <- out$mean$p0
  resO[i,6] <- out$mean$sigma
}

resO <- t(resO)

res <- rbind(resC,resO)

## calculate mean and stdeviation to check if too far apart
newdf <- matrix(NA, nrow = 12, ncol = 4, dimnames = list(c("MnN.c","Mnp01.c","Mnp02.c","Mnp03.c","Mnp04.c","sigma.c","MnN.o","Mnp01.o","Mnp02.o","Mnp03.o","Mnp04.o","sigma.o"),c("Mn","SD","2.5HDI","97.5HDI")))
newdf[,1] <- apply(res, 1, mean)
newdf[,2] <- apply(res, 1, sd)
newdf[,3] <- apply(res, 1, hdi)[1,]
newdf[,4] <- apply(res, 1, hdi)[2,]
newdf <- cbind(newdf, as.data.frame(c("Abundance.c","p01.c","p02.c","p03.c","p04.c","sigma.c","Abundance.o","p01.o","p02.o","p03.o","p04.o","sigma.o")))
newdf <- cbind(newdf, as.data.frame(rep(c("closed","open"), each = 6)))
newdf <- cbind(newdf, as.data.frame(c("Abundance",rep("p0",times=4),"sigma","Abundance",rep("p0",times=4),"sigma")))
colnames(newdf)[5:7] <- c("Par","Study Area Type","Parameter")

## Plot mean and 1 stdeviation within the mean to see if overlap between analyses

p1 <- ggplot(newdf, aes(x=Par, y = Mn, shape = `Study Area Type`, color = Parameter)) + 
  geom_point(size = 6, stroke = 2) +
  scale_shape_manual(values = c(19,21)) +
  facet_wrap(. ~ Parameter, scales = "free") +
  geom_errorbar(aes(ymin = `2.5HDI`, ymax = `97.5HDI`), width = 0.2, size = 0.9) +
  ylab("Mean Estimate") +
  scale_color_manual(values = c("#f48d79","#02adb5","#ab0902")) +
  theme(axis.text = element_text(size = 12), axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14), legend.text = element_text(size = 14),
        legend.title = element_text(size = 16), legend.position = "top",
        axis.title.x = element_blank(), axis.ticks.x = element_blank())

png(file="Open&Closedcomparison.png",width=15,height=8,units="in",res=600)
p1
dev.off()

