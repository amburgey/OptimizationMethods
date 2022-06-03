## Supplementary figure of visual surveys (VIS) and trapping surveys (TRAP)
## The purpose of this code is to create a quick plot to compare encounter probabilities of different snake size categories from real data analysis

library(ggplot2); library(jagsUI); library(coda);library(runjags)


## CREATE COMBINED SIGMA FROM MODEL RESULTS FROM VIS AND TRAP ANALYSIS

#### VISUAL SURVEY REAL DATA RESULTS ####
load("Real Data Analysis/Visual surveys/Results/NWFNVISALL_SCRpstar.Rdata")  # unified analysis of all datasets

## Truncate posterior in order to have longer burn-in, have to do separately for each chain and then recombine to single posterior for parameters of interest
mcmc.object <- as.mcmc.list(out$samples[[1]])
out1 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[2]])
out2 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[3]])
out3 <- window(mcmc.object, start=1499)
outALLV <- combine.mcmc(list(out1,out2,out3), thin = 1)
niter(outALLV)
postV <- as.matrix(outALLV)


#### TRAPPING SURVEY REAL DATA RESULTS.----
load("Real Data Analysis/Trapping/Results/NWFNTRAPALL_SCRpstar.Rdata")

## Truncate posterior in order to have longer burn-in, have to do separately for each chain and then recombine to single posterior for parameters of interest
mcmc.object <- as.mcmc.list(out$samples[[1]])
out1 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[2]])
out2 <- window(mcmc.object, start=1499)
mcmc.object <- as.mcmc.list(out$samples[[3]])
out3 <- window(mcmc.object, start=1499)
outALLT <- combine.mcmc(list(out1,out2,out3), thin = 1)
niter(outALLT)
postT <- as.matrix(outALLT)


##### PLOT TWO PANEL VIS AND TRAP ENCOUNTER PROBABILITIES.----

png(file="Real Data Analysis/Figures/TRAP&VISencounter.png",width=10,height=8,units="in",res=600)

par(mfrow = c(2, 1))

hist(postT[,1], col = alpha("#26580F",0.4), breaks = 25,
     xlim = c(0.0005,0.0055), ylim = c(0,5800), xlab = "Encounter Probability", ylab = "Frequency",
     main = NULL)
hist(postT[,2], col = alpha("#378805",0.4), add = TRUE, breaks = 25)
hist(postT[,3], col = alpha("#86DC3D",0.4), add = TRUE, breaks = 25)
hist(postT[,4], col = alpha("#C5E90B",0.4), add = TRUE, breaks = 25)
text(0.0017, 5700, expression("Trapping Surveys"), cex = 1.5)
text(0.0016, 5000, expression("< 850"))
text(0.0026, 3500, expression("850-950"))
text(0.004, 3100, expression("950-1150"))
text(0.0045, 3900, expression("> 1150"))

hist(postV[,1], col = alpha("#451A00",0.4), breaks = 25,
     xlim = c(0.0005,0.0055), ylim = c(0,5800), xlab = "Encounter Probability", ylab = "Frequency",
     main = NULL)
hist(postV[,2], col = alpha("#A66A2E",0.4), add = TRUE, breaks = 25)
hist(postV[,3], col = alpha("#C48A47",0.4), add = TRUE, breaks = 25)
hist(postV[,4], col = alpha("#F0C648",0.4), add = TRUE, breaks = 25)
text(0.0017, 5700, expression("Visual Surveys"), cex = 1.5)
text(0.003, 4300, expression("< 850"))
text(0.0042, 2900, expression("850-950"))
text(0.0036, 3800, expression("950-1150"))
text(0.0024, 2900, expression("> 1150"))
segments(0.0042, 2700, 0.0037, 2000,
         col = "black", lty = 1, lwd = 1)
segments(0.0024, 2700, 0.0026, 2000,
         col = "black", lty = 1, lwd = 1)

dev.off()
