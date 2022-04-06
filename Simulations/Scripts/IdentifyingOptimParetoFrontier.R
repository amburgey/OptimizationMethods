### Plotting monitoring scenarios by cost
### Highlighting Pareto frontier where a scenario results in a gain in the monitoring objective (i.e., reducing RMSE of abundance) with no gain or a reduction in cost

rm(list=ls())

library(ggplot2);library(tidyverse);library(dplyr);library(ggrepel);library(patchwork)


#### Costs of scenarios.----

cost <- read_csv('Simulations/simDat/AlternateScenarioCosts.csv')

### Measures of performance in achieving monitoring objective.----

perf <- read_csv('Simulations/Results/SimResultsSummary.csv')[,-1]

#### Plot initial/start-up cost comparison.----

lvls <- c("VIS", "TRAP Full","TRAP Half","VISTRAP Half","VISTRAP Third")

startup <- cost %>%
  dplyr::select("TotalStart","Scenario") %>%
  distinct() %>%
  add_row(TotalStart = c("$2,000"), Scenario = c("VIS")) %>%
  add_row(TotalStart = c("$18,900"), Scenario = c("VISTRAP Half")) %>%  ## create combined method cost
  add_row(TotalStart = c("$13,700"), Scenario = c("VISTRAP Third")) %>%  ## create combined method cost
  slice(-c(1,3,5:8),) %>%  # remove separate costs of combined method scenarios
  mutate(Code2 = factor(Scenario, levels = lvls)) %>%  #reorder levels
  mutate(Type = c("TRAP","TRAP","VIS","VISTRAP","VISTRAP")) %>%
  mutate(Cost = gsub(",","",TotalStart)) %>%
  mutate(Cost = as.numeric(gsub("\\$","",Cost)))
  


startplot <- ggplot(startup, aes(x = Scenario, y = Cost, fill = Type)) +
  geom_point(shape = 21, cex = 4) +
  ylab("Startup Cost") + xlab("Monitoring Scenario") +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("TRAP","VIS","VISTRAP"))

#### Combine performance and cost datasets to figure out value for estimate achieved.----

recur <- cost %>%
  dplyr::select("TotalRun","Scenario","Nights") %>%
  filter(Nights != 60) %>%  #not including 60 days anymore
  drop_na(TotalRun) %>%
  mutate(Type = ifelse(grepl("VIS Half +|VIS Third +", Scenario),"VISTRAP",
                       ifelse(grepl("VIS ", Scenario),"VIS",
                              ifelse(grepl("TRAP ", Scenario),"TRAP","9999")))) %>%
  mutate(Code = paste(Scenario,Nights,sep=".")) %>%
  mutate(Code = str_replace(Code, "VIS Half [+]", "VISTRAP Half")) %>%
  mutate(Code = str_replace(Code, "VIS Third [+]", "VISTRAP Third"))

perf <- perf %>%
  mutate(Type = ifelse(grepl("VISTRAP", Sim),"VISTRAP ",
                       ifelse(grepl("TRAP",Sim),"TRAP ","VIS "))) %>%
  mutate(Sched = ifelse(grepl("full", Sim),"Full.",
                       ifelse(grepl("half",Sim),"Half.","Third."))) %>%
  mutate(Nights = ifelse(grepl("l60|e60", Sim),"60",
                         ifelse(grepl("l30|e30",Sim),"30","14"))) %>%
  mutate(Code = paste(Type,Sched,Nights,sep="")) %>%
  mutate(Code2 = paste(Type,gsub("\\.","",Sched),sep="")) %>%
  mutate(Code2 = ifelse(grepl("VIS Full",Code2),"VIS",
                        ifelse(grepl("VIS Half",Code2),"VIS",Code2)))

perfcost <- inner_join(recur, perf, by = "Code")
perfcost2 <- perfcost %>%
  mutate(Cost = gsub(",","",TotalRun)) %>%
  mutate(Cost = as.numeric(gsub("\\$","",Cost))) %>%
  mutate(TypeALL = ifelse(grepl("closedTRAP",Sim),"closedTRAP",
                          ifelse(grepl("closedVIS",Sim),"closedVIS",
                                 ifelse(grepl("onewayTRAP",Sim),"onewayTRAP",
                                        ifelse(grepl("onewayVIS",Sim),"onewayVIS",
                                               ifelse(grepl("closedVISTRAP",Sim),"closedVISTRAP","onewayVISTRAP")))))) %>%
  mutate(Truth = ifelse(grepl("120",Sim),"120","60")) %>%
  mutate(Scenario = sub(".* ","",Code)) %>%
  mutate(Size = ifelse(grepl("large",Sim),"large",
                          ifelse(grepl("small",Sim),"small",NA)))


#### Plot recurring costs as compared to performance.----

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120 <- ggplot(subset(perfcost2, Truth == 120), aes(x = Cost, y = RMSE, fill = TypeALL, shape = Scenario, size = as.factor(Size))) + #, label=Sim
  geom_jitter(alpha = 0.75) + 
  xlab("Cost in USD") + ylab("Root Mean Square Error") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  # geom_text(position = position_jitter(width = 1, height = 1)) +
  # geom_text(
  # data=perfcost2 %>% filter(Truth == 120 & Code == c("VIS Half.14")), # Filter data first
  # aes(label=Sim), position = position_jitter(width = -2, height = 3)) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5))
pareto120

## For low BTS density
pareto60 <- ggplot(subset(perfcost2, Truth == 60), aes(x = Cost, y = RMSE, fill = TypeALL, shape = Scenario, size = as.factor(Size))) + #, label=Sim
  geom_jitter(alpha = 0.75) +
  xlab("Cost in USD") + ylab("Root Mean Square Error") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  # geom_text(position = position_jitter(width = 1, height = 1), check_overlap = T)
  # geom_text(
  # data=perfcost2 %>% filter(Truth == 60 & Code == c("VIS Full.14")), # Filter data first
  # aes(label=Sim), position = position_jitter(width = -2, height = 3)) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5))
pareto60



#### Plot start-up cost scenarios plus recurring costs as compared to performance.----

## Combine start-up and recurring costs
fullcost <- inner_join(startup, perfcost2, by = c("Code2"))
fullcost <- fullcost %>%
  mutate(TotalCost = Cost.x + fullcost$Cost.y)


## Scenario One: Full start-up costs (aka, no equipment already purchased) plus recurring costs

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120FULL <- ggplot(subset(fullcost, Truth == 120), aes(x = TotalCost/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5))+
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto120FULL

## For low BTS density
pareto60FULL <- ggplot(subset(fullcost, Truth == 60), aes(x = TotalCost/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60FULL


## Scenario Two: Half start-up costs (aka, some equipment/replacements needed) plus recurring costs

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120HALF <- ggplot(subset(fullcost, Truth == 120), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto120HALF

## For low BTS density
pareto60HALF <- ggplot(subset(fullcost, Truth == 60), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60HALF


## Scenario Three: No start-up costs (aka, all equipment purchased already) plus recurring costs (SAME AS FIRST PARETO GRAPHS)

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120NONE <- ggplot(subset(fullcost, Truth == 120), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto120NONE

## For low BTS density
pareto60NONE <- ggplot(subset(fullcost, Truth == 60), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60NONE


## Sacrifical plot for legend
paretoLEGEND <- ggplot(subset(fullcost, Truth == 60), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 0.75) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(5,3.5)) +
  guides(shape = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(size=5)))



#### Create combined panel figure of normal and low density

png(file="Simulations/Figures/NormalLowDensityParetoNoStartUp.png",width=12,height=7,units="in",res=600)
(pareto60NONE + pareto120NONE)
dev.off()

png(file="Simulations/Figures/NormalLowDensityParetoHalfStartUp.png",width=12,height=7,units="in",res=600)
(pareto60HALF + pareto120HALF)
dev.off()

png(file="Simulations/Figures/NormalLowDensityParetoFullStartUp.png",width=12,height=7,units="in",res=600)
(pareto60FULL + pareto120FULL)
dev.off()

## Manually place legend so able to better tweak/combine pieces
# png(file="Simulations/Figures/LEGEND.png",width=12,height=7,units="in",res=600)
# paretoLEGEND
# dev.off()



#### Plot start-up cost scenarios plus recurring costs as compared to performance averaged across snake size ratios (large and small).----

## For each sampling scenario, take the mean of RMSE and total cost across large and small snake situation
avsize <- fullcost %>%
  dplyr::select(Cost.x,Scenario.y,TypeALL,Truth,RMSE,Size,TotalCost) %>%
  mutate(AvScen = paste(TypeALL,Truth,Scenario.y,sep=".")) %>%
  group_by(AvScen,Cost.x,Scenario.y,Truth,TypeALL,TotalCost) %>%
  summarise(meanRMSE = mean(RMSE))

## Scenario Two: Half start-up costs (aka, some equipment/replacements needed) plus recurring costs

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120HALFAvSize <- ggplot(subset(avsize, Truth == 120), aes(x = (TotalCost-(Cost.x/2))/1000, y = meanRMSE, fill = TypeALL, shape = Scenario.y)) +
  geom_jitter(alpha = 0.75, cex = 5) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(avsize$meanRMSE)+5)) +
  xlim(0,80)
pareto120HALFAvSize

## For low BTS density
pareto60HALFAvSize <- ggplot(subset(avsize, Truth == 60), aes(x = (TotalCost-(Cost.x/2))/1000, y = meanRMSE, fill = TypeALL, shape = Scenario.y)) +
  geom_jitter(alpha = 0.75, cex = 5) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#C9F277","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayTRAP","onewayVIS")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60HALFAvSize



#### Create combined panel figure of normal and low density averaged across snake size ratios (large and small)

png(file="Simulations/Figures/NormalLowDensityParetoHalfStartUpAvSize.png",width=12,height=7,units="in",res=600)
(pareto60HALFAvSize + pareto120HALFAvSize)
dev.off()

