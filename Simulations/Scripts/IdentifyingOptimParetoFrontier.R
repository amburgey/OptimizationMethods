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
                                ifelse(grepl("onewayTRAP",Sim),"onewayTRAP",
                                    ifelse(grepl("closedVISTRAP",Sim),"closedVISTRAP",
                                        ifelse(grepl("onewayVISTRAP",Sim),"onewayVISTRAP",
                                               ifelse(grepl("closedVIS",Sim),"closedVIS","onewayVIS")))))) %>%
  mutate(Truth = ifelse(grepl("120",Sim),"120","60")) %>%
  mutate(Scenario = sub(".* ","",Code)) %>%
  mutate(Size = ifelse(grepl("large",Sim),"large",
                          ifelse(grepl("small",Sim),"small",NA))) %>%
  mutate(Days = ifelse(grepl("14",Sim),"14","30")) %>%
  mutate(Space = ifelse(grepl("halfhalf",Sim),"halfhalf",
                               ifelse(grepl("thirdthird",Sim),"thirdthird",
                                      ifelse(grepl("half",Sim),"half","full"))))



#### Plot start-up cost scenarios plus recurring costs as compared to performance.----

## Combine start-up and recurring costs
fullcost <- inner_join(startup, perfcost2, by = c("Code2"))
fullcost <- fullcost %>%
  mutate(TotalCost = Cost.x + fullcost$Cost.y)


## Scenario One: Full start-up costs (aka, no equipment already purchased) plus recurring costs

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120FULL <- ggplot(subset(fullcost, Truth == 120), aes(x = TotalCost/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  scale_size_manual(values=c(5,3.5))+
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto120FULL

## For low BTS density
pareto60FULL <- ggplot(subset(fullcost, Truth == 60), aes(x = TotalCost/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60FULL


## Scenario Two: Half start-up costs (aka, some equipment/replacements needed) plus recurring costs

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120HALF <- ggplot(subset(fullcost, Truth == 120), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto120HALF

## For low BTS density
pareto60HALF <- ggplot(subset(fullcost, Truth == 60), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60HALF


## Scenario Three: No start-up costs (aka, all equipment purchased already) plus recurring costs (SAME AS FIRST PARETO GRAPHS)

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120NONE <- ggplot(subset(fullcost, Truth == 120), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) + 
  xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Normal density (120 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto120NONE

## For low BTS density
pareto60NONE <- ggplot(subset(fullcost, Truth == 60), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
  scale_size_manual(values=c(5,3.5)) +
  ggtitle("Low density (60 snakes)") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80)
pareto60NONE


## Sacrifical plot for legend
paretoLEGEND <- ggplot(subset(fullcost, Truth == 60), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Scenario.y, size = as.factor(Size))) +
  geom_jitter(alpha = 1) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
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


#### Create combined panel figure of normal and low density averaged across snake size ratios (large and small)

png(file="Simulations/Figures/NormalLowDensityParetoHalfStartUpAvSize.png",width=12,height=7,units="in",res=600)
(pareto60HALFAvSize + pareto120HALFAvSize)
dev.off()



#### Plot start-up cost scenarios plus recurring costs with different panels for snake size (large and small).----

#### Identify Pareto frontier.----

getPareto <- function(x,y){
  pareto = 1:length(x)
  for(i in 1:length(x)){
    cond1 = y[i]!=min(y[which(x==x[i])])
    cond2 = x[i]!=min(x[which(y==y[i])])
    for(n in 1:length(x)){
      if((x[i]>x[n]  &  y[i]>y[n]) | (x[i]==x[n] & cond1) | (y[i]==y[n] & cond2)){
        pareto[i] = NA
        break
      }
    }
  }
  #All points not on the frontier should be marked as NA in the pareto variable
  return(pareto)
}


## Scenario Two: Half start-up costs (aka, some equipment/replacements needed) plus recurring costs

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy

## For normal BTS density and small snakes
normdenssmall <- subset(fullcost, Truth == 120 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
normdenssmall$TotalCost <- (normdenssmall$TotalCost-(normdenssmall$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenssmall$TotalCost, y = normdenssmall$RMSE)
#Subset to names of Sim on frontier
pfnormsmall <- subset(fullcost, Truth == 120 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("small") & Sim %in% pfnormsmall[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pfnormsmall[[1]]))


pareto120HALFsmall <- ggplot(subset(fullcost, Truth == 120 & Size == c("small")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_jitter(width = 1, height = 1.5, seed = 20)) + 
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  # geom_text_repel(data = highlight_df %>%
  #                   mutate(label = Num),
  #                 aes(label = label),
  #                 box.padding = 0.5,
  #                 point.padding = 0.5,
  #                 show.legend = FALSE,
  #                 size = 3) +
  annotate("rect", xmin = 54, xmax = 73, ymin = 140, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. onewayVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedTRAP30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedVISTRAP30halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 150, label = "8. closedVIS30full", size = 2.5, hjust = 0)
pareto120HALFsmall


## For normal BTS density and large snakes
normdenslarge <- subset(fullcost, Truth == 120 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
normdenslarge$TotalCost <- (normdenslarge$TotalCost-(normdenslarge$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenslarge$TotalCost, y = normdenslarge$RMSE)
#Subset to names of Sim on frontier
pfnormlarge <- subset(fullcost, Truth == 120 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("large") & Sim %in% pfnormlarge[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pfnormlarge[[1]]))


pareto120HALFlarge <- ggplot(subset(fullcost, Truth == 120 & Size == c("large")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_jitter(width = 1, height = 1.5, seed = 40)) + 
  xlab("Cost in Thousands of USD") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  # geom_text_repel(data = highlight_df %>% 
  #                   mutate(label = Num),
  #                 aes(label = label), 
  #                 box.padding = 0.85,
  #                 show.legend = FALSE,
  #                 size = 3,
  #                 point.padding = 0.25) +
  annotate("rect", xmin = 54, xmax = 73, ymin = 120, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. onewayVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedTRAP30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 150, label = "8. closedVISTRAP30halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 130, label = "9. closedVIS30full", size = 2.5, hjust = 0)
pareto120HALFlarge


## For low BTS density and small snakes
lowdenssmall <- subset(fullcost, Truth == 60 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
lowdenssmall$TotalCost <- (lowdenssmall$TotalCost-(lowdenssmall$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenssmall$TotalCost, y = lowdenssmall$RMSE)
#Subset to names of Sim on frontier
pflowsmall <- subset(fullcost, Truth == 60 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("small") & Sim %in% pflowsmall[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pflowsmall[[1]]))


pareto60HALFsmall <- ggplot(subset(fullcost, Truth == 60 & Size == c("small")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_jitter(width = 1, height = 1.5, seed = 72)) + 
  ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  # geom_text_repel(data = highlight_df %>% 
  #                   mutate(label = Num),
  #                 aes(label = label), 
  #                 box.padding = 0.85,
  #                 show.legend = FALSE,
  #                 size = 3,
  #                 point.padding = 0.25) +
  annotate("rect", xmin = 54, xmax = 73, ymin = 160, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedTRAP30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60HALFsmall


## For low BTS density and large snakes
lowdenslarge <- subset(fullcost, Truth == 60 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
lowdenslarge$TotalCost <- (lowdenslarge$TotalCost-(lowdenslarge$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenslarge$TotalCost, y = lowdenslarge$RMSE)
#Subset to names of Sim on frontier
pflowlarge <- subset(fullcost, Truth == 60 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("large") & Sim %in% pflowlarge[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pflowlarge[[1]]))


pareto60HALFlarge <- ggplot(subset(fullcost, Truth == 60 & Size == c("large")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_jitter(width = 1, height = 1.5, seed = 9)) + 
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  # geom_text_repel(data = highlight_df %>% 
  #                   mutate(label = Num),
  #                 aes(label = label), 
  #                 box.padding = 0.85,
  #                 show.legend = FALSE,
  #                 size = 3,
  #                 point.padding = 0.25) +
  annotate("rect", xmin = 54, xmax = 73, ymin = 180, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60HALFlarge

## Manually place legend so able to better tweak/combine pieces

paretoLEGEND <- ggplot(subset(fullcost, Truth == 120 & Size == c("small")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.5, height = 1)) + 
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80) +
  guides(shape = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(size=5)))

# png(file="Simulations/Figures/LEGEND.png",width=12,height=7,units="in",res=600)
# paretoLEGEND
# dev.off()


#### Create combined panel figure of normal and low density averaged across snake size ratios (large and small)

png(file="Simulations/Figures/NormalLowDensityParetoHalfStartUpSmallLarge.png",width=12,height=9,units="in",res=600)
((pareto60HALFsmall + pareto120HALFsmall) / (pareto60HALFlarge + pareto120HALFlarge))
dev.off()






###### Unused figures:

#### Plot recurring costs as compared to performance with no startup.----

# ## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
# ## For normal BTS density
# pareto120 <- ggplot(subset(perfcost2, Truth == 120), aes(x = Cost, y = RMSE, fill = TypeALL, shape = Scenario, size = as.factor(Size))) + #, label=Sim
#   geom_jitter(alpha = 1) + 
#   xlab("Cost in USD") + ylab("Root Mean Square Error") +
#   theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
#   scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
#   # geom_text(position = position_jitter(width = 1, height = 1)) +
#   # geom_text(
#   # data=perfcost2 %>% filter(Truth == 120 & Code == c("VIS Half.14")), # Filter data first
#   # aes(label=Sim), position = position_jitter(width = -2, height = 3)) +
#   scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
#   scale_size_manual(values=c(5,3.5))
# pareto120

# ## For low BTS density
# pareto60 <- ggplot(subset(perfcost2, Truth == 60), aes(x = Cost, y = RMSE, fill = TypeALL, shape = Scenario, size = as.factor(Size))) + #, label=Sim
#   geom_jitter(alpha = 1) +
#   xlab("Cost in USD") + ylab("Root Mean Square Error") +
#   theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
#   scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
#   # geom_text(position = position_jitter(width = 1, height = 1), check_overlap = T)
#   # geom_text(
#   # data=perfcost2 %>% filter(Truth == 60 & Code == c("VIS Full.14")), # Filter data first
#   # aes(label=Sim), position = position_jitter(width = -2, height = 3)) +
#   scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
#   scale_size_manual(values=c(5,3.5))
# pareto60



#### Plot start-up cost scenarios plus recurring costs as compared to performance averaged across snake size ratios (large and small).----

# ## For each sampling scenario, take the mean of RMSE and total cost across large and small snake situation
# avsize <- fullcost %>%
#   dplyr::select(Cost.x,Scenario.y,TypeALL,Truth,RMSE,Size,TotalCost) %>%
#   mutate(AvScen = paste(TypeALL,Truth,Scenario.y,sep=".")) %>%
#   group_by(AvScen,Cost.x,Scenario.y,Truth,TypeALL,TotalCost) %>%
#   summarise(meanRMSE = mean(RMSE))

## Scenario Two: Half start-up costs (aka, some equipment/replacements needed) plus recurring costs

# ## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
# ## For normal BTS density
# pareto120HALFAvSize <- ggplot(subset(avsize, Truth == 120), aes(x = (TotalCost-(Cost.x/2))/1000, y = meanRMSE, fill = TypeALL, shape = Scenario.y)) +
#   geom_jitter(alpha = 1, cex = 5) + 
#   xlab("Cost in Thousands of USD") + #ylab("Root Mean Square Error") +
#   theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
#   scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
#   scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
#   ggtitle("Normal density (120 snakes)") +
#   ylim(0,(max(avsize$meanRMSE)+5)) +
#   xlim(0,80)
# pareto120HALFAvSize

# ## For low BTS density
# pareto60HALFAvSize <- ggplot(subset(avsize, Truth == 60), aes(x = (TotalCost-(Cost.x/2))/1000, y = meanRMSE, fill = TypeALL, shape = Scenario.y)) +
#   geom_jitter(alpha = 1, cex = 5) +
#   xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
#   theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
#   scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
#   scale_shape_manual(values = c(21, 22, 24, 25, 23, 23)) +
#   ggtitle("Low density (60 snakes)") +
#   ylim(0,(max(fullcost$RMSE)+5)) +
#   xlim(0,80)
# pareto60HALFAvSize


## Bolded points?
#
# ggplot(subset(fullcost, Truth == 60 & Size == c("large")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
#   geom_point(alpha = 0.7, position = position_jitter(width = 0.5, height = 1)) + #, aes(colour=factor(Sim))) +
#   geom_point(data = highlight_df, aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE), color = "black", stroke = 1.5) +
#   theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
#   scale_fill_manual(values = c("#799D31","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
#   scale_shape_manual(values = c(21, 22, 24, 25)) +
#   scale_size_manual(values=c(3.5,5)) +
#   # scale_color_manual(values = cols) +
#   ggtitle("Low density (60) and majority large snakes") +
#   ylim(0,(max(fullcost$RMSE)+5)) +
#   xlim(0,80)

