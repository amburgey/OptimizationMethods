### Plotting monitoring scenarios by cost
### Highlighting Pareto frontier where a scenario results in a gain in the monitoring objective (i.e., reducing RMSE of abundance) with no gain or a reduction in cost

rm(list=ls())

library(ggplot2);library(tidyverse);library(dplyr)


#### Costs of scenarios.----

cost <- read_csv('Simulations/simDat/AlternateScenarioCosts.csv')

### Measures of performance in achieving monitoring objective.----

perf <- read_csv('Simulations/Results/SimResultsSummary.csv')[,-1]

#### Plot initial/start-up cost comparison.----

lvls <- c("VIS", "TRAP Full","TRAP Half","VISTRAP Half","VISTRAP Third")

startup <- cost %>%
  select("TotalStart","Scenario") %>%
  distinct() %>%
  add_row(TotalStart = c("$2,000"), Scenario = c("VIS")) %>%
  add_row(TotalStart = c("$18,900"), Scenario = c("VISTRAP Half")) %>%  ## create combined method cost
  add_row(TotalStart = c("$13,700"), Scenario = c("VISTRAP Third")) %>%  ## create combined method cost
  slice(-c(1,3,5:8),) %>%  # remove separate costs of combined method scenarios
  mutate(Scenario = factor(Scenario, levels = lvls)) %>%  #reorder levels
  mutate(Type = c("TRAP","TRAP","VIS","VISTRAP","VISTRAP")) %>%
  mutate(Cost = gsub(",","",TotalStart)) %>%
  mutate(Cost = as.numeric(gsub("\\$","",Cost)))
  


startplot <- ggplot(startup, aes(x = Scenario, y = Cost, fill = Type)) +
  geom_point(shape = 21, cex = 4) +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("TRAP","VIS","VISTRAP"))

#### Combine performance and cost datasets to figure out value for estimate achieved.----

recur <- cost %>%
  select("TotalRun","Scenario","Nights") %>%
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
  mutate(Code = paste(Type,Sched,Nights,sep=""))

perfcost <- inner_join(recur, perf, by = "Code")
perfcost2 <- perfcost %>%
  mutate(Cost = gsub(",","",TotalRun)) %>%
  mutate(Cost = as.numeric(gsub("\\$","",Cost))) %>%
  mutate(TypeALL = ifelse(grepl("closedTRAP",Sim),"closedTRAP",
                          ifelse(grepl("closedVIS",Sim),"closedVIS",
                                 ifelse(grepl("onewayTRAP",Sim),"onewayTRAP",
                                        ifelse(grepl("onewayVIS",Sim),"onewayVIS",
                                               ifelse(grepl("closedVISTRAP",Sim),"closedVISTRAP","onewayVISTRAP")))))) %>%
  mutate(Truth = ifelse(grepl("120",Sim),"120","60"))

### PLOT ALL TOGETHER OR SEPARATE DIFFERENT DENSITIES AND SIZE DISTRIBUTIONS ONTO DIFFERENT PANELS?

#### Plot recurring costs as compared to performance.----

## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy
## For normal BTS density
pareto120 <- ggplot(subset(perfcost2, Truth == 120), aes(x = Cost, y = RMSE, fill = TypeALL, label=Sim)) +
  geom_point(shape = 21, cex = 5) +
  xlab("Cost in USD") + ylab("Root Mean Square Error") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayVIS")) +
  geom_text(hjust=-0.1,vjust=0)
pareto120

## For low BTS density
pareto60 <- ggplot(subset(perfcost2, Truth == 60), aes(x = Cost, y = RMSE, fill = TypeALL, label=Sim)) +
  geom_point(shape = 21, cex = 5) +
  xlab("Cost in USD") + ylab("Root Mean Square Error") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("closedTRAP","closedVIS","onewayVIS")) +
  geom_text(hjust=-0.1,vjust=0)
pareto60

## Percent coefficient of Variation (%CV) - percent of variability of posterior around the mean, lower value indicates less variability
pareto2 <- ggplot(perfcost, aes(x = TotalRun, y = `%CV`, fill = Type.x)) +
  geom_point(shape = 21, cex = 3) +
  xlab("Cost in USD") + ylab("% Coefficient of Variation") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("TRAP","VIS","VISTRAP"))
pareto2

## Percent bias (%Bias) - each estimate subtracted from the true abundance, value closer to zero indicates less bias
pareto3 <- ggplot(perfcost, aes(x = TotalRun, y = `%Bias`, fill = Type.x)) +
  geom_point(shape = 21, cex = 3) +
  xlab("Cost in USD") + ylab("% Bias") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("TRAP","VIS","VISTRAP"))
pareto3

## Percent coverage (%Coverage) - percent of simulations where the highest density posterior interval (HDPI) included the true abundance, closer to 100% indicates better coverage
pareto4 <- ggplot(perfcost, aes(x = TotalRun, y = `%Coverage`, fill = Type.x)) +
  geom_point(shape = 21, cex = 3) +
  xlab("Cost in USD") + ylab("% Coverage") +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("TRAP","VIS","VISTRAP"))
pareto4

