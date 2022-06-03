### Plotting root mean square error (RMSE) of each monitoring scenario by the cost of conducting it
### Highlighting the Pareto frontier where a scenario results in a gain on the monitoring objective (i.e., reducing RMSE of abundance) with a reduction or no change in cost
### Density of transect sampling can be full (all transects) or half when using a single method (VIS or TRAP) or half and half or third and third when using combined methods (VISTRAP) 

rm(list=ls())

library(ggplot2);library(tidyverse);library(dplyr);library(patchwork)


#### COSTS OF SCENARIOS.----

cost <- read_csv('Simulations/simDat/AlternateScenarioCosts.csv')  # ignoring warning about duplicated column names


### METRICS OF PERFORMANCE (E.G., RMSE).----

perf <- read_csv('Simulations/Results/SimResultsSummary.csv')[,-1]


#### CREATE DATAFRAME OF STARTUP COSTS.----

lvls <- c("VIS", "TRAP Full","TRAP Half","VISTRAP Half","VISTRAP Third")

startup <- cost %>%
  dplyr::select("TotalStart","Scenario") %>%
  distinct() %>%
  add_row(TotalStart = c("$2,000"), Scenario = c("VIS")) %>%             # visual surveys don't change in base startup cost as equipment is the same
  add_row(TotalStart = c("$18,900"), Scenario = c("VISTRAP Half")) %>%   # create combined method startup cost as increased number of traps changes the cost
  add_row(TotalStart = c("$13,700"), Scenario = c("VISTRAP Third")) %>%  # create combined method startup cost as increased number of traps changes the cost
  slice(-c(1,3,5:8),) %>%                                                # remove separated costs of combined method scenarios to just retain the combined cost
  mutate(Code2 = factor(Scenario, levels = lvls)) %>%                    # reorder levels
  mutate(Type = c("TRAP","TRAP","VIS","VISTRAP","VISTRAP")) %>%          # type of sampling scenario
  mutate(Cost = gsub(",","",TotalStart)) %>%                             # remove comma from cost
  mutate(Cost = as.numeric(gsub("\\$","",Cost)))                         # remove dollar sign from cost
  
## Basic plot to see different costs in comparison
startplot <- ggplot(startup, aes(x = Scenario, y = Cost, fill = Type)) +
  geom_point(shape = 21, cex = 4) +
  ylab("Startup Cost") + xlab("Monitoring Scenario") +
  scale_fill_manual(values = c("#799D31","#B77A29","#FDEC9E"), breaks = c("TRAP","VIS","VISTRAP"))


#### CREATE DATAFRAME OF RECURRING COSTS.----

# Recurring costs based on monitoring scenario
recur <- cost %>%
  dplyr::select("TotalRun","Scenario","Nights") %>%
  filter(Nights != 60) %>%                                                            # not including 60 days anymore
  drop_na(TotalRun) %>%                                                               # remove separated costs of combined method scenarios to just retain the combined cost
  mutate(Type = ifelse(grepl("VIS Half +|VIS Third +", Scenario),"VISTRAP",           # rename scenarios to properly combine scenario names
                       ifelse(grepl("VIS ", Scenario),"VIS",
                              ifelse(grepl("TRAP ", Scenario),"TRAP","9999")))) %>%
  mutate(Code = paste(Scenario,Nights,sep=".")) %>%                                   # create code to describe method, number of transects, and number of days
  mutate(Code = str_replace(Code, "VIS Half [+]", "VISTRAP Half")) %>%                # rename scenarios to properly combine scenario names
  mutate(Code = str_replace(Code, "VIS Third [+]", "VISTRAP Third"))                  # rename scenarios to properly combine scenario names


#### SHAPE PERFORMANCE DATAFRAME FOR PLOTTING.----

perf <- perf %>%
  mutate(Type = ifelse(grepl("VISTRAP", Sim),"VISTRAP ",                     # separate out method from simulation name
                       ifelse(grepl("TRAP",Sim),"TRAP ","VIS "))) %>%
  mutate(Sched = ifelse(grepl("full", Sim),"Full.",                          # separate out number of transects from simulation name
                       ifelse(grepl("half",Sim),"Half.","Third."))) %>%
  mutate(Nights = ifelse(grepl("l60|e60", Sim),"60",                         # separate out number of days from simulation name
                         ifelse(grepl("l30|e30",Sim),"30","14"))) %>%
  mutate(Code = paste(Type,Sched,Nights,sep="")) %>%                         # create code to describe method, number of transects, and number of days
  mutate(Code2 = paste(Type,gsub("\\.","",Sched),sep="")) %>%                # separate out method and transects for plotting
  mutate(Code2 = ifelse(grepl("VIS Full",Code2),"VIS",                       # separate out method and transects for plotting
                        ifelse(grepl("VIS Half",Code2),"VIS",Code2)))        # get method without space for plotting


#### COMBINE PERFORMANCE METRIC (E.G., RMSE) AND COST.----

## Combine performance and cost by Code
perfcost <- inner_join(recur, perf, by = "Code")
## Mutate dataframe for plotting
perfcost2 <- perfcost %>%
  mutate(Cost = gsub(",","",TotalRun)) %>%                                   # remove comma from cost
  mutate(Cost = as.numeric(gsub("\\$","",Cost))) %>%                         # remove dollar sign from cost
  mutate(TypeALL = ifelse(grepl("closedTRAP",Sim),"closedTRAP",              # separate out method and barrier from simulation name
                                ifelse(grepl("onewayTRAP",Sim),"onewayTRAP",
                                    ifelse(grepl("closedVISTRAP",Sim),"closedVISTRAP",
                                        ifelse(grepl("onewayVISTRAP",Sim),"onewayVISTRAP",
                                               ifelse(grepl("closedVIS",Sim),"closedVIS","onewayVIS")))))) %>%
  mutate(Truth = ifelse(grepl("120",Sim),"120","60")) %>%                    # separate out number of days from simulation name
  mutate(Scenario = sub(".* ","",Code)) %>%                                  # get number of transects and days from simulation name
  mutate(Size = ifelse(grepl("large",Sim),"large",                           # get size of snakes from simulation name
                          ifelse(grepl("small",Sim),"small",NA))) %>%
  mutate(Days = ifelse(grepl("14",Sim),"14","30")) %>%                       # get number of days from simulation name
  mutate(Space = ifelse(grepl("halfhalf",Sim),"halfhalf",                    # get transect sampling scheme from simulation name
                               ifelse(grepl("thirdthird",Sim),"thirdthird",
                                      ifelse(grepl("half",Sim),"half","full"))))


#### CALCULATE FULL COST (STARTUP PLUS RECURRING COSTS).----

## Combine start-up and recurring costs to get full cost
fullcost <- inner_join(startup, perfcost2, by = c("Code2"))
fullcost <- fullcost %>%
  mutate(TotalCost = Cost.x + fullcost$Cost.y)


#### FIGURE - PLOT HALF START-UP COST SCENARIOS PLUS RECURRING COSTS AS COMPARED TO RMSE WITH DIFFERENT PANELS FOR SNAKE SIZE (SMALL AND LARGE).----

## Function to identify Pareto frontier

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


## Scenario: Half start-up costs (aka, some equipment/replacements needed) plus recurring costs
## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy


## PANEL ONE. LOW DENSITY (60), SMALL SNAKES
lowdenssmall <- subset(fullcost, Truth == 60 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
lowdenssmall$TotalCost <- (lowdenssmall$TotalCost-(lowdenssmall$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenssmall$TotalCost, y = lowdenssmall$RMSE)
#Subset to names of Sim on frontier
pflowsmall <- subset(fullcost, Truth == 60 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("small") & Sim %in% pflowsmall[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pflowsmall[[1]]))

pareto60HALFsmall <- ggplot(subset(fullcost, Truth == 60 & Size == c("small")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  ## Pareto efficient actions identified above
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


## PANEL TWO. NORMAL DENSITY (120), SMALL SNAKES
normdenssmall <- subset(fullcost, Truth == 120 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
normdenssmall$TotalCost <- (normdenssmall$TotalCost-(normdenssmall$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenssmall$TotalCost, y = normdenssmall$RMSE)
#Subset to names of simulations on the frontier
pfnormsmall <- subset(fullcost, Truth == 120 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("small") & Sim %in% pfnormsmall[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pfnormsmall[[1]]))

pareto120HALFsmall <- ggplot(subset(fullcost, Truth == 120 & Size == c("small")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  ## Pareto efficient actions identified above
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


## PANEL THREE. LOW DENSITY (60), LARGE SNAKES
lowdenslarge <- subset(fullcost, Truth == 60 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
lowdenslarge$TotalCost <- (lowdenslarge$TotalCost-(lowdenslarge$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenslarge$TotalCost, y = lowdenslarge$RMSE)
#Subset to names of Sim on frontier
pflowlarge <- subset(fullcost, Truth == 60 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("large") & Sim %in% pflowlarge[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pflowlarge[[1]]))

pareto60HALFlarge <- ggplot(subset(fullcost, Truth == 60 & Size == c("large")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 54, xmax = 73, ymin = 180, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60HALFlarge


## PANEL FOUR. NORMAL DENSITY (120), LARGE SNAKES
normdenslarge <- subset(fullcost, Truth == 120 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to half start up and divide by 1000 to put in terms of thousands of USD
normdenslarge$TotalCost <- (normdenslarge$TotalCost-(normdenslarge$Cost.x/2))/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenslarge$TotalCost, y = normdenslarge$RMSE)
#Subset to names of Sim on frontier
pfnormlarge <- subset(fullcost, Truth == 120 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("large") & Sim %in% pfnormlarge[[1]]) %>%
  arrange((TotalCost-(Cost.x/2))/1000) %>%
  mutate(Num = 1:length(pfnormlarge[[1]]))

pareto120HALFlarge <- ggplot(subset(fullcost, Truth == 120 & Size == c("large")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  xlab("Cost in Thousands of USD") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,75) +
  ## Pareto efficient actions identified above
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

## Create plot for legend to better tweak/combine pieces manually outside of R

paretoLEGEND <- ggplot(subset(fullcost, Truth == 120 & Size == c("small")), aes(x = (TotalCost-(Cost.x/2))/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_jitter(width = 0.5, height = 1)) + 
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(0,80) +
  guides(shape = guide_legend(override.aes = list(size=5)), fill = guide_legend(override.aes = list(size=5)))


#### WRITE PLOTS TO MULTI-PANEL PNG FIGURE.----

## Numeric labels and arrows of individual Pareto frontier points done manually outside R to prevent overlap

png(file="Simulations/Figures/NormalLowDensityParetoHalfStartUpSmallLarge.png",width=12,height=9,units="in",res=600)
((pareto60HALFsmall + pareto120HALFsmall) / (pareto60HALFlarge + pareto120HALFlarge))
dev.off()

#### WRITE PLOT TO PNG FIGURE IN ORDER TO HARVEST COMBINED LEGEND AND MANUALLY COMBINE.----

png(file="Simulations/Figures/LEGEND.png",width=12,height=7,units="in",res=600)
paretoLEGEND
dev.off()



#### APPENDIX - NO STARTUP.----

## Scenario: No start-up costs (aka, all equipment already purchased) plus recurring costs
## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy


## PANEL ONE. LOW DENSITY (60), SMALL SNAKES
lowdenssmall <- subset(fullcost, Truth == 60 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to no start up and divide by 1000 to put in terms of thousands of USD
lowdenssmall$TotalCost <- (lowdenssmall$TotalCost-lowdenssmall$Cost.x)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenssmall$TotalCost, y = lowdenssmall$RMSE)
#Subset to names of Sim on frontier
pflowsmall <- subset(fullcost, Truth == 60 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("small") & Sim %in% pflowsmall[[1]]) %>%
  arrange((TotalCost-Cost.x)/1000) %>%
  mutate(Num = 1:length(pflowsmall[[1]]))

pareto60NOsmall <- ggplot(subset(fullcost, Truth == 60 & Size == c("small")), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(10,75) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 54, xmax = 74, ymin = 140, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. closedTRAP14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedTRAP30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVISTRAP14halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 150, label = "8. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60NOsmall


## PANEL TWO. NORMAL DENSITY (120), SMALL SNAKES
normdenssmall <- subset(fullcost, Truth == 120 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to no start up and divide by 1000 to put in terms of thousands of USD
normdenssmall$TotalCost <- (normdenssmall$TotalCost-normdenssmall$Cost.x)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenssmall$TotalCost, y = normdenssmall$RMSE)
#Subset to names of Sim on frontier
pfnormsmall <- subset(fullcost, Truth == 120 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("small") & Sim %in% pfnormsmall[[1]]) %>%
  arrange((TotalCost-Cost.x)/1000) %>%
  mutate(Num = 1:length(pfnormsmall[[1]]))

pareto120NOsmall <- ggplot(subset(fullcost, Truth == 120 & Size == c("small")), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(10,75) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 54, xmax = 74, ymin = 100, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. onewayTRAP14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. onewayVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedTRAP30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVISTRAP14halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 150, label = "8. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 130, label = "9. closedVISTRAP30halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 110, label = "10. closedVIS30full", size = 2.5, hjust = 0)
pareto120NOsmall


## PANEL THREE. LOW DENSITY (60), LARGE SNAKES
lowdenslarge <- subset(fullcost, Truth == 60 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to no start up and divide by 1000 to put in terms of thousands of USD
lowdenslarge$TotalCost <- (lowdenslarge$TotalCost-lowdenslarge$Cost.x)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenslarge$TotalCost, y = lowdenslarge$RMSE)
#Subset to names of Sim on frontier
pflowlarge <- subset(fullcost, Truth == 60 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("large") & Sim %in% pflowlarge[[1]]) %>%
  arrange((TotalCost-Cost.x)/1000) %>%
  mutate(Num = 1:length(pflowlarge[[1]]))

pareto60NOlarge <- ggplot(subset(fullcost, Truth == 60 & Size == c("large")), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(10,75) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 54, xmax = 74, ymin = 140, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. closedTRAP14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedVISTRAP14halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 150, label = "8. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60NOlarge


## PANEL FOUR. NORMAL DENSITY (120), LARGE SNAKES
normdenslarge <- subset(fullcost, Truth == 120 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to no start up and divide by 1000 to put in terms of thousands of USD
normdenslarge$TotalCost <- (normdenslarge$TotalCost-normdenslarge$Cost.x)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenslarge$TotalCost, y = normdenslarge$RMSE)
#Subset to names of Sim on frontier
pfnormlarge <- subset(fullcost, Truth == 120 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("large") & Sim %in% pfnormlarge[[1]]) %>%
  arrange((TotalCost-Cost.x)/1000) %>%
  mutate(Num = 1:length(pfnormlarge[[1]]))

pareto120NOlarge <- ggplot(subset(fullcost, Truth == 120 & Size == c("large")), aes(x = (TotalCost-Cost.x)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.8)) +
  xlab("Cost in Thousands of USD") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(10,75) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 54, xmax = 74, ymin = 80, ymax = 320, alpha = 0.2) +
  annotate("text", x = 55, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 55, y = 290, label = "1. closedTRAP14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 270, label = "2. onewayVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 250, label = "3. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 230, label = "4. closedTRAP30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 210, label = "5. closedVISTRAP14halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 190, label = "6. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 170, label = "7. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 150, label = "8. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 130, label = "9. closedTRAP30full", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 110, label = "10. closedVISTRAP30halfhalf", size = 2.5, hjust = 0) +
  annotate("text", x = 55, y = 90, label = "11. closedVIS30full", size = 2.5, hjust = 0)
pareto120NOlarge


#### WRITE PLOTS TO MULTI-PANEL PNG FIGURE.----

## Numeric labels and arrows of individual Pareto frontier points done manually outside R to prevent overlap
## Use previous legend from manuscript figure and manually add to this as well

png(file="Simulations/Figures/FigureParetoNoCost.png",width=12,height=9,units="in",res=600)
((pareto60NOsmall + pareto120NOsmall) / (pareto60NOlarge + pareto120NOlarge))
dev.off()



#### APPENDIX - FULL STARTUP.----

## Scenario: Full start-up costs (aka, all equipment needed) plus recurring costs
## Root Mean Square Error (RMSE) - square root of the mean of the square of all error, lower value indicates more accuracy


## PANEL ONE. LOW DENSITY (60), SMALL SNAKES
lowdenssmall <- subset(fullcost, Truth == 60 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to full start up and divide by 1000 to put in terms of thousands of USD
lowdenssmall$TotalCost <- (lowdenssmall$TotalCost)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenssmall$TotalCost, y = lowdenssmall$RMSE)
#Subset to names of Sim on frontier
pflowsmall <- subset(fullcost, Truth == 60 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("small") & Sim %in% pflowsmall[[1]]) %>%
  arrange((TotalCost)/1000) %>%
  mutate(Num = 1:length(pflowsmall[[1]]))

pareto60FULLsmall <- ggplot(subset(fullcost, Truth == 60 & Size == c("small")), aes(x = (TotalCost)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.6)) +
  ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 10), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,85) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 64, xmax = 85, ymin = 160, ymax = 320, alpha = 0.2) +
  annotate("text", x = 65, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 65, y = 290, label = "1. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 250, label = "3. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 230, label = "4. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 210, label = "5. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 190, label = "6. closedVIS30full", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 170, label = "7. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60FULLsmall


## PANEL TWO. NORMAL DENSITY (120), SMALL SNAKES
normdenssmall <- subset(fullcost, Truth == 120 & Size == c("small"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to full start up and divide by 1000 to put in terms of thousands of USD
normdenssmall$TotalCost <- (normdenssmall$TotalCost)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenssmall$TotalCost, y = normdenssmall$RMSE)
#Subset to names of Sim on frontier
pfnormsmall <- subset(fullcost, Truth == 120 & Size == c("small"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("small") & Sim %in% pfnormsmall[[1]]) %>%
  arrange((TotalCost)/1000) %>%
  mutate(Num = 1:length(pfnormsmall[[1]]))

pareto120FULLsmall <- ggplot(subset(fullcost, Truth == 120 & Size == c("small")), aes(x = (TotalCost)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.6)) +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority small snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,85) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 64, xmax = 85, ymin = 200, ymax = 320, alpha = 0.2) +
  annotate("text", x = 65, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 65, y = 290, label = "1. onewayVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 250, label = "3. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 230, label = "4. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 210, label = "5. closedVIS30full", size = 2.5, hjust = 0)
pareto120FULLsmall


## PANEL THREE. LOW DENSITY (60), LARGE SNAKES
lowdenslarge <- subset(fullcost, Truth == 60 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to full start up and divide by 1000 to put in terms of thousands of USD
lowdenslarge$TotalCost <- (lowdenslarge$TotalCost)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = lowdenslarge$TotalCost, y = lowdenslarge$RMSE)
#Subset to names of Sim on frontier
pflowlarge <- subset(fullcost, Truth == 60 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 60 & Size == c("large") & Sim %in% pflowlarge[[1]]) %>%
  arrange((TotalCost)/1000) %>%
  mutate(Num = 1:length(pflowlarge[[1]]))

pareto60FULLlarge <- ggplot(subset(fullcost, Truth == 60 & Size == c("large")), aes(x = (TotalCost)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.6)) +
  xlab("Cost in Thousands of USD") + ylab("Root Mean Square Error") +
  theme(legend.position = "none", axis.title = element_text(size = 12), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Low density (60) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,85) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 64, xmax = 85, ymin = 180, ymax = 320, alpha = 0.2) +
  annotate("text", x = 65, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 65, y = 290, label = "1. closedVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 250, label = "3. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 230, label = "4. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 210, label = "5. closedVIS30full", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 190, label = "6. closedVISTRAP30halfhalf", size = 2.5, hjust = 0)
pareto60FULLlarge


## PANEL FOUR. NORMAL DENSITY (120), LARGE SNAKES
normdenslarge <- subset(fullcost, Truth == 120 & Size == c("large"))[,c("TotalCost","Cost.x","RMSE","Sim")]
#Correct total cost to full start up and divide by 1000 to put in terms of thousands of USD
normdenslarge$TotalCost <- (normdenslarge$TotalCost)/1000
#Get the row number of all Pareto efficient actions
numSim <- getPareto(x = normdenslarge$TotalCost, y = normdenslarge$RMSE)
#Subset to names of Sim on frontier
pfnormlarge <- subset(fullcost, Truth == 120 & Size == c("large"))[numSim[!is.na(numSim)],c("Sim")]

## Flag points in dataframe on frontier in order to label
highlight_df <- fullcost %>%
  filter(Truth == 120 & Size == c("large") & Sim %in% pfnormlarge[[1]]) %>%
  arrange((TotalCost)/1000) %>%
  mutate(Num = 1:length(pfnormlarge[[1]]))

pareto120FULLlarge <- ggplot(subset(fullcost, Truth == 120 & Size == c("large")), aes(x = (TotalCost)/1000, y = RMSE, fill = TypeALL, shape = Space, size = as.factor(Days))) +
  geom_point(alpha = 0.7, position = position_dodge(width = 0.6)) +
  xlab("Cost in Thousands of USD") +
  theme(legend.position = "none", axis.title.x = element_text(size = 12), axis.text.x = element_text(size = 10), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values = c("#567022","#B77A29","#0066cc","#C9F277","#FDEC9E","#00ccff"), breaks = c("closedTRAP","closedVIS","closedVISTRAP","onewayTRAP","onewayVIS","onewayVISTRAP")) +
  scale_shape_manual(values = c(21, 22, 24, 25)) +
  scale_size_manual(values=c(3,5)) +
  ggtitle("Normal density (120) and majority large snakes") +
  ylim(0,(max(fullcost$RMSE)+5)) +
  xlim(15,85) +
  ## Pareto efficient actions identified above
  annotate("rect", xmin = 64, xmax = 85, ymin = 180, ymax = 320, alpha = 0.2) +
  annotate("text", x = 65, y = 310, label = "Pareto efficient actions", size = 3.5, hjust = 0) +
  annotate("text", x = 65, y = 290, label = "1. onewayVIS14half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 270, label = "2. closedVISTRAP14thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 250, label = "3. closedVIS14full", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 230, label = "4. closedVIS30half", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 210, label = "5. closedVISTRAP30thirdthird", size = 2.5, hjust = 0) +
  annotate("text", x = 65, y = 190, label = "6. closedVIS30full", size = 2.5, hjust = 0)
pareto120FULLlarge


#### WRITE PLOTS TO MULTI-PANEL PNG FIGURE.----

## Numeric labels and arrows of individual Pareto frontier points done manually outside R to prevent overlap
## Use previous legend from manuscript figure and manually add to this as well

png(file="Simulations/Figures/FigureParetoFullCost.png",width=12,height=9,units="in",res=600)
((pareto60FULLsmall + pareto120FULLsmall) / (pareto60FULLlarge + pareto120FULLlarge))
dev.off()
