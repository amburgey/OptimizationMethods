### Plotting snake reproductive pulses ###

library(tidyverse); library(lubridate); library(ggplot2)

data <- read_csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VIScaptures.csv")

alldata <- data %>%
  dplyr::select(Date,SITE,PITTAG,NEW,SVL,SEX,TRANSECT,LOCATION,TRAPTYPE) %>%
  mutate_if(is.character, str_replace_all, pattern = "HMUI1", replacement = "HMUI") %>%
  mutate_if(is.character, str_replace_all, pattern = "HMUI2", replacement = "HMUI") %>%
  mutate_if(is.character, str_replace_all, pattern = "HMUI3", replacement = "HMUI") %>%
  mutate_if(is.character, str_replace_all, pattern = "HMUI4", replacement = "HMUI") %>%
  mutate_if(is.character, str_replace_all, pattern = "HMUI5", replacement = "HMUI") %>%
  filter(TRAPTYPE == "V", SITE == c("NWFN","HMUR","HMUI","ASOI","NCRI","NCRR")) %>%  #, SVL < 700
  mutate(Date = dmy(Date), Year = year(Date), Month = month(Date), Time = paste(Year,Month,sep=".")) %>%
  filter(Year < 2017) %>%
  mutate(Time = fct_relevel(Time, paste(rep(2004:2016,each=12),".", rep(1:12,each=1), sep="")))

puls1 <- ggplot(alldata, aes(x=Time, y=SVL)) + geom_point() +
  facet_grid(rows = vars(SITE), scales = "free_x") +
  geom_hline(yintercept = 700, lty = 2, col = "red")

test <- expand.grid(Year = unique(alldata$Year), Month = unique(alldata$Month))
CPdata <- alldata %>%
  filter(SITE == "NWFN" | SITE == "NA")

CPdata <- merge(CPdata, test, by = c("Year","Month"), all = TRUE)
CPdata <- CPdata[order(CPdata$Year, CPdata$Month),]
CPdata$Group <- ifelse(is.na(CPdata$Time),"None","Survey")
CPdata$Time2 <- NA
for(i in 1:nrow(CPdata)){
  CPdata[i,14] <- paste(CPdata[i,1],".",CPdata[i,2],sep="")
}

CPdata$SVL[is.na(CPdata$SVL)]<-0

pulsNWFN <- ggplot(CPdata, aes(x=Time2, y=SVL, color=factor(Month), shape=Group)) + 
  geom_point() +
  geom_hline(yintercept = 700, lty = 2, col = "red") + 
  theme(axis.text.x = element_text(angle = 90))
