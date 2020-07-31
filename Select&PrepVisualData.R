#### July 29, 2020; S. Amburgey ####
## Read in data (survey info and captures) from visual surveys from numerous projects
## Select what projects to use
## Format for SCR analysis

allsurv <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/VISinfo.csv")[,-c(21:62)]

## Study areas deemed suitable
subsurv <- subset(allsurv, SITE %in% c("NWFN","NWFR","HMUI","HMUR","HMUI1","HMUI2","HMUI2B","HMUI3","HMUI4","HMUI5","HMUI5B","NCRI","NCRR"))

## Projects deemed suitable (marked snakes, enough recaptures, closed period of time, spatial information, etc.)
subsurv <- subset(subsurv, PROJECTCODE %in% c("NWFO VIS 1","NWFO VIS 2","NWFO VIS HL 1","NWFO VIS HL 2","NWFN VIS 1","NWFO VIS 2","NWFN VIS 2","NWFN VISPACE", "NWFN VISTRAP VIS", "POST KB VIS", "POST KB VIS 1", "POST KB VIS 2",  "POST KB VIS 3", "POST KB VIS 3 EXTRA", "PRE BT2 VIS", "EDGE EFFECT VIS", "LOWDENS VIS"))

