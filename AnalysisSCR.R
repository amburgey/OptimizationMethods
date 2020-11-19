##### Begin analysis small - start with one site and build from there #####
## Create dataframes of captures (PITTAG by Date) and activity (TRANSECT & LOCATION by Date)

source("Select&PrepVisualData.R")
source("SpecifyStateSpace.R")

library(tibble)


## Start with HMU as that state space is a little tricky
## HMU EDGE EFFECT
caps <- tibble(subcap)
caps$seen <- c(1)
survs <- tibble(subsurv)

## Snake dataframe (PITTAG by Date)
hmuEsnks <- as.data.frame(caps %>% 
  filter(SITE == c("HMUI", "HMUR") & PROJECTCODE == c("EDGE EFFECT VIS")) %>%
  dplyr::select(PITTAG,Date,TRANSECT, LOCATION, seen) %>% 
  pivot_wider(names_from = Date, values_from = seen, values_fill = 0))  ## add traps

## how to handle capture info if not gridded? Overlay on grid and assign? How to assign trap location if just coordinates/meter markers?

## Active/inactive dataframe
hmuEsurvs <- as.data.frame(survs %>% 
  filter(SITE == c("HMUI", "HMUR") & PROJECTCODE == c("EDGE EFFECT VIS")) %>%
  dplyr::select(TRANSECT,Date,DISTANCE) %>% 
  pivot_wider(names_from = Date, values_from = seen, values_fill = 0))


SX[i] ~ dunif(0, xu) #priors for the activity centre for each individual
SY[i] ~ dunif(0, yu) #  xu, yu are the upper limits of x, y
pOK[i] <- habmat[trunc(SX[i]+1), trunc(SY[i]+1)] # habitat check
OK[i] ~ dbern(pOK[i]) # OK[i] = 1, the ones trick