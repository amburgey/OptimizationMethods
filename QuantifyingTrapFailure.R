#### Summary of trap failure throughout all trap projects used in optimization work ####
## Trap failure = denotes a trap was functionally unable to or imperfect in its ability to capture a snake on a given evening

rm(list=ls())

source("Select&PrepTrapData.R")  ## Creation of subcap and subsurv

## Master database of nonfunctional traps but this just lists specific traps that needed maintenance (not necessarily all traps that were nonfunctional, see mismatch in number of traps being run when compared to survey data and trap counts)
fail <- read.csv("/Users/Staci Amburgey/Documents/USGS/BrownTreesnakes/Optim Monitoring Methods/Data/TrapNonfunctional.csv")

## Surveys where some specific note was made
commented <- subset(subsurv, COMMENTS != "")
## 274 traps have some kind of comment on them
## Some traps were not secured well and COULD HAVE resulted in a snake escaping, but some of these still catch snakes
## Failure often that door was stuck open, traps not properly clasped, mouse chamber blocking entry door, mouse missing or dead, non-functional (unspecified), or pig attacked trap

trapnums <- table(subsurv$PROJECTCODE, subsurv$COUNT)
## No trap failure for SWIFT TRAP ARRAY sampling or not included in the master database either in TRAP file or nonfunctional file

trapfail <- subsurv
trapfail$trapgoal <- ifelse(trapfail$PROJECTCODE == "GNWR SSP MOUSE", 18, 
                            ifelse(trapfail$PROJECTCODE == "NWFN TRAP 1", 13,
                            ifelse(trapfail$PROJECTCODE == "NWFN TRAP 2 LINVIS", 13,
                            ifelse(trapfail$PROJECTCODE == "NWFN TRAP 3", 13,
                            ifelse(trapfail$PROJECTCODE == "NWFN TRAP 4 LCM", 13,
                            ifelse(trapfail$PROJECTCODE == "NWFN VISTRAP TRAP", 13, 
                            ifelse(trapfail$PROJECTCODE == "POST BT2 TRAP", 13,
                            ifelse(trapfail$PROJECTCODE == "POST KB TRAP 1", 13,
                            ifelse(trapfail$PROJECTCODE == "POST KB TRAP 2", 13,
                            ifelse(trapfail$PROJECTCODE == "PRE BT1 TRAP", 13,
                            ifelse(trapfail$PROJECTCODE == "SWIFT TRAP ARRAY", NA, 
                            ifelse(trapfail$PROJECTCODE == "TOX DROP TRAP", 10, -99999))))))))))))

trapfail$fail <- ((trapfail$trapgoal - trapfail$COUNT)/trapfail$trapgoal)*100

max(trapfail$fail, na.rm = TRUE)
min(trapfail$fail, na.rm = TRUE)

mean(trapfail$fail, na.rm = TRUE)
median(trapfail$fail, na.rm = TRUE)
dim(subset(trapfail, fail > 0))[1]  # transects with any failure

