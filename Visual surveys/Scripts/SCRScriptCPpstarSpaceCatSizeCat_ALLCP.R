#### Unified analysis of all CP (NWFN) datasets ####

rm(list=ls())

library(secr); library(reshape2); library(jagsUI)

projects <- c("Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VIS2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISHL1.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISHL2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPREBT2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTBT2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB1.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB2.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISPOSTKB3.R",
        "Visual surveys/Scripts/SCRScriptCPpstarSpaceCatSizeCat_VISVISTRAP.R")

for(i in 1:length(projects)){
        ## Read in prep code for project
        source(projects[1])
        ## Create and add to a combined dataset
        
        
}

