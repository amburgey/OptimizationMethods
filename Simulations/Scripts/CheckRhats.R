### Script to check rhats of saved results to identify simulations that need to be re-run
### list.files will run through all files in the specified directory
### Script will stop and throw error when it finds a failed simulation - have to manually restart after removing the failed simulation

library(jagsUI)

#create a list of the files from your target directory
file_list <- list.files(path="D:/Simulations/Results")


for(i in 1:length(file_list)){
  
  load(file=paste("D:/Simulations/Results/",file_list[i],sep=""))
    
  print(i)
    
  ifelse(any(unlist(out$Rhat) > 1.1095) == TRUE, stop("did not converge"), 1)
  
}
  
