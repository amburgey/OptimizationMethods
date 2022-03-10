### Check rhats on saved results to identify simulations that need to be re-run

library(jagsUI)

#create a list of the files from your target directory
file_list <- list.files(path="D:/Results")


for(i in 1:length(file_list)){
    
  load(file=paste("D:/Results/",file_list[i],sep=""))
    
  print(i)
    
  ifelse(any(unlist(out$Rhat) > 1.1) == TRUE, stop("did not converge"), 1)
  
}
  
