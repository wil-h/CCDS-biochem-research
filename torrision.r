setwd("~/") #moving back into the original directory
#for debug purposes.

#install dependancies.
#TODO, check if packages are already installed.
install.packages("bio3d", dependencies=TRUE) 
install.packages("devtools") 

#import dependancies
library(devtools) 
library(bio3d) 
#install_bitbucket("Grantlab/bio3d", subdir = "ver_devel/bio3d/")  




dir = getwd() 
path = readline("enter relative path: ") #TODO support absolute and relative paths

print(paste(dir, "/", path, sep=''))
setwd(paste(dir, "/", path, sep=''))



#prompts the user to type in the console a three letter resid ID, such a PRO for proline
resid <- readline(prompt="Enter 3 letter residue ID: ")

pdbfiles <- list.files(pattern="*.pdb", full.names=TRUE) # Rerieving .pdb files 

#for loop that does the function for each pdb in pdbfiles vector
for (i in pdbfiles){
  #uses the read.pdb function in bio3d to store the pdb file from the working directory as a usable vector matrix with all the atoms and coordinates
  pdb <- read.pdb(i)
  
  #removes all data except date for the residue that was typed into console
  index <- trim(pdb, resid = resid, inds = NULL, sse = TRUE)
  
  #not quite sure what this does but there's errors otherwise... i think it removes null values
  index <- index[lapply(index,length)>0]
  
  #prints the pdb file name
  print(i)
  
  #prints the table summary of the output of torsion.pdb function, the 'try' part means it only prints if the function runs without error
  try(print(torsion.pdb(index)$tbl))
  excel_data = torsion.pdb(index)$tbl
  y <- as.character(excel_data) 
  write.csv(excel_data, file = 'tor.csv', row.names = F) 
}
