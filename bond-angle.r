#setwd("~/")
#.libPaths(c()) place library paths here.

library(MASS)
library(hash)
library(devtools) 
library(bio3d) 


get_atom <- function (current_resno, starting_index, atom.data, elety)
{
  i = starting_index
  if (substring(elety,nchar(elety)) == "-") {
    i = i - 1
    if (i == 0) {return ("N/A")}
    while(as.numeric(current_resno) - 1 == as.numeric(atom.data[i,3])) {
      if (atom.data[i,4] == substring (elety,1,nchar (elety)-1)) {return (atom.data[i,5:7])}
      i = i-1
      if (i == 0) {return ("N/A")}
    }
    return ("N/A")
  }
  else if (substring(elety,nchar(elety)) == "+")
  {
    while(current_resno == atom.data[i,3]) {
      i = i + 1
      if (i == length(atom.data)/7 + 1) {return ("N/A")}
    }
    while (as.numeric(current_resno) + 1 == as.numeric(atom.data[i,3]))
    {
      if (atom.data[i,4] == substring (elety,1,nchar (elety)-1)) {return (atom.data[i,5:7])}
      i = i+1
      if (i == length(atom.data)/7 + 1) {return ("N/A")}
    }
  }
  else {
    while(current_resno == atom.data[i,3]) 
    {
      if (atom.data[i,4] == elety && atom.data[i,3] == current_resno) {return (atom.data[i,5:7])}
      i = i + 1
      if (i == length(atom.data)/7 + 1) {return ("N/A")}
    }
  }
  return ("N/A")
}


dir = getwd() 

resid <- readline(prompt="Enter 3 letter residue ID: ")
path = readline("path of output directory: ")
filename = readline("enter file containing pdb accession codes (leave blank if files are already downloaded): ")

if (filename != "") {
  pdb_ids = scan(filename,what='character',sep=",")
  get.pdb(pdb_ids,path=paste(dir, "/", path, sep=''))
}

print(paste(dir, "/", path, sep=''))
setwd(paste(dir, "/", path, sep=''))

colpaste <- function(x, col.names = colnames(x)) {
  apply(x, 1, function(row) paste(row[col.names], collapse = ","))
}

#prompts the user to type in the console a three letter resid ID, such a PRO for proline


pdbfiles <- list.files(pattern="*.pdb", full.names=TRUE)

result_matrix = matrix(data = c("pdb code","chain","residue number", "CO-N-CA","N-CA-CO","N-CA-CB","CA-CB-CG","CB-CG-CDa","CB-CG-CDb","CB-CA-CO","CA-CO-N"), nrow = 1, ncol = 11)
atoms = c ("C-","C", "CA","N","CB","CD1","CD2", "CG","N+")
bonds = c ("C-","N","CA","N","CA","C","N","CA","CB","CA","CB","CG","CB","CG","CD1","CB","CG","CD2","CB","CA","C","CA","C","N+")

for(current_file in pdbfiles){
  
  pdb <- read.pdb(current_file)
  index <- trim(pdb, resid = resid, inds = NULL, sse = TRUE)
  index <- index[lapply(index,length)>0]
  atom.data <- colpaste(pdb$atom, c("resid","chain","resno","elety","x", "y","z"))
  atom.data <- matrix(unlist(strsplit(atom.data, "\\,")),ncol = 7, byrow = TRUE)
  
  #print(atom.data)
  i = 1
  
  while(resid != atom.data[i,1]) {
    i = i + 1
    if(i == length(atom.data)/7) {
      break #needed here in the unlikely event that a file has zero of the target
      #residue present.
    }
  }
  
  while (i < length(atom.data)/7) {
    current_resno = atom.data[i,3]
    
    coords <- hash()
    for (atom in atoms) {
      coords[[atom]] = get_atom (current_resno, i, atom.data, atom)
    }
    
    angle_vector = c(current_file, atom.data[i,2],current_resno)
    for (bond_index in 1:(length(bonds)/3))
    {
      angle = angle.xyz(c(as.numeric(get(bonds[bond_index*3-2], coords)),as.numeric(get(bonds[bond_index*3-1], coords)),as.numeric(get(bonds[bond_index*3], coords))))
      angle_vector = c(angle_vector,angle)
    }
    result_matrix = rbind(result_matrix, angle_vector)
    repeat { #since r does not have a "do while" loop, we need to use this instead
      i = i + 1
      if(i == length(atom.data)/7 || (current_resno != atom.data[i,3] && resid == atom.data[i,1])) {
        break
      }
    }
  }
}

write.matrix(result_matrix, file = paste("angle_",resid,".csv",sep=""),sep=",") 
