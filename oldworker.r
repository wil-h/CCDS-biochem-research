#this is very similar to the other worker, but it does not account for splitting with occupancy values.
#if there are two different conformations from a single pdb, it will just choose one to put in the spreadsheet. 
#if your AA does not have any splitting, then it is okay to use this version.
options(warn=-1)
args = commandArgs(trailingOnly=TRUE)
.libPaths(c(args[3],args[4]))


library(MASS)
library(hash)
library(devtools)
library(bio3d)
library(rjson)

args = commandArgs(trailingOnly=TRUE)

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
      if (i == length(atom.data)/9 + 1) {return ("N/A")}
    }
    while (as.numeric(current_resno) + 1 == as.numeric(atom.data[i,3]))
    {
      if (atom.data[i,4] == substring (elety,1,nchar (elety)-1)) {return (atom.data[i,5:7])}
      i = i+1
      if (i == length(atom.data)/9 + 1) {return ("N/A")}
    }
  }
  else {
    while(current_resno == atom.data[i,3]) 
    {
      if (atom.data[i,4] == elety && atom.data[i,3] == current_resno) {return (atom.data[i,5:7])}
      i = i + 1
      if (i == length(atom.data)/9 + 1) {return ("N/A")}
    }
  }
  return ("N/A")
}

colpaste <- function(x, col.names = colnames(x)) {
  apply(x, 1, function(row) paste(row[col.names], collapse = ","))
}

calc_bondlength <- function (array1, array2)
{
  
  x1 <- c(as.numeric(array1[1]))
  x2 <- c(as.numeric(array2[1]))
  y1 <- c(as.numeric(array1[2]))
  y2 <- c(as.numeric(array2[2]))
  z1 <- c(as.numeric(array1[3]))
  z2 <- c(as.numeric(array2[3]))
  length = c(sqrt(sum((((x2-x1))^2), (((y2-y1))^2), (((z2-z1))^2))))
  
  return (length) 
}

############################
#boring initializaton stuff#
############################

config <- fromJSON(file = args[2])

atoms = config$atoms
resid = config$resid

angle_header = config$angle$header
angle_matrix = matrix(data =angle_header, nrow = 1, ncol = length(angle_header))
angle_bonds = config$angle$bonds
angle_row2 = config$angle$header2
angle_matrix = rbind(angle_matrix, angle_row2)

torsion_header = config$torsion$header
torsion_matrix = matrix(data = torsion_header, nrow = 1, ncol = length(torsion_header))
torsion_row2 = config$torsion$header2
torsion_matrix = rbind(torsion_matrix, torsion_row2)
torsion_bonds = config$torsion$bonds

length_header = config$length$header
length_matrix = matrix(data = length_header, nrow = 1, ncol = length(length_header))
length_row2 = config$length$header2
length_matrix = rbind(length_matrix, length_row2)
length_bonds = config$length$bonds

filename = paste("worker_", args[1] ,".txt", sep='')
pdbfiles = scan(filename,what='character',sep=",")

dir = getwd() 
setwd(paste(dir, "/", config$output_dir, sep=''))

#ok, that's done, now for the actual program.
for(file_index in 1:length(pdbfiles)){
  current_file = pdbfiles[file_index]
  if (file_index %% 5 == 0) {print(paste("thread ",args[1]," has processed [", file_index, "/",length(pdbfiles), "] files",sep = ""))}


  pdb <- read.pdb(current_file, rm.alt=FALSE)
  index <- trim(pdb, resid = resid, inds = NULL, sse = TRUE)
  index <- index[lapply(index,length)>0]
  atom.data <- colpaste(pdb$atom, c("resid","chain","resno","elety","x", "y","z","o","alt"))
  atom.data <- matrix(unlist(strsplit(atom.data, "\\,")),ncol = 9, byrow = TRUE)
  
  print(atom.data)

  i = 1
  
  while(resid != atom.data[i,1]) {
    i = i + 1
    if(i == length(atom.data)/9) {
      break #needed here in the unlikely event that a file has zero of the target
      #residue present.
    }
  }
  
  while (i < length(atom.data)/9) {
    current_resno = atom.data[i,3]
    coords <- hash()
    for (atom in atoms) {
      coords[[atom]] = get_atom (current_resno, i, atom.data, atom)
    }
    angle_vector = c(current_file, atom.data[i,2],current_resno)
    torsion_vector = c(current_file, atom.data[i,2],current_resno)
    length_vector = c(current_file, atom.data[i,2],current_resno)
    
    for (bond_index in 1:(length(angle_bonds)/3))
    {
      tryCatch (
        expr={
          angle_res = angle.xyz(c(as.numeric(get(angle_bonds[bond_index*3-2], coords)),as.numeric(get(angle_bonds[bond_index*3-1], coords)),as.numeric(get(angle_bonds[bond_index*3], coords))))
          },
        error=function (e){angle_res="NA"})
      angle_vector = c(angle_vector,angle_res)
    }
    angle_matrix = rbind(angle_matrix, angle_vector)
    
    for (bond_index in 1:(length(length_bonds)/2))
    {
      tryCatch(
        expr = {
          array1= get(length_bonds[bond_index * 2 - 1],coords)
          array2= get(length_bonds[bond_index * 2], coords)
          len = calc_bondlength(array1,array2)},
        error=function (e){len="NA"})
      length_vector = c(length_vector,len)
    }
    length_matrix = rbind(length_matrix, length_vector)

    for (bond_index in 1:(length(torsion_bonds)/4))
    {
      tryCatch(
        expr = {
          array1 = get(torsion_bonds[bond_index * 4 - 3],coords)
          array2= get(torsion_bonds[bond_index * 4 - 2], coords)
          array3= get(torsion_bonds[bond_index * 4 - 1], coords)
          array4= get(torsion_bonds[bond_index * 4], coords)
          torsional_angle = torsion.xyz(c(as.numeric(array1),as.numeric(array2),as.numeric(array3),as.numeric(array4)))

        },
        error=function (e){torsional_angle="NA"})
      torsion_vector = c(torsion_vector,torsional_angle)#adds an angle each time
      bondOccupancy= atom.data[i,8]
      torsion_vector = c(torsion_vector,bondOccupancy)
      #print(c("angle:",torsional_angle))
      #print(c("occupancy:",bondOccupancy))
      #print(torsion_vector)
    }
    torsion_matrix = rbind(torsion_matrix, torsion_vector)
    
    repeat { #since r does not have a "do while" loop, we need to use this instead
      i = i + 1
      if(i == length(atom.data)/9 || (current_resno != atom.data[i,3] && resid == atom.data[i,1])) {
        break
      }
    }
  }
}


write.matrix(angle_matrix, file = paste("angle_",args[1],".csv",sep=""),sep=",") 
write.matrix(length_matrix, file = paste("length_",args[1],".csv",sep=""),sep=",") 
write.matrix(torsion_matrix, file = paste("torsion_",args[1],".csv",sep=""),sep=",") 
