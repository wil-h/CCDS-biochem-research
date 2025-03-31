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

closest_integer <- function(item, chain, resnos){
  closest=10000000
  k=1
  while(k<=length(resnos)){
    if(abs(as.integer(item)-closest)>abs(as.integer(item)-as.integer(resnos[k])) && chain==resnos[k+1]){
      closest=as.integer(resnos[k])
    }
    k=k+2
  }
  return(closest)
}

#find the valid resnos to use the filler data
valid_pdb_resnos <- function(atomdata,resid){
  i=1
  while(i<=nrow(atomdata)){
    if(atomdata[i,1]!=resid){
      i=i+1
    }
    else{
      break
    }
  }
  valid_resnos <- c()
  while(i<=nrow(atomdata)){
    value <- atomdata[i,3]
    chain <- atomdata[i,2]
    if(atomdata[i,1]==resid){
      prevvalue=valid_resnos[length(valid_resnos)-1]
      if(is.null(prevvalue)){
        prevvalue=0
      }
      if (value != "torsion vector" && prevvalue != value) {
        valid_resnos <- c(valid_resnos, value)  
        valid_resnos <- c(valid_resnos, chain) 
      }
    }
    i=i+1
  }
  return(valid_resnos)
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
filler_torsion_matrix = matrix(data = torsion_header, nrow = 1, ncol = length(torsion_header))
torsion_row2 = config$torsion$header2
torsion_matrix = rbind(torsion_matrix, torsion_row2)
filler_torsion_matrix = rbind(filler_torsion_matrix, torsion_row2)
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
  tryCatch(expr={ 
    current_file = pdbfiles[file_index]
    if (file_index %% 5 == 0) {print(paste("thread ",args[1]," has processed [", file_index, "/",length(pdbfiles), "] files",sep = ""))}


    pdb <- read.pdb(current_file, rm.alt=FALSE)
    index <- trim(pdb, resid = resid, inds = NULL, sse = TRUE)
    index <- index[lapply(index,length)>0]
    atom.data <- colpaste(pdb$atom, c("resid","chain","resno","elety","x", "y","z","o","alt"))
    atom.data <- matrix(unlist(strsplit(atom.data, "\\,")),ncol = 9, byrow = TRUE)
    atom.data.torsion <- matrix(ncol=9,byrow=TRUE)

    #ensure there are no preceding spaces in the atom data resnos
    for(i in 1:(length(atom.data)/9)){
      atom.data[i,3]=gsub(" ","",atom.data[i,3])
    }

    #find split branches and reformat them for torsion
    subnum=0
    atom_index=1
    while(atom_index<=(length(atom.data)/9)){
      if(atom.data[atom_index,1]==resid && atom.data[atom_index,9]!="NA"){
        affectedResno=atom.data[atom_index,3]
        affectedStrand=atom.data[atom_index,2]
        resAtoms=matrix(ncol=9,byrow=TRUE)
        for(l in 1:(length(atom.data)/9)){
          if(atom.data[l,3]==affectedResno && atom.data[l,2]==affectedStrand){
            resAtoms=rbind(resAtoms,atom.data[l,])
          }
        }
        resAtoms=resAtoms[-1,]
        aRes=matrix(ncol=9,byrow=TRUE)
        bRes=matrix(ncol=9,byrow=TRUE)
        for(l in 1:(length(resAtoms)/9)){
          if(resAtoms[l,9]=="NA"){
            aRes=rbind(aRes,resAtoms[l,])
            bRes=rbind(bRes,resAtoms[l,])
          }
          if(resAtoms[l,9]=="A"){
            aRes=rbind(aRes,resAtoms[l,])
          }
          if(resAtoms[l,9]=="B"){
            bRes=rbind(bRes,resAtoms[l,])
          }
        }
        aRes=aRes[-1,]
        bRes=bRes[-1,]
        for (i in 1:nrow(bRes)){
          bRes[i, 3]=as.character(as.integer(bRes[i, 3]) + 1) 
        }

        atom.data.torsion = rbind(atom.data.torsion, unname(aRes))
        atom.data.torsion = rbind(atom.data.torsion, unname(bRes))
        subnum=subnum+1
        while(atom.data[atom_index,3]==affectedResno){
          atom_index=atom_index+1
        }

      }
      else{
        row=atom.data[atom_index,]
        row[3]=as.character(as.integer(row[3])+subnum)
        atom.data.torsion = rbind(atom.data.torsion, unname(row))
      }
      atom_index=atom_index+1
    }
    atom.data.torsion=atom.data.torsion[-1,]#remove the blank first row

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
        torsion_vector = c(torsion_vector,"OCC")#adds an occupancy placeholder
      }
      filler_torsion_matrix = rbind(filler_torsion_matrix, torsion_vector)
      
      repeat { #since r does not have a "do while" loop, we need to use this instead
        i = i + 1
        if(i == length(atom.data)/9 || (current_resno != atom.data[i,3] && resid == atom.data[i,1])) {
          break
        }
      }
    }
    #repeat for torsion with different atom data
      i = 1
    while(resid != atom.data[i,1]) {
      i = i + 1
      if(i == length(atom.data)/9) {
        break #needed here in the unlikely event that a file has zero of the target
        #residue present.
      }
    }

    while (i < length(atom.data.torsion)/9) {
      current_resno = atom.data.torsion[i,3]
      coords <- hash()
      for (atom in atoms) {
        coords[[atom]] = get_atom (current_resno, i, atom.data.torsion, atom)
      }
      torsion_vector = c(current_file, atom.data.torsion[i,2],current_resno)

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
        bondOccupancy= atom.data.torsion[i+bond_index,8]
        torsion_vector = c(torsion_vector,bondOccupancy)
      }
      torsion_matrix = rbind(torsion_matrix, torsion_vector)

      #now go through with the filler data and fill in any blanks that were caused by the modified atom data
      p=1
      while(p<nrow(torsion_matrix)){
        if(torsion_matrix[p,1]==current_file){
          break
        }
        else{
          p=p+1
        }
      }
      valid_resnos=valid_pdb_resnos(atom.data,resid)
      for(row in p:nrow(torsion_matrix)){
        for(column in 4:ncol(torsion_matrix)){
          if(torsion_matrix[row,column]==0){
            closest = closest_integer(torsion_matrix[row,3],torsion_matrix[row,2],valid_resnos)
            #closest working, now just need to change the filling logic below
            k=1
            while(k<=nrow(filler_torsion_matrix)){
              if(filler_torsion_matrix[k,1]==current_file && as.integer(filler_torsion_matrix[k,3])==as.integer(closest)){
                break
              }
              else{
                k=k+1
              }
            }
            torsion_matrix[row,column]=filler_torsion_matrix[k,column]
          }
        }
      } 
      
      repeat { #since r does not have a "do while" loop, we need to use this instead
        i = i + 1
        if(i == length(atom.data.torsion)/9 || (current_resno != atom.data.torsion[i,3] && resid == atom.data.torsion[i,1])) {
          break
        }
      }
    }
  }, error=function(e){
    print(paste("Error processing file: ", current_file, " in thread ", args[1], sep=""))})
}

write.matrix(angle_matrix, file = paste("angle_",args[1],".csv",sep=""),sep=",") 
write.matrix(length_matrix, file = paste("length_",args[1],".csv",sep=""),sep=",") 
write.matrix(torsion_matrix, file = paste("torsion_",args[1],".csv",sep=""),sep=",") #there are some issues with the torsion output, but they are cleaned up in collect_data.py