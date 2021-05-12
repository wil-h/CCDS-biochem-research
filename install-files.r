library(bio3d)

dir = getwd()

path = readline("path of output directory: ")
filename = readline("enter file containing pdb accession codes: ")

pdb_ids = scan(filename,what='character',sep=",")
get.pdb(pdb_ids,path=paste(dir, "/", path, sep=''))
