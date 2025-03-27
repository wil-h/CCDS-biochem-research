from subprocess import Popen
import multiprocessing
import json
import time
import sys
import csv
import os

start_time = time.time()

if (len(sys.argv) == 1):
    print("please specify config file")
    exit()

config_data = json.loads(open(sys.argv[1]).read())
output_dir = config_data['output_dir']

if  (config_data['cores'] > 0):
    cores = int(config_data['cores'])

elif (config_data['cores'] < 0):
    cores = multiprocessing.cpu_count() - int(config_data['cores'])
    
else:
    cores = multiprocessing.cpu_count()

cw_dir = os.getcwd().replace("\\", "/")
out_dir = cw_dir + "/" + config_data["output_dir"]

os.chdir(out_dir)
try:
    os.remove("bond_length" + config_data['resid'] + ".csv")
except:
    None
try:
    os.remove("bond_angle" + config_data['resid'] + ".csv")
except:
    None
try:
    os.remove("torsion" + config_data['resid'] + ".csv")
except:
    None

files = os.listdir()
pdbs = []
for filename in files:
    if (filename[-4:] == ".pdb"):
        pdbs += [filename]

tasks = [[] for i in range(cores)]
for file_index in range(len(pdbs)): tasks[file_index % cores] += [pdbs[file_index]]

os.chdir(cw_dir)

for array_num in range(cores):
    pdb_chunk = ','.join(tasks[array_num])
    task_file = open("worker_" + str(array_num + 1) + ".txt", "w")
    task_file.write(pdb_chunk)
    task_file.close()

processes = [Popen('"' + config_data["rscript_path"] + '" --vanilla worker.r '+ str(i + 1) + ' ' + sys.argv[1] + ' ' + config_data["libpaths"][0]+ ' ' + config_data["libpaths"][1], shell=False) for i in range(cores)]
exitcodes = [p.wait() for p in processes]   

os.chdir(out_dir)
files = os.listdir()
csv_delimeters = ["a","l","t"]
csvs = [[],[],[]]
matrices = [[],[],[]]

#clean weird first row in some csvs:
for filepath in files:
    if filepath[-4:] == ".csv" and "torsion" in filepath:
        with open(filepath, mode='r', newline="") as file:
            reader=list(csv.reader(file))
            firstval=reader[0][0]
        if "pdb code" not in firstval:
            modified = reader[1:] 
            with open(filepath, mode='w', newline="") as file:
                writer = csv.writer(file)
                writer.writerows(modified)

for delimeter in csv_delimeters:
    for filename in files:
        if (filename[-4:] == ".csv" and filename[0] == delimeter):
            csvs[csv_delimeters.index(delimeter)] += [filename]

for index in range(3):
    matrices[index] = list(csv.reader(open(csvs[index][0], "rt",encoding="utf8"), delimiter=","))
    for csv_file in csvs[index][1:]:
        matrices[index] += list(csv.reader(open(csv_file, "rt",encoding="utf8"), delimiter=","))[2:]

with open("bond_angle"+ config_data['resid'] + ".csv","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(matrices[0])

with open("bond_length" + config_data['resid'] + ".csv","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(matrices[1])

with open("torsion" + config_data['resid'] + ".csv","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(matrices[2])

for category in csvs:
    for filename in category:
        os.remove(filename)

#clean issues in the torsion output:
def valid_resnos(pdbname):
    with open("bond_length" + config_data['resid'] + ".csv","r") as lengthCsv:
        reader = csv.reader(lengthCsv)
        lengthmatrix = [row for row in reader]
    validresnos=[]    
    i=0
    while(lengthmatrix[i][0][:7]!=pdbname[:7]):
        i+=1
    for row in range(i,len(lengthmatrix)):
        chain=lengthmatrix[row][1][0]
        resno=int(lengthmatrix[row][2])
        if (chain, resno) not in validresnos:
            validresnos.append((chain, resno))
        else:
            break
    return validresnos

def closest_integer(item, resnos):
    closest = 10000000
    chain=item[0]
    for k in resnos:
        if (abs(int(item[1]) - closest) > abs(int(item[1]) - k[1])) and chain==k[0]:
            closest = k[1]
    return closest

def lowestRowOccupancyValue(row):
    lowest=1
    for i in range(4,16,2):
        occ=float(row[i])
        if occ<lowest:
            lowest=float(row[i])
    return lowest

def align_matrix(matrix):
    col_widths = [max(len(str(row[i])) for row in matrix) for i in range(len(matrix[0]))]
    aligned_matrix = [
        [str(row[i]).ljust(col_widths[i]) for i in range(len(row))]
        for row in matrix
    ]
    return aligned_matrix

with open("torsion" + config_data['resid'] + ".csv","r") as torsionCsv:
    reader = csv.reader(torsionCsv)
    torsionmatrix = [row for row in reader]
subtractor=0
prevpdb=""
for row in range(3,len(torsionmatrix)):
    resno=(torsionmatrix[row][1][0], int(torsionmatrix[row][2]))
    pdbname=torsionmatrix[row][0]
    validresnos=valid_resnos(pdbname)

    if resno not in validresnos:
        closestResno=closest_integer(resno, validresnos)
        torsionmatrix[row][2]=closestResno


prevresno=0
for index in range(3,len(torsionmatrix)):
    row=torsionmatrix[index]
    try:
        if int(row[2])==int(prevresno):
            torsionmatrix[index-1][2]=str(row[2])+"A"
            torsionmatrix[index][2]=str(row[2])+"B"

            #make sure the occupancy values line up:
            occ1=lowestRowOccupancyValue(torsionmatrix[index-1])
            occ2=lowestRowOccupancyValue(torsionmatrix[index])
            if occ1+occ2!=1:
                occ1=1-occ2
            
            #apply the occupancy to the whole row
            for i in range(4,18,2):
                torsionmatrix[index-1][i]=round(float(occ1),2)
            for i in range(4,18,2):
                torsionmatrix[index][i]=round(float(occ2),2)
    except:
        None

    prevresno=torsionmatrix[index][2]

#write the corrected torsion matrix
torsionmatrix=align_matrix(torsionmatrix)
os.remove("torsion" + config_data['resid'] + ".csv")
with open("torsion" + config_data['resid'] + ".csv","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(torsionmatrix)

os.chdir(cw_dir)
for i in range(cores):
    os.remove("worker_" + str(i + 1) + ".txt")

print("--- %s seconds ---" % (time.time() - start_time))
