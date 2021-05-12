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

for delimeter in csv_delimeters:
    for filename in files:
        if (filename[-4:] == ".csv" and filename[0] == delimeter):
            csvs[csv_delimeters.index(delimeter)] += [filename]

for index in range(3):
    matrices[index] = list(csv.reader(open(csvs[index][0], "rt",encoding="utf8"), delimiter=","))
    for csv_file in csvs[index][1:]:
        matrices[index] += list(csv.reader(open(csv_file, "rt",encoding="utf8"), delimiter=","))[2:]

with open("bond_angle","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(matrices[0])

with open("bond_length","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(matrices[1])

with open("bond_length","w",newline='') as my_csv:
    newarray = csv.writer(my_csv,delimiter=',')
    newarray.writerows(matrices[2])

for category in csvs:
    for filename in category:
        os.remove(filename)

os.chdir(cw_dir)
for i in range(cores):
    os.remove("worker_" + str(i + 1) + ".txt")

print("--- %s seconds ---" % (time.time() - start_time))
