import pandas as pd
import numpy as np
import os

DIR_add = "/u/scratch/x/xinjunzh/slim/additive/stat/"
DIR_part = "/u/scratch/x/xinjunzh/slim/partial/stat/"
DIR_rec = "/u/scratch/x/xinjunzh/slim/recessive/stat/"



DIRs = [DIR_add,DIR_part,DIR_rec]
tags = ["additive","partial","recessive"]

whichbatch = "batch1"

for t in range(0,3):#range(0,3)
	os.chdir(DIRs[t])
	genes = [name for name in os.listdir(".") if os.path.isdir(name)]
	n = 0
	print(tags[t])
	for g in range(0,len(genes)):
		gene = genes[g]
		DIR = DIRs[t]
		tag = tags[t]    
		dir_folder = DIR+gene + "/"
		os.chdir(dir_folder)
		command = "awk ' FNR==1 && NR!=1 { while (/^rep/) getline; } 1 {print}' *_"+whichbatch[5]+".txt > ../"+tag+"-"+gene+"_1ksegs.txt"
		os.system(command)
		
		n+=1
		print(n)
        #os.system('rm *.txt')


for t in range(0,3):
	DIR = DIRs[t]
	tag = tags[t]
	dir_folder = DIR
	os.chdir(dir_folder)
	command = "awk ' FNR==1 && NR!=1 { while (/^rep/) getline; } 1 {print}' *.txt > /u/scratch/x/xinjunzh/slim/"+tag+"_1kseg.txt"
	os.system(command)

print("DONE COMPILE")

for t in range(0,3):#range(0,3)
	tag = tags[t]
	name = "/u/scratch/x/xinjunzh/slim/"+tag+"_1kseg.txt"
	command = "gzip -c "+name +">"+ name+".gz"
	os.system(command)

print("DONE COMPRESS")


#os.chdir("/u/scratch/x/xinjunzh/slim/"+whichbatch)
#command = "awk ' FNR==1 && NR!=1 { while (/^rep/) getline; } 1 {print}' *1kseg.txt > all_sweep_"+whichbatch+".txt"
#os.system(command)



