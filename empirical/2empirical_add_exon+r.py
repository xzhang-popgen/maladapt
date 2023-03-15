
import itertools
import os
import numpy as np
import pandas as pd
import gzip

#awk ' FNR==1 && NR!=1 { while (/^chr/) getline; } 1 {print}' *.txt >autosomes_windows_overlap.txt.added.txt

DIR_stat = "/u/project/kirk-bigdata/xinjunzh/stats_overlap_19pops/"

os.chdir(DIR_stat)

allfiles = os.listdir()
allfiles = [f for f in allfiles if ".csv" in f]
annotation = "../annotations/autosomes_windows_overlap.txt.added.txt"

info = pd.read_csv(annotation,sep="\t")

for raw_path in allfiles:
    dataframe = pd.read_csv (raw_path)
    names = list(dataframe.columns.values)
    for name in names:
        if name[0:7] == "Unnamed":
            dataframe = dataframe.drop([name], axis=1)   
    dataframe=dataframe.sort_values(["chr",'start'],ascending=[True,True])    
    new = dataframe.merge(info,how='left',on=["chr","start","end","length","n_snps"])    
    new.to_csv("stats_overlapp-19pops-added/"+raw_path.split(".")[0]+"_added.csv",index=False)
    print(raw_path)



######concat chromosomes
allpop = ['CHB','KHV','GIH','PJL','CHS','CEU','GBR','IBS','PEL','TSI','JPT','CDX','FIN','MXL','PUR','CLM','BEB','STU','ITU']
os.chdir("/u/project/kirk-bigdata/xinjunzh/stats_overlap_9pops/stats_overlapp-19pops-added/")

all = os.listdir()
all = [f for f in all if ".csv" in f]

for p in allpop:
	popnea = [f for f in all if "nea"+p in f]
	df1 = pd.read_csv(popnea[0], header = 0,sep=",")
	names = list(df1.columns.values)
	for name in names:
		if name[0:7] == "Unnamed":
			df1 = df1.drop([name], axis=1)
	df=df1
	for file in popnea[1:]:
		df2 = pd.read_csv(file, header = 0,sep=",")
		names = list(df2.columns.values)
		for name in names:
			if name[0:7] == "Unnamed":
				df2 = df2.drop([name], axis=1)
			df = pd.concat([df,df2])
	df.to_csv("stats_overlapp-19pops-added_allchrom/"+"nea"+p+"_ALL_added.csv",index=False)
	print(p)

for p in allpop:
	popden = [f for f in all if "den"+p in f]
	df1 = pd.read_csv(popden[0], header = 0,sep=",")
	names = list(df1.columns.values)
	for name in names:
		if name[0:7] == "Unnamed":
			df1 = df1.drop([name], axis=1)
	df=df1
	for file in popden[1:]:
		df2 = pd.read_csv(file, header = 0,sep=",")
		names = list(df2.columns.values)
		for name in names:
			if name[0:7] == "Unnamed":
				df2 = df2.drop([name], axis=1)
			df = pd.concat([df,df2])
	df.to_csv("stats_overlapp-19pops-added_allchrom/"+"nea"+p+"_ALL_added.csv",index=False)
	print(p)
