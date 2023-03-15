
import pandas as pd
import numpy as np
import os,argparse
import gzip
import random



tags = ["additive","partial","recessive"]

parser = argparse.ArgumentParser(description="A script for computing summary statistics in 50kb windows across given chromosome, between modern human and archaic human.")
parser.add_argument('-t', '--tag', action="store", dest="tag_id",
                        help="which dominance tag, default: 1",
                        default=1, type=int)

args = parser.parse_args()
tag_this = tags[int(args.tag_id)-1]


def read_largedata (data_path):
    data_iterator = pd.read_csv(data_path, chunksize=1000000,sep="\t",compression='gzip')
    chunk_list = []  
    for data_chunk in data_iterator:  
        filtered_chunk = data_chunk
        chunk_list.append(filtered_chunk)    
    filtered_data = pd.concat(chunk_list)
    return filtered_data

def remove_unnecessary (dataframe,downsized2_path): 
    names = list(dataframe.columns.values)
    print(dataframe.shape)
    for name in names:
        if name[0:7] == "Unnamed":
            dataframe = dataframe.drop([name], axis=1)
            print(dataframe.shape)
    all_0s = dataframe.index[dataframe['classifier'] == 0].tolist()
    count_1s = (dataframe['classifier'] == 1).sum()
    all_1s = dataframe.index[dataframe['classifier'] == 1].tolist()    
    keep = random.sample(all_0s,count_1s*2) + all_1s #keep 2 times more of class 0 observations than class 1
    dataframe = dataframe[dataframe.index.isin(keep)]   
    dataframe.to_csv(downsized2_path, compression="gzip")
    print(dataframe.shape)




#########################################################
DIR = "/u/scratch/x/xinjunzh/slim/batch1/"

os.chdir(DIR)


newfile_path = "added_"+tag_this+"_1kseg.txt.gz"


downsized2_path = "downsized2-1_"+newfile_path #create 2:1 ratio downsampled simulation pile

dataframe = read_largedata (newfile_path)
print("done reading \n")

remove_unnecessary (dataframe,downsized2_path)
print("done processing \n")






