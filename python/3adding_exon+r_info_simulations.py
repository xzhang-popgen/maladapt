
import pandas as pd
import numpy as np
from sklearn import preprocessing
import matplotlib.pyplot as plt
from inspect import signature
import seaborn as sns
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2
from sklearn.feature_selection import RFE
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso
import os,argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import label_binarize
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from itertools import cycle
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from itertools import cycle
import pickle

from sklearn.svm import SVC
from sklearn.decomposition import PCA

import gzip
import random

whichbatch = "batch1"

DIR = "/u/scratch/x/xinjunzh/slim/"+whichbatch+"/"
os.chdir(DIR)

tags = ["additive","partial","recessive"]

parser = argparse.ArgumentParser(description="A script for computing summary statistics in 50kb windows across given chromosome, between modern human and archaic human.")
parser.add_argument('-t', '--tag', action="store", dest="tag_id",
                        help="which dominance tag, default: 1",
                        default=1, type=int)

args = parser.parse_args()
tag_this = [tags[int(args.tag_id)-1]]


tag_this=1

#for tag in tag_this:
for tag in tags:
	data_path0 = tag+"_1kseg.txt.gz" 
	DIR_info =  "/u/scratch/x/xinjunzh/slim/1ksegments/info_exon+r/"
	newfile_path = "added_"+data_path0 #+ ".txt"
	with gzip.open(data_path0, 'rt') as f:
		count = 0
		with gzip.open(newfile_path, "wt") as new_file:
			for line in f:
				count +=1
				field = line.split()
				#field=[x.decode('utf-8') for x in field]
				if count ==1:
					newfield = field+["exon","r"]
					new_file.writelines(i+"\t" for i in newfield)
					new_file.write("\n")				
				else:
					if field[0] != "rep":
						seg = field[2]
						start = field[10]
						end = field[11]				
						info_file = DIR_info + "exon+r_sim_seq_info_" + seg + ".txt"
						with open(info_file,"rt") as info:
							for l in info:
								these = l.split()
								if (str(these[0]) == start) & (str(these[1]) ==end):
									newfield = field+ [str(these[2]),str(these[3])]
									break
						new_file.writelines(i+"\t" for i in newfield)
						new_file.write("\n")
	
print("DONE ADDING")
