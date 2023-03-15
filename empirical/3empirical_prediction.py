
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

def read_largedata (data_path):
    data_iterator = pd.read_csv(data_path, chunksize=1000000,sep=",",compression='gzip')
    chunk_list = []  
# Each chunk is in dataframe format
    for data_chunk in data_iterator:  
        #filtered_chunk = chunk_preprocessing(data_chunk)
        filtered_chunk = data_chunk
        chunk_list.append(filtered_chunk)    
    filtered_data = pd.concat(chunk_list)
    return filtered_data

def read_largecsv (data_path):
    data_iterator = pd.read_csv(data_path, chunksize=1000000,sep=",",compression='gzip')
    chunk_list = []  
# Each chunk is in dataframe format
    for data_chunk in data_iterator:  
        #filtered_chunk = chunk_preprocessing(data_chunk)
        filtered_chunk = data_chunk
        chunk_list.append(filtered_chunk)    
    filtered_data = pd.concat(chunk_list)
    return filtered_data



#########################################################

DIR = "/u/project/kirk-bigdata/xinjunzh/stats_overlap_19pops_ready/"

os.chdir(DIR)
data_list = os.listdir()
data_list = [f for f in data_list if "_added.csv" in f]
data_list_nea = [f for f in data_list if "nea" in f]
data_list_den = [f for f in data_list if "den" in f]

#data_list = data_list_nea
data_list = data_list_nea+data_list_den
feature_list = ["set4-allbutQ.txt"]

DIR_model = "/u/scratch/x/xinjunzh/slim/pretrained_model/"


for feature in feature_list:
	feature_file = "feature/"+feature
	model_name = DIR_model+"MaLAdapt_pretrained_model.sav"
	ML = pickle.load(open(model_name, "rb"))
	savename = "trained-model"

	for datapath in data_list:
		dataframe = pd.read_csv(datapath)    
		names = list(dataframe.columns.values)
		for name in names:
			if name[0:7] == "Unnamed":
				dataframe = dataframe.drop([name], axis=1)
		names = list(dataframe.columns.values)[5:] 

		with open(feature_file,"rt") as file:
			these = file.readline().split()
	
		exclude = [i for i in names if i not in these]    #remove unused features
		dataframe = dataframe.drop(exclude, axis=1)        
		dataframe = dataframe.replace(float("inf"), 0) 
		dataframe = dataframe.fillna(0)
		dataframe = dataframe.replace(str("Nan"), 0)
		array = dataframe.values
		test_X = array[:,:]
		scaler = StandardScaler()
		scaler.fit(test_X[:,5:dataframe.shape[1]])
		test_X[:,5:dataframe.shape[1]] = scaler.transform(test_X[:,5:dataframe.shape[1]])
	
		pred = ML.predict(test_X[:,5:dataframe.shape[1]].astype(float))
		y_score = ML.predict_proba(test_X[:,5:dataframe.shape[1]].astype(float))    
		actual_Y = pred
		z = np.zeros((test_X.shape[0],1))
		all_data = np.append(test_X,z,1)
		all_data[:,all_data.shape[1]-1] = actual_Y
		all_data = np.concatenate((all_data,y_score),axis=1)
	
		header = list(dataframe.columns.values)+["classifier","prob0","prob1"]
		pd.DataFrame(all_data).to_csv("/u/scratch/x/xinjunzh/predictions/predicted_AI+sweep_"+savename+"_"+datapath,header=header)
		print(datapath)
	
