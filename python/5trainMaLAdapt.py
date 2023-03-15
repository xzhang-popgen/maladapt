
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

parser = argparse.ArgumentParser(description="A script for computing summary statistics in 50kb windows across given chromosome, between modern human and archaic human.")
parser.add_argument('-m', '--model', action="store", dest="model_id",
                        help="which simulation batch, default: 1",
                        default=1, type=int)

args = parser.parse_args()

whichmodel = int(args.model_id) #now sim id is going to refer to a 


def read_largedata (data_path):
    data_iterator = pd.read_csv(data_path, chunksize=1000000,sep="\t",compression='gzip')
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



###############################################################

def prec_recall_mod (cm, num_classes): #this is to get precision recall for multiple flanking classes
    prec_list = []
    reca_list = []    
    for i in range(0,num_classes):
        prec_i = round(cm[i,i]/sum(cm[:,i]),4)
        reca_i = round(cm[i,i]/sum(cm[i,:]),4)
        prec_list.append(prec_i)
        reca_list.append(reca_i)
    prec_list = [0 if np.isnan(x) else x for x in prec_list]
    reca_list = [0 if np.isnan(x) else x for x in reca_list]    
    print(prec_list)
    print(reca_list)    
    return prec_list,reca_list

def plot_prcurve_average_mod(testY,yscore,savepath,num_classes):
    y_score = yscore
    ydata = testY #binarized
    precision = dict()
    recall = dict()
    average_precision = dict()
    thresh = dict()
    label_list = [1,0]
    label_name = ["non-AI","AI"]
    for i in range(num_classes):
        precision[i], recall[i], thresh[i] = precision_recall_curve(ydata[:, i],y_score[:, i],pos_label=1)
        average_precision[i] = average_precision_score(ydata[:, i], y_score[:, i])
    precision["micro"], recall["micro"], thresh["micro"] = precision_recall_curve(ydata.ravel(),y_score.ravel())
    average_precision["micro"] = average_precision_score(ydata, y_score,average="micro")
    # setup plot details
    #colors = cycle(['blue', 'red', 'yellow'])
    colors = ["b","r"] + sns.color_palette("Greens_r",num_classes)[2:num_classes]
    #sns.palplot(colors)
    #colors = sns.palplot('blue', 'red',sns.color_palette("Oranges_r",7)[2:7])
    #choose several shades of yellows
    plt.figure(figsize=(7, 8))
    f_scores = np.linspace(0.2, 0.8, num=4)
    lines = []
    labels = []
    for f_score in f_scores:
        x = np.linspace(0.01, 1)
        y = f_score * x / (2 * x - f_score)
        l, = plt.plot(x[y >= 0], y[y >= 0], color='gray', alpha=0.2)
    for i, color in zip(range(num_classes), colors):
        l, = plt.plot(recall[i], precision[i], color=color, lw=2)
        lines.append(l)
        labels.append('Precision-Recall for class {0} (area = {1:0.2f})'
                  ''.format(label_name[i], average_precision[i]))
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall curve of all classes')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=10))
    fig.savefig(savepath+'.png')
    plt.show()



def plot_roc_mod(testY,yscore,savepath,num_classes):
    y_score = yscore
    ydata = testY #binarized
    fpr = dict()
    tpr = dict()
    thresh = dict()
    label_list = [1,0]
    label_name = ["non-AI","AI"]
    for i in range(num_classes):
        fpr[i], tpr[i], thresh[i] = roc_curve(ydata[:, i],y_score[:, i],pos_label=1)
    colors = ["b","r"] + sns.color_palette("Greens_r",num_classes)[2:num_classes]
    plt.figure(figsize=(7, 8))
    lines = []
    labels = []
    for i, color in zip(range(num_classes), colors):
        plt.plot([0, 1], [0, 1], 'k--')
        l, = plt.plot(fpr[i], tpr[i], color=color, lw=2)
        lines.append(l)
        labels.append('ROC for class {0}'
                  ''.format(label_name[i]))
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25)
    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('FPR')
    plt.ylabel('TPR')
    plt.title('ROC curve of all classes')
    plt.legend(lines, labels, loc=(0, -.38), prop=dict(size=10))    
    fig.savefig(savepath+'.png')    
    plt.show()


def write_importance (outfile,feature_list,score_list):
    with open(outfile, "w") as feature_file:
        feature_name = feature_list
        score = score_list
        feature_file.writelines(i+"\t" for i in feature_name)
        feature_file.write("\n")
        feature_file.writelines(str(i)+"\t" for i in score)
        feature_file.write("\n")

def plot_featurescore (scores_file,savefile):
    with open(scores_file,"rt") as file:
        names = file.readline()
        scores = file.readline()        
    names = names.split("\t")
    names = names[0:len(names)-1]
    scores = scores.split("\t")
    scores = scores[0:len(scores)-1]  
    scores   = [float(i) for i in scores]    
    scores, names = zip(*sorted(zip(scores, names)))    
    fig = plt.gcf()
    fig.subplots_adjust(bottom=0.25,left=0.25)
    plt.barh(range(len(names)),scores, align="center")
    plt.yticks(range(len(names)),names,fontsize=8)
    plt.xlabel("Feature Importance",fontsize=16)
    plt.ylabel("Features",fontsize=16)
    plt.title("Feature Importance Ranking")
    plt.show()
    fig.savefig(savefile,orientation='landscape') 


def ML_cm(dataframe,whichML,savename):
    ETC = ExtraTreesClassifier(n_estimators=100, random_state=0,max_features="sqrt",min_samples_split=10,min_samples_leaf=10) #saga, sag, lbfgs,liblinear
    RF = RandomForestClassifier(n_estimators=100, random_state=0,max_depth=2) #saga, sag, lbfgs,liblinear
    L0LR = LogisticRegression(solver = 'saga', multi_class='multinomial',penalty='l1',tol=0.01,C=1e10)
    L1LR = LogisticRegression(solver = 'saga', multi_class='multinomial',penalty='l1')
    L2LR = LogisticRegression(solver = 'lbfgs', multi_class='multinomial',penalty='l2')
    MLP = MLPClassifier(solver='lbfgs', alpha=1e-5, hidden_layer_sizes=(512, 256, 128), random_state=1)    
    RBF = SVC(kernel='rbf',class_weight='balanced',probability=True) #kernal, penalized    
    if (whichML=="ETC"):
        ML = ETC
    elif (whichML=="RF"):
        ML = RF
    elif (whichML=="L0LR"):
        ML = L0LR
    elif (whichML=="L1LR"):
        ML = L1LR       
    elif (whichML=="L2LR"):
        ML = L2LR        
    elif (whichML=='MLP'):
        ML = MLP   
    elif (whichML=='RBF'):
        ML = RBF
    names = list(dataframe.columns.values)
    for name in names:
        if name[0:7] == "Unnamed":
            dataframe = dataframe.drop([name], axis=1)
    array = dataframe.values
    X = array[:,:]
    scaler = StandardScaler()
    scaler.fit(X[:,13:dataframe.shape[1]])
    X[:,13:dataframe.shape[1]] = scaler.transform(X[:,13:dataframe.shape[1]])
    Y = array[:,0]
    Y=Y.astype(int)
    train_X, test_X, train_Y, test_Y = train_test_split( X, Y, test_size=1/8.0, random_state=0)
    ML.fit(train_X[:,13:X.shape[1]].astype(float), train_Y)
    filename = DIR_save+whichML+'feature-'+savename+'_finalized_model.sav'
    pickle.dump(ML, open(filename, 'wb'))    
    pred = ML.predict(test_X[:,13:X.shape[1]].astype(float))
    ML.score(test_X[:,13:X.shape[1]].astype(float), test_Y) 
    cm = metrics.confusion_matrix(test_Y, pred)
    print(cm)
    misclassified = pred!=test_Y
    misclassified_X = test_X[misclassified,:]
    misclassified_Y = test_Y[misclassified]
    z = np.zeros((misclassified_X.shape[0],1))
    mis_all = np.append(misclassified_X,z,1)
    mis_all[:,mis_all.shape[1]-1] = misclassified_Y        
    actual_Y = pred[misclassified]
    z = np.zeros((mis_all.shape[0],1))
    mis_all = np.append(mis_all,z,1)
    mis_all[:,mis_all.shape[1]-1] = actual_Y
    np.savetxt(DIR_save+whichML+"/feature-misclassified_all-all_"+savename+".csv", mis_all, delimiter=",", fmt='%s')
    okclassified = pred==test_Y
    okclassified_X = test_X[okclassified,:]
    okclassified_Y = test_Y[okclassified]        
    z = np.zeros((okclassified_X.shape[0],1))
    ok_all = np.append(okclassified_X,z,1)
    ok_all[:,ok_all.shape[1]-1] = okclassified_Y        
    actual_Y = pred[okclassified]
    z = np.zeros((ok_all.shape[0],1))
    ok_all = np.append(ok_all,z,1)
    ok_all[:,ok_all.shape[1]-1] = actual_Y        
    np.savetxt(DIR_save+whichML+"/feature-okclassified_all-all_"+savename+".csv", ok_all, delimiter=",", fmt='%s')     
    np.savetxt(DIR_save+whichML+"/feature-confusionmatrix_"+savename+".csv", cm, delimiter=",") #save confusion matrix   
    prec,reca = prec_recall_mod(cm,2)
    prec_reca_mat = np.row_stack((prec,reca))
    np.savetxt(DIR_save+whichML+"/feature-precision-recall_"+savename+".csv", prec_reca_mat, delimiter=",")    
    y_score = ML.predict_proba(test_X[:,13:X.shape[1]].astype(float))
    ydata = label_binarize(test_Y, classes=[*range(3)])[:,0:2]    
    plot_prcurve_average_mod(ydata,y_score,DIR_save+whichML+"/feature-precision-recall_"+savename,2)      
    plot_roc_mod(ydata,y_score,DIR_save+whichML+"/feature-roc_"+savename+"_01",2)
    result = np.concatenate((ydata,y_score),axis=1)[:,1:]
    header = ["classifier","prob0","prob1"]
    pd.DataFrame(result).to_csv(DIR_save+whichML+"/prediction-prob_"+whichML+savename+".txt.gz",header=header,sep="\t",index=False,compression=("gzip"))
    if whichML == "ETC":
        print(ML.feature_importances_)
        score = ML.feature_importances_.tolist()
        feature_names = list(dataframe.columns.values)[13:X.shape[1]]
        features = feature_names#[1:]
        write_importance (DIR_save+whichML+"/feature-feature_importance_"+savename+".txt",features,score)        
        plot_featurescore(DIR_save+whichML+"/feature-feature_importance_"+savename+".txt",DIR_save+whichML+"/feature_importance_"+savename+".png")

    y_mod_0 = np.zeros((ydata.shape[0]))#,1))
    y_mod_0[y_score[:, 1]>=0.9] = 1
    real_Y = test_Y    
    TP = sum((y_mod_0==1) & (real_Y==1))
    TN = sum((y_mod_0==0) & (real_Y==0))
    FP = sum((y_mod_0==1) & (real_Y==0))
    FN = sum((y_mod_0==0) & (real_Y==1))    
    rec = TPR =TP/(TP+FN)
    FPR = FP/(FP+TN)
    prec = TP/(TP+FP)
    FDR=1-prec
    fscore = (2 * prec * rec) / (prec + rec)
    with open(DIR_save+whichML+"/metrics_summary_"+whichML+savename+".txt","w") as f:
        f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(FPR,TPR,prec,rec,FDR,fscore))


#########################################################
DIR = "/u/scratch/x/xinjunzh/slim_nonAIsweep/data_AI+sweep/"

os.chdir(DIR)


downsized2_path = "downsized2-1_AI-sweep.txt.gz" #load the downsampled simulation data; class ratio 2:1

#whichML_list = ["ETC","RF","L0LR","L1LR","L2LR"]
whichML_list = ["ETC"]


#feature_list = ["set4-allbutQ.txt","set5-allbutQmax.txt","set6-all.txt","set1-highrank.txt","set2-highrank_noQmax.txt","set3-midrank.txt"]
feature_list = ["feature_list/set4-allbutQ.txt"] #use the selected features

for feature in feature_list:
    datapath = downsized2_path
    dataframe = read_largecsv (datapath)
        
    names = list(dataframe.columns.values)[14:]
    
    feature_file = "../feature/"+feature
    
    with open(feature_file,"rt") as file:
        these = file.readline().split()
    
    exclude = [i for i in names if i not in these]
    
    dataframe = dataframe.drop(exclude, axis=1)
    
    print(list(dataframe.columns.values))
        
    dataframe = dataframe.replace(float("inf"), 0) #there shouldn't be Inf values, but if there are, double check if filling in 0s is ok
    dataframe = dataframe.fillna(0) #fill the nan with 0; ok because this applies to only the Q/U stats
    dataframe = dataframe.replace(str("Nan"), 0) #fill the "Nan"s with 0; ok because this would apply only to ZnS
    dataframe = dataframe.sample(frac=1) #shuffle the table
    dataframe = dataframe[dataframe.classifier != "classifier"]
    DIR_save = "../model_AI+sweep/"
    whichML=whichML_list[whichmodel-1]
    savename = "AI+sweep_"+feature.split(".")[0]
    ML_cm(dataframe,whichML,savename)


