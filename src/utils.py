# coding:utf-8

import numpy as np
import pandas as pd
import os
import subprocess
import scipy.io
import pickle
from time import sleep

import requests
from bs4 import BeautifulSoup

import pubchempy as pcp

###########################
## VERBOSE               ##
###########################

def print_dataset(ratings, user_col, item_col, rating_col):
    '''
        @param\tratings\tPandas Dataframe: 3 columns (user_col x item_col x rating_col) that lists all non zero values {-1,1}.
        @param\tuser_col\tPython character string: how to name the column associated with users.
        @param\titem_col\tPython character string: how to name the column associated with items.
        @param\trating_col\tPython character string: how to name the column associated with ratings.
    '''
    ratings2 = ratings.copy()
    ratings2[item_col] = ratings2[item_col].astype(str)
    args = [len(np.unique(ratings2[item_col])), len(np.unique(ratings2[user_col]))]
    args += [ratings2.loc[ratings2[rating_col]==v].shape[0] for v in [1,-1,0]]
    assert args[2]+args[3]+args[4]==ratings2.shape[0]
    args[-1] = args[-1] if (args[-1]>0) else args[0]*args[1]-args[2]-args[3]
    dataset_str = "Ndrugs=%d\tNdiseases=%d\n%d positive\t%d negative\t%d unknown matchings"
    print(dataset_str % tuple(args))
    
def compute_sparsity(df):
    '''
        @param\tdf\tPandas Dataframe: size #items x #users containing values {-1,0,1}.
        @return\tsparsity\tPython int: 100 times the number of non zero values divided by the total number of values in @df.
    '''
    return 100*(np.sum(df.values!=0)/(df.shape[0]*df.shape[1]))

#############################
## HANDLING RATINGS MATRIX ##
#############################

def matrix2ratings(df, user_col="user", item_col="item", rating_col="rating"):
    '''
        @param\tdf\tPandas Dataframe: size #items x #users containing values {-1,0,1}.
        @param\tuser_col\tPython character string: how to name the column associated with users.
        @param\titem_col\tPython character string: how to name the column associated with items.
        @param\trating_col\tPython character string: how to name the column associated with ratings.
        @return\tres_df\tPandas Dataframe: 3 columns (user_col x item_col x rating_col) that lists all non zero values.
    '''
    assert all([a in [-1,0,1] for a in np.unique(df.values.flatten())])
    non_missing = np.argwhere(df.values!=0)
    res_df = pd.DataFrame([], index=range(non_missing.shape[0]))
    res_df[user_col] = [df.columns[x] for x in list(non_missing[:, 1].flatten())]
    res_df[item_col] = [df.index[x] for x in list(non_missing[:,0].flatten())]
    res_df[rating_col] = [v for v in df.values[non_missing[:,0], non_missing[:,1]].flatten()]
    return res_df

def load_dataset(model_name, save_folder="./"):
    '''
        @param\tmodel_name\tPython character string: Name of the dataset to load.
        @param\tsave_folder\tPython character string[default="./"]: (Relative) path to which the dataset should be stored.
        @return\tdi\tPython dictionary: {"A": drug-disease matrix Nd x Np, "P": disease feature matrix Fp x Np, "S": drug feature matrix Fd x Nd}
        P and S are real-valued matrices
        A has values in {-1,0,1}. -1 means negative matching, 1: positive matching, 0: unknown matching.
    '''
    assert model_name in ["Gottlieb", "Cdataset_Aonly", "indep", "Fdataset", "DNdataset", "Cdataset", "TRANSCRIPT", "PREDICT", "LRSSL", "PREDICT_Gottlieb"]
    if (model_name in ["LRSSL", "PREDICT_Gottlieb", "OrphanDrug"]):
        print("Warning: this dataset has no drug/disease names!")
        url_dda_skf = "https://github.com/GCQ2119216031/DDA-SKF/raw/master/data/"
        dda_skf_dataset_path = save_folder+"DDA_SKF/data/"
        mmodel_name = model_name if (model_name!="PREDICT_Gottlieb") else "PREDICT"
        if (not os.path.exists(dda_skf_dataset_path+mmodel_name+".mat")):
            subprocess.call(" ".join(['mkdir', '-p', dda_skf_dataset_path]), shell=True)
            subprocess.call(" ".join(["wget", "-qO", dda_skf_dataset_path, url_dda_skf+mmodel_name+".mat"]), shell=True)
        ddt = scipy.io.loadmat(dda_skf_dataset_path+mmodel_name+".mat")
        if (mmodel_name=="LRSSL"):
            A = ddt["lrssladmatdgc"]
            S_chemical = ddt["lrsslsimmatdcchemical"]
            S_sideeffects = ddt["lrsslsimmatdcgo"]
            S = np.concatenate((S_chemical,S_sideeffects), axis=0)
            S.index = [("chemical--" if (iss < S_chemical.shape[0]) else "se--")+str(s) for iss, s in enumerate(list(S.index))]
            P = ddt["lrsslsimmatdg"]
        elif (mmodel_name=="PREDICT"):
            A = ddt['predictAdMatdgc'].T
            S_chemical = ddt['predictSimMatdcChemical']
            S_domain = ddt['predictSimMatdcDomain']
            S_GO = ddt['predictSimMatdcGo']
            S = np.concatenate((S_chemical,S_domain,S_GO), axis=0)
            S.index = [("chemical--" if (iss < S_chemical.shape[0]) else ("domain--" if (iss < S_chemical.shape[0]+S_domain.shape[0]) else "go--"))+str(s) for iss, s in enumerate(list(S.index))]
            P = ddt['predictSimMatdg']
        else:
            raise ValueError("Undefined dataset '%s'" % mmodel_name)
        P = pd.DataFrame(P, index=range(P.shape[0]), columns=["disease%d" % (i+1) for i in range(P.shape[1])])
        S = pd.DataFrame(S, index=range(S.shape[0]), columns=["drug%d" % (i+1) for i in range(S.shape[1])])
        A = pd.DataFrame(A, index=S.columns, columns=P.columns)
    if (model_name in ["TRANSCRIPT", "PREDICT"]):
        path=save_folder+model_name+"/" 
        fnames = {"A": "ratings_mat.csv", "P": "users.csv", "S": "items.csv"} 
        if (not os.path.exists(path+fnames["A"])):
            if (model_name == "TRANSCRIPT"):
                url_dataset = "https://zenodo.org/record/7982976/files/TRANSCRIPT_dataset_v2.0.0.zip"
            else:
                print("Run the notebook or use the publicly available data (by default, download the latest public dataset)")
                url_dataset = "https://zenodo.org/record/7983090/files/PREDICT_dataset_v2.0.0.zip"
            subprocess.call(" ".join(["wget", "-qO", save_folder+model_name+".zip", "\""+url_dataset+"\""]), shell=True)  
            subprocess.call(" ".join(["unzip", save_folder+model_name+".zip", "&&", "mv", model_name+"_dataset_v2.0.0", save_folder+model_name, "&&", "rm", "-f", save_folder+model_name+".zip"]), shell=True)   
        if (model_name == "TRANSCRIPT"):
            assert os.path.exists(path+fnames["A"])
            assert os.path.exists(path+fnames["P"])
            assert os.path.exists(path+fnames["S"])
            A, P, S = [pd.read_csv(path+fnames[k], engine="python", index_col=0) for k in ["A", "P", "S"]]
        else:
            assert os.path.exists(path+fnames["A"])
            if (os.path.exists(path+fnames["P"]) and os.path.exists(path+fnames["S"])):
                A, P, S = [pd.read_csv(path+fnames[k], engine="python", index_col=0) for k in ["A", "P", "S"]]
            else: ## use publicly available data 
                A = pd.read_csv(path+fnames["A"], engine="python", index_col=0)
                P_phenotype = pd.read_csv(path+"disease_phenotype_PREDICT_matrix.csv", engine="python", index_col=0)
                P_semantic = pd.read_csv(path+"disease_semantic_PREDICT_matrix.csv", engine="python", index_col=0)
                S_se = pd.read_csv(path+"se_PREDICT_matrix.csv", engine="python", index_col=0)
                S_signature = pd.read_csv(path+"signature_PREDICT_matrix.csv", engine="python", index_col=0)
                P = P_phenotype.T.join(P_semantic.T, how="outer").T
                P.index = [("phenotype--" if (iss < P_phenotype.shape[0]) else "semantic--")+str(s) for iss, s in enumerate(list(P.index))]
                S = S_se.T.join(S_signature.T, how="outer").T
                S.index = [("se--" if (iss < S_se.shape[0]) else "signature--")+str(s) for iss, s in enumerate(list(S.index))]
                A = A.loc[S.columns][P.columns]
    if (model_name == "Gottlieb"):
        url_mbirw = "https://raw.githubusercontent.com/bioinfomaticsCSU/MBiRW/master/Datasets/"
        gottlieb_dataset_path = save_folder+"Gottlieb_dataset/MBiRW_files/"
        if (not os.path.exists(gottlieb_dataset_path+"DiDrAMat")):
            subprocess.call(" ".join(['mkdir', '-p', gottlieb_dataset_path]), shell=True)
            for fname in ["DiDrAMat","DiseaseSimMat","DiseasesName","DrugSimMat","DrugsName"]:
                subprocess.call(" ".join(["wget", "-qO", gottlieb_dataset_path+fname, url_mbirw+fname]), shell=True)
        fnames = {"A": "DiDrAMat", "P": "DiseaseSimMat", "S": "DrugSimMat"}
        names_fnames = {"drugs": "DrugsName", "diseases": "DiseasesName"}
        P, S = [pd.read_csv(gottlieb_dataset_path+fnames[k], sep=" ", header=None) for k in ["P", "S"]]
        A = pd.read_csv(gottlieb_dataset_path+fnames["A"], sep="\t", header=None).iloc[:,:-1].T
        drug_names = pd.read_csv(gottlieb_dataset_path+names_fnames["drugs"], header=None).T.values.tolist()
        drug_names = [x for x in drug_names]
        disease_names = pd.read_csv(gottlieb_dataset_path+names_fnames["diseases"], header=None).T.values.tolist()
        disease_names = [x for x in disease_names]
        A.index = drug_names
        S.index = drug_names
        S.columns = drug_names
        S.index = [x[0] for x in S.index]
        S.columns = [x[0] for x in S.columns]
        A.columns = disease_names
        A.index = [x[0] for x in A.index]
        A.columns = [x[0] for x in A.columns]
        P.index = disease_names
        P.columns = disease_names
        P.index = [x[0] for x in P.index]
        P.columns = [x[0] for x in P.columns]
    if (model_name == "indep"):
        url_mbirw = "https://raw.githubusercontent.com/bioinfomaticsCSU/MBiRW/master/Datasets/Datasets_indep/"
        indep_dataset_path = save_folder+"Dataset_indep/MBiRW_files/"
        if (not os.path.exists(indep_dataset_path+"DiDrMat.mat")):
            subprocess.call(" ".join(['mkdir', '-p', indep_dataset_path]), shell=True)
            for fname in ["DiDrMat.mat","R_Wdname","R_Wrname"]:
                subprocess.call(" ".join(["wget", "-qO", indep_dataset_path+fname, url_mbirw+fname]), shell=True)
        A = scipy.io.loadmat(indep_dataset_path+"DiDrMat.mat")['R_Wdr']
        A = pd.DataFrame(A, index=range(A.shape[0]), columns=range(A.shape[1])).T
        drug_names = pd.read_csv(indep_dataset_path+"R_Wrname", header=None).T.values.tolist()[0]
        disease_names = pd.read_csv(indep_dataset_path+"R_Wdname", header=None).T.values.tolist()[0]
        A.index = drug_names
        A.columns = disease_names
        P = pd.DataFrame([], index=A.columns).T
        S = pd.DataFrame([], index=A.index).T
    if (model_name == "Cdataset_Aonly"):
        url_mbirw = "https://raw.githubusercontent.com/bioinfomaticsCSU/MBiRW/master/Datasets/CDatasets/"
        cdataset_dataset_path = save_folder+"Cdatasets/MBiRW_files/"
        if (not os.path.exists(cdataset_dataset_path+"DiDrMat")):
            subprocess.call(" ".join(['mkdir', '-p', cdataset_dataset_path]), shell=True)
            for fname in ["DrDiMat","DiseasesName","DrugsName"]:
                subprocess.call(" ".join(["wget", "-qO", cdataset_dataset_path+fname, url_mbirw+fname]), shell=True)
        A = pd.read_csv(cdataset_dataset_path+"DrDiMat", sep="\t", header=None)
        A = A.T.dropna(how="any").T
        drug_names = pd.read_csv(cdataset_dataset_path+"DrugsName", header=None).T.values.tolist()[0]
        disease_names = pd.read_csv(cdataset_dataset_path+"DiseasesName", header=None).T.values.tolist()[0]
        A.index = drug_names
        A.columns = disease_names
        P = pd.DataFrame([], index=A.columns).T
        S = pd.DataFrame([], index=A.index).T
    elif (model_name in ["Cdataset", "Fdataset", "DNdataset"]):
        drrs_dataset_path = save_folder+"Cdatasets/DRRS_files/"
        url_drrs = "http://bioinformatics.csu.edu.cn/resources/softs/DrugRepositioning/DRRS/soft/"
        subprocess.call(" ".join(["mkdir", "-p", drrs_dataset_path]), shell=True)
        if (not os.path.exists(drrs_dataset_path+model_name+"s/DiDrA.txt")):
            if (not os.path.exists(drrs_dataset_path+model_name+"s.zip")):
                subprocess.call(" ".join(["wget", "-qO", drrs_dataset_path+model_name+"s.zip", url_drrs+model_name+"s.zip"]), shell=True)
            subprocess.call(" ".join(["unzip", "-d", drrs_dataset_path, drrs_dataset_path+model_name+"s.zip"]), shell=True)
        if (model_name not in ["Cdataset", "Fdataset"]):
            print("Warning: this dataset has no drug/disease names!")
        path = drrs_dataset_path+model_name+"s/"
        fnames = {"S": "DrugSim.txt", "P": "DiseaseSim.txt", "A": "DiDrA.txt"}
        A, P, S = [pd.read_csv(path+fnames[k], sep="\t", header=None) for k in ["A", "P", "S"]]
        A = A.T
        if (model_name in ["Cdataset", "Fdataset"]):
            model_name2 = "Cdataset_Aonly" if (model_name=="Cdataset") else "Gottlieb"
            A_dataset = load_dataset(model_name2, save_folder=save_folder)["ratings_mat"]
            drug_names = list(A_dataset.index)
            disease_names = list(A_dataset.columns)
            A.index = drug_names
            S.index = drug_names
            S.columns = drug_names
            A.columns = disease_names
            P.index = disease_names
            P.columns = disease_names
    assert A.shape[0] == S.shape[1]
    assert A.shape[1] == P.shape[1]
    return {"ratings_mat": A.fillna(0).astype(int), "users": P.astype(float), "items": S.astype(float)}

def merge_ratings(rating_dfs, user_col, item_col, rating_col, verbose=True):
    '''
        @param\trating_dfs\tPython list of Pandas dataframes: Pandas dataframes size #items x #users containing values {-1,0,1} to merge
        @param\tuser_col\tPython character string: the name of the column associated with users.
        @param\titem_col\tPython character string: the name of the column associated with items.
        @param\trating_col\tPython character string: the name of the column associated with ratings.
                We correct inconsistent outcomes as follows:
                - If there is at least one negative outcome (-1) reported, then it is a negative outcome (-1)
                 (in order to be conservative with respect to drug recommendations).
                - If there is at least one positive outcome (1) and no negative outcome (-1) reported, then it is a positive outcome (1).
        @return\tres_df\tPandas Dataframe: merged Pandas Dataframe of ratings containing values {-1,0,1}
    '''
    ratings = pd.concat(tuple([df.fillna(0)[[user_col, item_col, rating_col]] for df in rating_dfs]), axis=0)
    ratings.index = ["--".join(list(ratings.iloc[idx][[user_col,item_col]].astype(str))) for idx in range(len(ratings.index))]
    ratings[rating_col] = ratings[[rating_col]].groupby(level=0).apply(lambda x : ((-1)**((x==-1).any().any()))*int((x!=0).any().any())).astype(int)
    ratings.index = range(ratings.shape[0])
    return ratings

def ratings2matrix(ratings, user_col, item_col, rating_col):
    '''
        @param\tratings\tPandas Dataframe: 3 columns (user_col x item_col x rating_col) that lists all non zero values {-1,1}.
        @param\tuser_col\tPython character string: how to name the column associated with users.
        @param\titem_col\tPython character string: how to name the column associated with items.
        @param\trating_col\tPython character string: how to name the column associated with ratings.
        @return\tres_df\tPandas Dataframe: size #items x #users containing values {-1,0,1}.
    '''
    res_df = ratings.pivot_table(index=item_col, columns=user_col, values=rating_col).fillna(0).astype(int)
    return res_df