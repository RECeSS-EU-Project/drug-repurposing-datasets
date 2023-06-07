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

## From https://github.com/RECeSS-EU-Project/drug-repurposing-datasets

###########################
## VERBOSE               ##
###########################

def print_dataset(ratings, user_col, item_col, rating_col):
    '''
    Prints values of a drug repurposing dataset

    ...

    Parameters
    ----------
    ratings : pandas.DataFrame of shape (n_ratings, 3)
        the list of ratings with columns user_col, item_col, rating_col
    user_col : str
        column denoting users
    item_col : str
        column denoting items
    rating_col : str
        column denoting ratings in {-1, 0, 1}

    Returns
    -------
    None

    Prints
    -------
    The number of items/drugs, users/diseases, and the number of positive (1), negative (-1) and unknown (0) matchings.
    '''
    assert ratings.shape[1]==3
    assert all([c in ratings.columns for c in [user_col, item_col, rating_col]])
    assert all([a in [-1,0,1] for a in np.unique(ratings[[rating_col]].values)])
    ratings2 = ratings.copy()
    ratings2[item_col] = ratings2[item_col].astype(str)
    args = [len(np.unique(ratings2[item_col])), len(np.unique(ratings2[user_col]))]
    args += [ratings2.loc[ratings2[rating_col]==v].shape[0] for v in [1,-1,0]]
    assert args[2]+args[3]+args[4]==ratings2.shape[0]
    args[-1] = args[-1] if (args[-1]>0) else args[0]*args[1]-args[2]-args[3]
    dataset_str = "Ratings: %d drugs\t%d diseases\n%d positive, %d negative, %d unknown matchings"
    print(dataset_str % tuple(args))
    
def compute_sparsity(df):
    '''
    Computes the sparsity number of a collaborative filtering dataset

    ...

    Parameters
    ----------
    df : pandas.DataFrame of shape (n_items, n_users)
        the matrix of ratings where unknown matchings are denoted with 0

    Returns
    -------
    sparsity : float
        the percentage of non missing values in the matrix of ratings
    '''
    return 100*(np.sum(df.values!=0)/(df.shape[0]*df.shape[1]))

#############################
## HANDLING RATINGS MATRIX ##
#############################

def matrix2ratings(df, user_col="user", item_col="item", rating_col="rating"):
    '''
    Converts a matrix into a list of ratings

    ...

    Parameters
    ----------
    df : pandas.DataFrame of shape (n_items, n_users)
        the matrix of ratings in {-1, 1, 0} where unknown matchings are denoted with 0
    user_col : str
        column denoting users
    item_col : str
        column denoting items
    rating_col : str
        column denoting ratings in {-1, 0, 1}

    Returns
    -------
    ratings : pandas.DataFrame of shape (n_ratings, 3)
        the list of known ratings where the first column correspond to users, second to items, third to ratings
    '''
    assert all([a in [-1,0,1] for a in np.unique(df.values)])
    non_missing = np.argwhere(df.values!=0)
    res_df = pd.DataFrame([], index=range(non_missing.shape[0]))
    res_df[user_col] = [df.columns[x] for x in list(non_missing[:, 1].flatten())]
    res_df[item_col] = [df.index[x] for x in list(non_missing[:, 0].flatten())]
    res_df[rating_col] = [df.values[i,j] for i,j in non_missing.tolist()]
    return res_df[[user_col,item_col,rating_col]]

def ratings2matrix(ratings, user_col, item_col, rating_col):
    '''
    Converts a list of ratings into a matrix

    ...

    Parameters
    ----------
    ratings : pandas.DataFrame of shape (n_ratings, 3)
        the list of known ratings where the first column (user_col) correspond to users, second (item_col) to items, third (rating_col) to ratings in {-1,0,1}
    user_col : str
        column denoting users
    item_col : str
        column denoting items
    rating_col : str
        column denoting ratings in {-1, 0, 1}

    Returns
    -------
    df : pandas.DataFrame of shape (n_items, n_users)
        the matrix of ratings in {-1, 1, 0} where unknown matchings are denoted with 0
    '''
    assert ratings.shape[1]==3
    assert all([c in ratings.columns for c in [item_col, user_col, rating_col]])
    assert all([a in [-1,0,1] for a in np.unique(ratings[[rating_col]].values)])
    res_df = ratings.pivot_table(index=item_col, columns=user_col, values=rating_col).fillna(0).astype(int)
    return res_df

#############################
## LOADING DATASETS        ##
#############################

def load_dataset(model_name, save_folder="./", sep_feature="-"):
    '''
    Loads a drug repurposing dataset

    ...

    Parameters
    ----------
    model_name : str
        the name of the dataset to load. Should belong to the following list: ["Gottlieb", "DNdataset", "Cdataset", "LRSSL", "PREDICT_Gottlieb", "TRANSCRIPT", "PREDICT"]
    save_folder : str
        the path to the folder where dataset-related files are or will be stored

    Returns
    -------
    dataset_di : dictionary
        a dictionary where key "ratings_mat" contains the drug-disease matching pandas.DataFrame of shape (n_drugs, n_diseases) (where missing values are denoted by 0), key "users" correspond to the disease pandas.DataFrame of shape (n_disease_features, n_diseases), and "items" correspond to the drug feature pandas.DataFrame of shape (n_drug_features, n_drugs)
    '''
    assert model_name in ["Gottlieb", "Cdataset_Aonly", "indep", "Fdataset", "DNdataset", "Cdataset", "TRANSCRIPT", "PREDICT", "LRSSL", "LRSSL2", "PREDICT_Gottlieb"]
    if (model_name == "LRSSL"):
        url_lrssl = "https://raw.githubusercontent.com/LiangXujun/LRSSL/master/"
        lrssl_dataset_path = save_folder+"LRSSL/"
        fnames = {
            "A": "drug_dis_mat.txt", "P": "disease_similarity.txt", "S1": "drug_from_drugbank_without_ind_dommat.txt",
            "S2": "drug_pubchem_mat.txt", "S3": "drug_target_domain_mat.txt", "N4": "drug_target_go_mat.txt", 
            "N5": "drug_without_ind_dommat_new.txt", "S6": "drug_without_ind_gomat_new.txt", 
            "N7": "drug_without_ind_pubchem_mat_new.txt",
        }
        if (not os.path.exists(lrssl_dataset_path+fnames["A"])):
            subprocess.call(" ".join(['mkdir', '-p', lrssl_dataset_path]), shell=True)
            for fn in fnames:
                subprocess.call(" ".join(["wget", "-qO", lrssl_dataset_path+fnames[fn], url_lrssl+fnames[fn]]), shell=True)
        A = pd.read_csv(lrssl_dataset_path+fnames["A"], index_col=0, sep="\t")
        P = pd.read_csv(lrssl_dataset_path+fnames["P"], index_col=0, sep="\t")
        A.columns = P.columns
        S1 = pd.read_csv(lrssl_dataset_path+fnames["S2"], index_col=0, sep="\t")
        S2 = pd.read_csv(lrssl_dataset_path+fnames["S3"], index_col=0, sep="\t")
        S = S1.join(S2, how="outer").T[A.index]
    if (model_name in ["LRSSL2", "PREDICT_Gottlieb", "OrphanDrug"]):
        print("Warning: this dataset has no drug/disease names!")
        url_dda_skf = "https://github.com/GCQ2119216031/DDA-SKF/raw/master/data/"
        dda_skf_dataset_path = save_folder+"DDA_SKF/data/"
        mmodel_name = model_name if (model_name=="OrphanDrug") else ("PREDICT" if (model_name=="PREDICT_Gottlieb") else "LRSSL")
        if (not os.path.exists(dda_skf_dataset_path+mmodel_name+".mat")):
            subprocess.call(" ".join(['mkdir', '-p', dda_skf_dataset_path]), shell=True)
            subprocess.call(" ".join(["wget", "-qO", dda_skf_dataset_path+mmodel_name+".mat", url_dda_skf+mmodel_name+".mat"]), shell=True)
        ddt = scipy.io.loadmat(dda_skf_dataset_path+mmodel_name+".mat")
        if (mmodel_name=="LRSSL"):
            A = ddt["lrssladmatdgc"]
            S_chemical = ddt["lrsslsimmatdcchemical"]
            S_sideeffects = ddt["lrsslsimmatdcgo"]
            S = pd.DataFrame(np.concatenate((S_chemical,S_sideeffects), axis=0), index=range(S_chemical.shape[0]+S_sideeffects.shape[0]), columns=["drug%d" % (i+1) for i in range(S_chemical.shape[1])])
            S.index = [("chemical" if (iss < S_chemical.shape[0]) else "se")+sep_feature+("drug%d" % s+1) for iss, s in enumerate(list(S.index))]
            P = ddt["lrsslsimmatdg"]
            P = pd.DataFrame(P, index=["disease%d" % (i+1) for i in range(P.shape[0])], columns=["disease%d" % (i+1) for i in range(P.shape[1])])
        elif (mmodel_name=="PREDICT"):
            A = ddt['predictAdMatdgc'].T
            S_chemical = ddt['predictSimMatdcChemical']
            S_domain = ddt['predictSimMatdcDomain']
            S_GO = ddt['predictSimMatdcGo']
            S = pd.DataFrame(np.concatenate((S_chemical,S_domain,S_GO), axis=0), index=range(S_chemical.shape[0]+S_domain.shape[0]+S_GO.shape[0]), columns=["drug%d" % (i+1) for i in range(S_chemical.shape[1])])
            S.index = [("chemical" if (iss < S_chemical.shape[0]) else ("domain" if (iss < S_chemical.shape[0]+S_domain.shape[0]) else "go"))+sep_feature+("drug%d" % (s+1)) for iss, s in enumerate(list(S.index))]
            P = ddt['predictSimMatdg']
            P = pd.DataFrame(P, index=["disease%d" % (i+1) for i in range(P.shape[0])], columns=["disease%d" % (i+1) for i in range(P.shape[1])])
        else:
            raise ValueError("Undefined dataset '%s'" % mmodel_name)
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
                P.index = [("phenotype" if (iss < P_phenotype.shape[0]) else "semantic")+sep_feature+str(s) for iss, s in enumerate(list(P.index))]
                S = S_se.T.join(S_signature.T, how="outer").T
                S.index = [("se" if (iss < S_se.shape[0]) else "signature")+sep_feature+str(s) for iss, s in enumerate(list(S.index))]
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
    A = A.fillna(0).astype(int)
    P = P.astype(float)
    S = S.astype(float)
    assert A.shape[0] == S.shape[1]
    assert A.shape[1] == P.shape[1]
    assert all([a in [-1, 0, 1] for a in np.unique(A).tolist()])
    return {"ratings_mat": A, "users": P, "items": S}

#############################
## MERGING RATINGS         ##
#############################

def merge_ratings(rating_dfs, user_col, item_col, rating_col):
    '''
    Merges rating lists from several sources by solving conflicts. Conflicting ratings are resolved as follows: if there is at least one negative rating (-1) reported for a (drug, disease) pair, then the final rating is negative (-1); if there is at least one positive rating (1) and no negative rating (-1) reported, then the final rating is positive (1)

    ...

    Parameters
    ----------
    rating_dfs : list of pandas.DataFrame of shape (n_ratings, 3)
        the list of rating lists where one column (of name user_col) is associated with users, one column (of name item_col) is associated with items, and one column (of name rating_col) is associated with ratings in {-1, 0, 1}
    user_col : str
        column denoting users
    item_col : str
        column denoting items
    rating_col : str
        column denoting ratings in {-1, 0, 1}
    verbose : bool
       
    Returns
    -------
    rating_df : pandas.DataFrame of shape (n_ratings, 3)
        the list of rating lists where one column (of name user_col) is associated with users, one column (of name item_col) is associated with items, and one column (of name rating_col) is associated with ratings in {-1, 0, 1}
    '''
    for ratings in ratings_dfs:
        assert ratings.shape[1]==3
        assert all([c in ratings.columns for c in [item_col, user_col, rating_col]])
        assert all([a in [-1,0,1] for a in np.unique(ratings[[rating_col]].values)])
    ratings = pd.concat(tuple([df.fillna(0)[[user_col, item_col, rating_col]] for df in rating_dfs]), axis=0)
    ratings.index = ["--".join(list(ratings.iloc[idx][[user_col,item_col]].astype(str))) for idx in range(len(ratings.index))]
    ratings[rating_col] = ratings[[rating_col]].groupby(level=0).apply(lambda x : ((-1)**((x==-1).any().any()))*int((x!=0).any().any())).astype(int)
    ratings.index = range(ratings.shape[0])
    return ratings