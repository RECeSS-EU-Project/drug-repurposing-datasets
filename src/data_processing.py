#coding: utf-8

############################
## BioMaRt                ##
############################

import subprocess as sb

def get_biomart(pr, gene):
    '''
        @param\tpr\tPython integer: number (can be ignored)
        @param\tgene\tPython character string: gene to convert
        @return\tentrez_id\tPython integer: gene EntrezID
    '''
    cmd='wget -q -O - \'http://www.ensembl.org/biomart/martservice?query=<?xml'
    cmd+=' version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName'
    cmd+=' = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" '
    cmd+='datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" '
    cmd+='interface = "default" ><Filter name = "hgnc_symbol" value = "'+gene+'"/>'
    cmd+='<Attribute name = "hgnc_symbol" /><Attribute name = "entrezgene_id"/>'
    cmd+='</Dataset></Query>\''
    try:
        vals = sb.check_output(cmd, shell=True).decode("utf-8")
    except:
        return "-"
    vals = vals.split("\n")[:-1]
    vals = [x for x in vals if (len(x.split("\t")[-1])>0)]
    if (len(vals)==0):
        return "-"
    val = vals[0].split("\t")[-1]
    if (len(val)==0):
        return "-"
    return int(val)

############################
## LINCS L1000            ##
############################

import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from tqdm import tqdm
import pickle

from NORDic.UTILS.LINCS_utils import *

def pubchem_cid2pert_iname(cid, user_key):
    '''
        @param\tcid\tPython character string: PubChem CID
        @param\tuser_key\tPython character string: LINCS L1000 clue.io API key
        @return\tiname\tPython character string: corresponding pert_iname identifier in LINCS L1000
    '''
    endpoint = "perts"
    method = "filter"
    where = {"pert_type": "trt_cp", "pubchem_cid": int(cid)}
    params = {
            "where": where,
            "fields": ["pert_iname"]
    }
    request_url = build_url(endpoint, method, params=params, user_key=user_key)
    data_pert_ = post_request(request_url, quiet=True)
    if (len(data_pert_)>0):
        pert_iname = data_pert_[0]["pert_iname"]
        return pert_iname
    return None

def pubchem_cid2brds(cid, user_key):
    '''
        @param\tcid\tPython character string: PubChem CID
        @param\tuser_key\tPython character string: LINCS L1000 clue.io API key
        @return\tlst\tPython list of character strings: corresponding BRD identifiers
    '''
    endpoint = "perts"
    method = "filter"
    pert_iname = pubchem_cid2pert_iname(cid, user_key)
    if (pert_iname is None):
        return []
    params = {
        "where": {"pert_type": "trt_cp", "pert_iname": pert_iname},
        "fields": ["pert_id"]
    }
    request_url = build_url(endpoint, method, params=params, user_key=user_key)
    data_pert_ = post_request(request_url, quiet=True)
    if (len(data_pert_)==0):
        return []
    return list(set([dt["pert_id"] for dt in data_pert_]))

def brd2bestcellline(brd, user_key):
    '''
        @param\tbrd\tPython character string: BRD drug identifier
        @param\tuser_key\tPython character string: LINCS L1000 clue.io API key
        @return\tcell_line\tPython character string: cell id which corresponds 
                to the cell line which maximizes the Transcriptional Activity 
                Score (TAS) for the drug associated with the BRD identifier
    '''
    endpoint = "perts"
    method = "filter"
    where = {"pert_type": "trt_cp", "pert_id": brd}
    params = {
            "where": where,
            "fields": ["cell_tas"]
    }
    request_url = build_url(endpoint, method, params=params, user_key=user_key)
    data_pert_ = post_request(request_url, quiet=True)
    assert len(data_pert_) > 0
    data_pert = pd.DataFrame(data_pert_[0]['cell_tas']).sort_values('tas')
    cell_line = list(data_pert["cell_id"])[-1]
    return cell_line

def get_profiles(drug_ids, gene_list, path_to_lincs, user_key, save_folder="./", nsigs=2, njobs=1, verbose=0):
    '''
        Check out https://clue.io/developer-resources#apisection
        @param\tdrug_list\tPython list of integers: PubChemID identifiers
        @param\tgene_list\tPython list of integers: Entrez ID identifiers
        @param\tpath_to_lincs\tPython character string: path where the LINCS L1000 files are stored
        @param\tuser_key\tPython character string: LINCS L1000 clue.io API key
        @param\tsave_folder\tPython character string: path where the intermediary files are stored
        @param\tnsigs\tPython integer[default=2]: minimal number of replicates
        @param\tnjobs\tPython integer[default=1]: use parallelism if >1
        @param\tverbose\tPython integer[default=0]: prints
        @returns\tsignatures\tPandas DataFrame: Characteristic Direction-binarized drug signatures
    '''

    def single_run(di, drug):
        print("%d/%d"%(di+1,len(drug_ids)))
        if (not os.path.exists(save_folder+"sig_ids_%d.pck" % drug)):
            sig_ids = {}
        else:
            with open(save_folder+"sig_ids_%d.pck" % drug, "rb") as f:
                sig_ids = pickle.load(f)
            return sig_ids
        brds = pubchem_cid2brds(drug, user_key)
        endpoint = "sigs"
        method = "filter"
        ############# GET TREATED PROFILE
        ## Select perturbation type, drug
        data_pert_ = []
        msg = "%s (pubchem %s) not found." % (drug,drug)
        if (len(brds)==0):
            print(msg+" (0, BRDs)")
            with open(save_folder+"sig_ids_%d.pck" % drug, "wb") as f:
                pickle.dump(None, f)
            return None
        for brd in brds:
            cell_line = brd2bestcellline(brd, user_key)
            where = {"pert_type": "trt_cp", "pert_id": brd, "cell_id": cell_line}
            ## Return distil id (unique to profile), brew prefix (unique to set of replicates)
            params = {
                    "where": where,
                    "fields": ["distil_id", "brew_prefix", "distil_cc_q75", 
                               "pct_self_rank_q25", "distil_ss", "cell_id"]
            }
            ## Create URL
            request_url = build_url(endpoint, method, params=params, user_key=user_key)
            ## Get results from request
            data_pert_ += post_request(request_url, quiet=True)
        if (len(data_pert_) == 0):
            print(msg+" (0)")
            with open(save_folder+"sig_ids_%d.pck" % drug, "wb") as f:
                pickle.dump(None, f)
            return None
        ## "Gold" profiles: https://clue.io/connectopedia/glossary
        data_pert = [d for d in data_pert_ if ((d["distil_cc_q75"]>=0.2) and (d["pct_self_rank_q25"]<=0.05))]  
        if (len(data_pert_) == 0):
            print(msg+" (1)")
            with open(save_folder+"sig_ids_%d.pck" % drug, "wb") as f:
                pickle.dump(None, f)
            return None
        ## Remove profiles with less than nsigs replicates
        data_pert = [d for d in data_pert_ if (len(d["distil_id"])>=nsigs)]
        if (len(data_pert) == 0):
            print(msg+" (2)")
            with open(save_folder+"sig_ids_%d.pck" % drug, "wb") as f:
                pickle.dump(None, f)
            return None
        ## Get signature with maximal value of selection measure 
        dt_pert = data_pert[np.argmax([d["distil_ss"] for d in data_pert])]
        ############# GET SAME-PLATE CONTROL (VEHICLE) PROFILE
        params = {
                "where": {"pert_type": "ctl_vehicle", "brew_prefix": dt_pert["brew_prefix"][0]},
                "fields": ["distil_id", "distil_ss", "cell_id"]
        }
        ## Create URL
        request_url = build_url(endpoint, method, params=params, user_key=user_key)
        ## Get results from request
        data_ctrl_ = post_request(request_url, quiet=(not verbose))
        ## Remove profiles with less than nsigs replicates
        data_ctrl = [d for d in data_ctrl_ if (len(d["distil_id"])>=nsigs)]
        ## Get signature with maximal value of selection measure
        dt_ctrl = data_ctrl[np.argmax([d["distil_ss"] for d in data_ctrl])]
        assert dt_ctrl["cell_id"]==dt_pert["cell_id"]
        ############# GET PROFILES FROM LINCS
        sig_ids = {drug : {"cell_line": dt_ctrl["cell_id"], "pert": dt_pert["distil_id"], "ctrl": dt_ctrl["distil_id"]}}
        with open(save_folder+"sig_ids_%d.pck" % drug, "wb") as f:
            pickle.dump(sig_ids, f)
        return sig_ids

    if (njobs>1):
        sigs_ids_list = Parallel(n_jobs=njobs, backend='loky')(delayed(single_run)(di, drug) for di, drug in tqdm(enumerate(drug_ids)))
    else:
        sigs_ids_list = [single_run(di, drug) for di, drug in tqdm(enumerate(drug_ids))]   

    sigs_list, seen_drugs = [], []
    for si, sig_ids in enumerate(sigs_ids_list):
        print("* Signature %d/%d" % (si+1,len(sigs_ids_list)))
        if (sig_ids is None):
            continue
        drug = list(sig_ids.keys())[0]
        if (drug in seen_drugs):
            continue
        seen_drugs.append(drug)
        msg = "%s (pubchem %s) not found." % (drug,drug)
        cell_line = sig_ids[drug]["cell_line"]
        pert, ctrl = sig_ids[drug]["pert"], sig_ids[drug]["ctrl"]
        sig_ids = pert+ctrl
        profiles = create_restricted_drug_signatures(sig_ids, gene_list, path_to_lincs, which_lvl=[3], strict=False)
        if (profiles is None):
            print(msg+" (3)")
            continue
        brew_suffix = [s.split("_")[-1] for s in profiles.columns]
        pert_brew_suffix = [s.split("_")[-1] for s in pert]
        ctrl_brew_suffix = [s.split("_")[-1] for s in ctrl]
        samples = [2 if (s in pert_brew_suffix) else (1 if (s in ctrl_brew_suffix) else 0) for s in brew_suffix]
        assert all([s in list(set(samples)) for s in [1,2]])
        assert all([s>0 for s in samples])
        ############# TURN PROFILES INTO SIGNATURES ("DIFFERENTIAL ANALYSIS")
        sigs = binarize_via_CD(profiles, samples, binarize=0, nperm=10000, quiet=(verbose==0))
        sigs.columns = [drug]
        sleep(1)
        sigs_list.append(sigs) 
        print("Total %d signatures" % len(sigs_list))
    return sigs_list[0].join(sigs_list[1:]).fillna(0)

############################
## CHEMICAL SIMILARITY    ##
############################

## Compute hashed fingerprints & similarity scores
## https://ctr.fandom.com/wiki/Report_the_similarity_between_two_structures
from indigo import Indigo, IndigoException

import pubchempy as pcp

def smiles_similarity(s1, s2):
    '''
        @param\ts1,s2\ttwo SMILES-encoded fingerprints
        @return\tsim\tTanimoto score similarity value
    '''
    indigo_module = Indigo()
    try:
        m1 = indigo_module.loadMolecule(s1)
        m2 = indigo_module.loadMolecule(s2)
    except IndigoException as err:
        if ("chirality not possible on atom" in str(err)):
            return np.nan
        elif ("does not match pending bond description" in str(err)):
            return np.nan
        raise
    # Aromatize molecules because second molecule is not in aromatic form
    m1.aromatize()
    m2.aromatize()
    # Calculate similarity between "similarity" fingerprints
    fp1 = m1.fingerprint("sim")
    fp2 = m2.fingerprint("sim")
    sim = indigo_module.similarity(fp1, fp2, "tanimoto")
    return sim

def get_pubchem_smiles(pubchem_id):
    '''
        @param\tpubchem_id\ta PubChem CID
        @return\tsmiles\tSMILES fingerprint from PubChem
    '''
    cpd = pcp.Compound.from_cid(int(pubchem_id))
    return cpd.canonical_smiles

############################
## SIDE EFFECT SIMILARITY ##
############################

from selenium import webdriver
from selenium.webdriver import Firefox
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.common.keys import Keys

from time import sleep

from sklearn.metrics import pairwise_distances
from Bio import Align

def get_drugbank_ids(drugnames, driver_folder="./", verbose=True):
    '''
        @param\tdrugnames\tPython list of character strings: SIDER provided drug names
        @param\tdriver_folder\tPython character string: where to save the driver for browsing
        @param\tverbose\tPython boolean: print out the result
        @return\tdrugbank_cids\tPython list of integers: corresponding DrugBank CIDs
    '''
    if (not os.path.exists(driver_folder+"geckodriver")):
        driver_url="https://github.com/mozilla/geckodriver/releases/download/v0.30.0/geckodriver-v0.30.0-linux64.tar.gz"
        sb.Popen(["wget", "-qO-", driver_url+" | tar -xvz -C "+driver_folder])
        sleep(3)
        os.environ["PATH"] += os.pathsep + driver_folder

    options = Options()
    prefs = {
        "browser.download.folderList": 2,
        "browser.download.useDownloadDir": True,
        "browser.download.manager.showWhenStarting": False,
        "browser.download.dir": driver_folder,
        "browser.helperApps.neverAsk.saveToDisk": "text/csv",
    }
    for opt in prefs:
        options.set_preference(opt, prefs[opt])
    options.headless = True
    driver = webdriver.Firefox(executable_path=driver_folder+"geckodriver", options=options)
    drugbank_cids = []

    missing_lst = []

    for idr, drugname in enumerate(drugnames):
        driver.get("http://en.wikipedia.org/wiki/Special:Search?search="+drugname)
        driver.implicitly_wait(10)
        sleep(1)
        if (verbose):
            print((drugname, driver.current_url))

        try:
            drugbank=driver.find_element(By.XPATH,"//span[@title='www.drugbank.ca']").click()
            driver.implicitly_wait(10)
            sleep(1)
            drugbank_cid = driver.current_url.split("/")[-1]
        except:
            missing_cids = {
		'diflorasone':"DB00223", 'APAs': None,
		'mTHPC': None, 'Locorten': None,
		'polythiazide': None, 
            }
            missing_lst.append(drugname)
            drugbank_cid = None#missing_cids[drugname]

        drugbank_cids.append(drugbank_cid)
        if (verbose):
            print("%d/%d %s %s" % (idr+1, len(drugnames), drugname, str(drugbank_cid)))
    driver.close()
    print(missing_lst)
    return drugbank_cids

def get_pubchemcid_from_SIDER(drugnames, driver_folder="./", verbose=False):
    '''
        @param\tdrugnames\tPython list of character strings: SIDER provided drug names
        @param\tdriver_folder\tPython character string: where to save the driver for browsing
        @param\tverbose\tPython boolean: print out the result
        @return\tpubchem_cids\tPython list of integers: corresponding PubChem CIDs
    '''
    url_sider = "http://sideeffects.embl.de/drugs/"

    if (not os.path.exists(driver_folder+"geckodriver")):
        driver_url="https://github.com/mozilla/geckodriver/releases/download/v0.30.0/geckodriver-v0.30.0-linux64.tar.gz"
        sb.Popen(["wget", "-qO-", driver_url+" | tar -xvz -C "+driver_folder])
        sleep(3)
        os.environ["PATH"] += os.pathsep + driver_folder

    options = Options()
    prefs = {
        "browser.download.folderList": 2,
        "browser.download.useDownloadDir": True,
        "browser.download.manager.showWhenStarting": False,
        "browser.download.dir": driver_folder,
        "browser.helperApps.neverAsk.saveToDisk": "text/csv",
    }
    for opt in prefs:
        options.set_preference(opt, prefs[opt])
    options.headless = True
    driver = webdriver.Firefox(executable_path=driver_folder+"geckodriver", options=options)
    driver.get(url_sider)
    driver.implicitly_wait(10)
    sleep(1)

    pubchem_cids = []

    for idr, drugname in enumerate(drugnames):
        drug_form = driver.find_element(By.XPATH, '//*[@id="searchBox"]')
        drug_form.clear()
        drug_form.send_keys(drugname)
        driver.implicitly_wait(10)

        submit_button = driver.find_element(By.XPATH, "/html/body/nav/div/div[2]/form/div/div/button").click()
        driver.implicitly_wait(10)
        sleep(1)

        if (verbose):
            print((drugname, driver.current_url))
        try:
            pubchem_cid = int(driver.current_url.split("/")[-2])
        except:
            try:
                first_res = driver.find_element(By.XPATH, "/html/body/div[2]/ul/li[1]/a").click()
                driver.implicitly_wait(10)
                sleep(1)
                if (verbose):
                    print((2,drugname, driver.current_url))
                pubchem_cid = int(driver.current_url.split("/")[-2])
            except:
                try:
                    first_res = driver.find_element(By.XPATH, "/html/body/div[3]/div[2]/table/tbody/tr/td[1]/ul/li[1]/a").click()
                    driver.implicitly_wait(10)
                    sleep(1)
                    if (verbose):
                        print((3,drugname, driver.current_url))
                    pubchem_cid = int(driver.current_url.split("/")[-2])
                except:
                    pubchem_cid = None
                    driver.get(url_sider)
        pubchem_cids.append(pubchem_cid)

        if (verbose):
            print("%d/%d %s %s" % (idr+1, len(drugnames), drugname, str(pubchem_cid)))

    driver.close()
    return pubchem_cids

def sideeffect_similarity(se):
    '''
        @param\tse\tPandas DataFrame: index=drugs, column=side effect SIDER identifier
        @return\tse_based\tJaccard score similarity matrix
    '''
    se = se.groupby(level=0).apply(lambda x : ",".join(list(sorted(set(list(x.values.flatten()))))))
    ## one-hot encoding of drugs according to their side effects
    vals = list(set([y for x in list(se.values) for y in x.split(",")]))
    se = se.T.apply(lambda x : [int(v in x) for v in vals])
    se = pd.DataFrame(np.matrix([np.array(x) for x in se.values], dtype=bool), index=se.index, columns=vals)
    se_based = pd.DataFrame(1-pairwise_distances(se.values, metric='jaccard'), index=se.index, columns=se.index)
    return se_based

############################
## TARGET SEQ SIMILARITY  ##
############################

from joblib import Parallel, delayed

def seq_function(i, j, target_ids, sequences_di):
    '''
        @param\ti,j\tPython integers: indices in target_ids
        @param\ttarget_ids\tPython list of character strings: drug or disease identifiers
        @param\tsequences_di\tPython dictionary: sequences indexed by corresponding drug or disease identifiers
        @return\tscore\tPython float: average alignment score between sequences associated with targets_ids[i] and target_ids[j] (=0 if i>j)
    '''
    sequence_aligner = Align.PairwiseAligner()
    if (i>j):
        return 0.
    target_i, target_j = target_ids[i], target_ids[j]
    ## all sequences in which drug_i appears
    target_seq_i = sequences_di[target_i]
    if (i==j):
        align_scores = [sequence_aligner.score(s1, s2) for s1 in target_seq_i for s2 in target_seq_i]
    else:
        ## all sequences in which drug_j appears
        target_seq_j = sequences_di[target_j]
        align_scores = [sequence_aligner.score(s1, s2) for s1 in target_seq_i for s2 in target_seq_j]
    ## take the average score across all sequences
    score = np.mean(align_scores) if (len(align_scores)>0) else np.nan
    return score

def sequence_similarity(target_ids, target_ids_per_seq, target_sequences, njobs=1, save_folder="./"):
    '''
        @param\ttarget_ids\tPython list of character strings: drug or disease identifiers
        @param\ttarget_ids_per_seq\tPython list of Python list of character strings: drug or disease identifiers associated with a sequence in target_sequences
        @param\ttarget_sequences\tPython list of sequences: list of sequences (same length as @target_ids_per_seq)
        @return\tsequence_similarity_df\tPandas DataFrame: alignment score matrix which contains the average alignment score between sequences associated with targets_ids[i] and target_ids[j]
    '''
    sequence_aligner = Align.PairwiseAligner()
    parallel = Parallel(n_jobs=njobs, backend='loky')
    ndrugs=len(drug_ids)
    sequences_di = {
	drug: [
            target_sequences[idx] for idx in [ix for ix, drug_ls in enumerate(drug_ids_per_seq) if (drug in drug_ls)]
            ] for drug in drug_ids
    }
    if (os.path.exists(save_folder+"sequence_similarity.out")):
        sequence_similarity_arr = np.loadtxt(save_folder+"sequence_similarity.out")
        start_ = [[i,j] for i,j in np.argwhere(sequence_similarity_arr==0).tolist() if (j>=i)]
        if (len(start_)==0):
                start_i = ndrugs
        else:
                start_i = start_[0][0]
    else:
        sequence_similarity_arr = np.zeros((ndrugs, ndrugs))
        start_i = 0
    for i in range(start_i, ndrugs):
        print("%d/%d" % (i+1, ndrugs))
        Sj = range(ndrugs)#range(max(i, start_j), ndrugs)
        if (njobs==1):
            sequence_similarity = [seq_function(i, j, drug_ids, sequences_di) for j in Sj]
        else:
            sequence_similarity = parallel(delayed(seq_function)(i, j, drug_ids, sequences_di) for j in Sj)
        sequence_similarity_arr[i,:] = sequence_similarity
        np.savetxt(save_folder+"sequence_similarity.out", sequence_similarity_arr)
    sequence_similarity = np.array(sequence_similarity_arr).reshape((ndrugs, ndrugs))
    sequence_similarity_df = pd.DataFrame(sequence_similarity_arr, index=drug_ids, columns=drug_ids)
    sequence_similarity_df_diag = [sequence_similarity_df.iloc[i,i] for i in range(ndrugs)]
    sequence_similarity_df = sequence_similarity_df.fillna(0)+sequence_similarity_df.fillna(0).T
    sequence_similarity_df.iloc[range(ndrugs), range(ndrugs)] = sequence_similarity_df_diag
    return sequence_similarity_df