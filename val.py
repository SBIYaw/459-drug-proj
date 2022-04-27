# NOTE, needs drug_overlap_repodb.csv, drug_overlap_broad.csv to run

import pandas as pd
import numpy as np
import sys

def changeHeader(df):
    new_header = df.iloc[0] 
    df = df[1:]
    df.columns = new_header
    return df

def getData():
    ddr = pd.read_csv('validation data/DDR_full_net.csv',header=None)
    ddr = changeHeader(ddr)

    repodb = pd.read_csv('validation data/drug_overlap_repodb.csv',header=None)
    repodb = changeHeader(repodb)

    broad = pd.read_csv('validation data/drug_overlap_broad.csv',header=None)
    broad = changeHeader(broad)
    
    return ddr, repodb, broad

# Validations Methods (to be called by the calc_relative_risk script)
def validate(prim_drug: str, candidate_drug: str):
    ddr, repodb, broad = getData()
    # the 1st approach - checking disease overlaps
    cand_approved_diseases_repodb = checkDiseaseOverlap_repodb(prim_drug, candidate_drug, repodb)
    cand_approved_diseases_broad = checkDiseaseOverlap_broad(prim_drug, candidate_drug, broad)

    return cand_approved_diseases_repodb + cand_approved_diseases_broad

def checkDiseaseOverlap_repodb(prim_drug: str, candidate_drug: str, repodb):
    # getting the diseases assoc with the primary drug and the candidate drug
    prim_drug_diseases = repodb[repodb['drug'] == prim_drug]['ind_name'].to_numpy()
    prim_drug_diseases = list(map(lambda x: x.lower(), prim_drug_diseases))

    cand_drug_status = repodb[repodb['drug'] == candidate_drug]['status'].to_numpy()
    cand_drug_ind = repodb[repodb['drug'] == candidate_drug]['ind_name'].to_numpy()
    cand_drug_diseases = [(cand_drug_ind[i], cand_drug_status[i]) for i in range(0, len(cand_drug_ind))]

    # check whether the candidate drug can be used
    approved_for_diseases = []
    for i in cand_drug_diseases:
        if i[1] == "Approved" and i[0].lower() in prim_drug_diseases:
            approved_for_diseases.append(i[0])
    
    '''
     returning diseases the candidate drug is approved for, which are also
     the diseases assoc with the primary drug.
     
     If empty is list is returned, the candidate drug is a FP
    '''
    return approved_for_diseases

def checkDiseaseOverlap_broad(prim_drug: str, candidate_drug: str, broad):
    # getting the diseases assoc with the primary drug and the candidate drug
    prim_drug_diseases = np.asarray(broad[broad['drug'] == prim_drug]['disease_area'])[0].split('|')
    prim_drug_diseases = list(map(lambda x: x.lower(), prim_drug_diseases))
    
    cand_drug_status = broad[broad['drug'] == candidate_drug]['clinical_phase'].to_numpy()
    
    cand_drug_ind = np.asarray(broad[broad['drug'] == candidate_drug]['disease_area'])[0]
    if (cand_drug_ind != cand_drug_ind):
        print('candidate drug ' + candidate_drug + ' has no disease information')
        return []
    cand_drug_ind = cand_drug_ind.split('|')
    cand_drug_ind = list(map(lambda x: x.lower(), cand_drug_ind))
    cand_drug_status = str(cand_drug_status[0])

    # check whether the candidate drug can be used
    approved_for_diseases = []
    if cand_drug_status == 'Launched' or cand_drug_status == 'Phase 3':
        for i in cand_drug_ind:
            if i in prim_drug_diseases:
                approved_for_diseases.append(i)
    
    '''
     returning diseases the candidate drug is approved for, which are also
     the diseases assoc with the primary drug.

     If empty is list is returned, the candidate drug is a FP
    '''
    return approved_for_diseases