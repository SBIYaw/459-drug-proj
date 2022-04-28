# NOTE, needs drug_overlap_repodb.csv, drug_overlap_broad.csv to run

import pandas as pd
import numpy as np
import sys

def changeHeader(df):
    new_header = df.iloc[0] 
    df = df[1:]
    df.columns = new_header
    return df

def get_DD_relations():
    ddr = pd.read_csv('validation data/DDR_full_net.csv',header=None)
    ddr = changeHeader(ddr)
    return ddr

def get_repoDB():
    repodb = pd.read_csv('validation data/drug_overlap_repodb.csv',header=None)
    repodb = changeHeader(repodb)
    return repodb

def get_broad():
    broad = pd.read_csv('validation data/drug_overlap_broad.csv',header=None)
    broad = changeHeader(broad)
    return broad

# Validations Methods (to be called by the calc_relative_risk script)
def validate(prim_drug: str, candidate_drug: str):
    repodb = get_repoDB()
    broad = get_broad()

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
    prim_drug_dis = np.asarray(broad[broad['drug'] == prim_drug]['disease_area'])[0].split('|')
    prim_drug_dis = list(map(lambda x: x.lower(), prim_drug_dis))
    
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
            if i in prim_drug_dis:
                approved_for_diseases.append(i)
    
    '''
     returning diseases the candidate drug is approved for, which are also
     the diseases assoc with the primary drug.

     If empty is list is returned, the candidate drug is a FP
    '''
    return approved_for_diseases

def checkNeighbors_repodb(prim_drug: str, ddr, repodb):
    candidate_drugs = []

    drug_df = ddr[ddr['drug_a'] == prim_drug]    

    prim_drug_diseases = repodb[repodb['drug'] == prim_drug]['ind_name'].to_numpy()
    prim_drug_diseases = list(map(lambda x: x.lower(), prim_drug_diseases))

    candidate_drug = [(x.lower(),y,str(np.asarray(z)).split('|')) for x,y,z in zip(drug_df['drug_b'],drug_df['clinical_phase_b'],drug_df['disease_area_b'])]
    for drug in candidate_drug:
        if drug[1] == 'Launched' or drug[1] == 'Phase 3':
            common_diseases = [i for i in drug[2] if i in prim_drug_diseases]
            if len(common_diseases) != 0:
                print('LGTM! ',drug)
                candidate_drugs.append((drug[0],common_diseases))
    
    return candidate_drugs

def checkNeigbors(prim_drug: str, potential_drugs: list):
    neighboring_repurposable_drugs = []
    ddr = get_DD_relations()
    repo_db = get_repoDB()

    neighboring_cand_drugs = checkNeighbors_repodb(prim_drug, ddr, repo_db)
    
    #check if candidate drugs are within the potential drugs (drugs which appear in the rank calculation)
    for n_drug in neighboring_cand_drugs:
        for p_drug in potential_drugs:
            if n_drug[0].lower() == p_drug.lower():
                print('catch!')
                neighboring_repurposable_drugs.append(n_drug[0], n_drug[1]) #returning the drug and its assoc diseases

    return neighboring_repurposable_drugs