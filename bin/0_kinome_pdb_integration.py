#!/usr/bin/env python

import pandas as pd
import numpy as np
import os,sys

def integrate_kinome_pdb(kinome_filename,kincon_filename):
    #read the two tables

    df_kinome = pd.read_csv(kinome_filename,sep = '\t')
    df_pdb = pd.read_csv(kincon_filename,sep='\t')
    df_pdb = df_pdb.drop(df_pdb.columns[[0,2,3,5]], axis=1)
    
    #creat a dictionary for {uniprot:pdb1 ,pdb2...} 
    dic_key = set()
    dic = {}
    for i in df_pdb['uniprot']:
        dic_key.add(i)
    
    dic = {}
    for index,row in df_pdb.iterrows():
        dic.setdefault(row['uniprot'],set()).add(row['PDB'])
    
    dic_detect = {} 
    for i in df_kinome['uniprot'] :
        if i in dic.keys():
            dic_detect[i]=dic[i]
        else:
            dic_detect.setdefault(i,[]).append('Nah')
    #put the dictionary to the kinome dataframe
    sLength = len(df_kinome['RANK'])
    df_kinome['pdb_detect'] = np.random.randn(sLength) #creat a empty list
    
    pdb_detect = []
    for index,row in df_kinome.iterrows():
        if row['uniprot'] in dic_detect.keys():
            row['pdb_detect'] = list(dic_detect[row['uniprot']])
        else:
            row['pdb_detect'] = 'Nah'
        pdb_detect.append(row['pdb_detect'])
    df_kinome['pdb_detect'] = list(pdb_detect)
    
    #convert it to str
    pdb_detect_str = []
    for index,row in df_kinome.iterrows():
        if type(row['pdb_detect']) == list:
            row= ','.join(row['pdb_detect'])
        else:
            pass
        pdb_detect_str.append(row)
    df_kinome['pdb_detect'] = pdb_detect_str
    
    df_kinome.to_csv('df_kinome.tsv',sep='\t')
    return df_kinome

if __name__ == "__main__":
    kinome_filename = sys.argv[1]
    kincon_filename = sys.argv[2]
    integrate_kinome_pdb(kinome_filename,kincon_filename)