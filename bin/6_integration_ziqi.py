#!/usr/bin/env python

import pandas as pd
import os,sys
import shutil



bin_folder = '.'
input_data_folder = '../inputdata/'
result_folder = '../results/'
refined_path = result_folder + 'refined_ouput/'

os.chdir(bin_folder)

all_results_path = result_folder +'all_results/'
contacts_path = result_folder +'contacts/'

os.chdir(contacts_path) #change working directory to contacts forlder to concat the files

"""organize the contacts report to a readable file"""
def makefile(all_files):
    columns = ['peptide_model','peptide_name','Atom1', 'Atom2','Clash','Covalent','VdW Clash','VdW','Proximal','Hydrogenbond',
                        'Weak_Hydrogen_Bond','Halogen_Bond','Ionic','Metal_Complex','Aromatic','Hydrophobic',
                        'Carbonyl','Polar','Weak_Polar','Interacting_entities']
    fulltable = pd.DataFrame(columns=columns)

    for filename in all_files: #filename = chk1_2ydimodel_2phkmuta_LLRLCTW_atms_top_0001.contacts (canoncila name)
        df = pd.read_table(filename,sep = '\t', names = ['Atom1', 'Atom2','Clash','Covalent','VdW Clash','VdW','Proximal','Hydrogenbond',
                        'Weak_Hydrogen_Bond','Halogen_Bond','Ionic','Metal_Complex','Aromatic','Hydrophobic',
                        'Carbonyl','Polar','Weak_Polar','Interacting_entities']) #put the header for each contacts_results
        peptide_name = filename[filename.find('muta')+5:filename.find('.')] #get the peptidename
        peptide_model = filename[:filename.find('_')] #get peptide model name, ex. chk1,chk2
        df['peptide_name'] = peptide_name #add column of peptide_name
        df['peptide_model'] = peptide_model #add column of peptide model name
        b = df.pop('peptide_name')
        df.insert(0,'peptide_name',b)
        a = df.pop('peptide_model')
        df.insert(0,'peptide_model',a)
        fulltable = pd.concat([fulltable,df])
    fulltable.to_csv('total_contacts',sep='\t',index = False)
    return fulltable

"""calculate the sum of each contacts"""
def make_calculation(fulltable):
    unipeptide = set(fulltable['peptide_name'])
     #get the dataframe
    columns = ['peptide_model','peptide_name','Clash','Covalent','VdW Clash','VdW','Proximal','Hydrogenbond',
                        'Weak_Hydrogen_Bond','Halogen_Bond','Ionic','Metal_Complex','Aromatic','Hydrophobic',
                        'Carbonyl','Polar','Weak_Polar']
    peptide_feature_contacts = pd.DataFrame(columns=columns)
    for peptide in unipeptide:
        peptide_feature_contacts.loc[peptide] =fulltable.loc[fulltable['peptide_name'] == peptide].iloc[:,4:19].sum()
        peptide_feature_contacts.loc[peptide,'peptide_name'] = peptide
        peptide_feature_contacts.loc[peptide,'peptide_model'] = fulltable.loc[fulltable['peptide_name'] == peptide]['peptide_model'][0]
        #df.loc['Total'] = df.iloc[:,4:19].sum()
        #peptide_feature_contacts = pd.concat([peptide_feature_contacts,df])
    
    peptide_feature_contacts[['Clash','Covalent','VdW Clash','VdW','Proximal','Hydrogenbond',
                        'Weak_Hydrogen_Bond','Halogen_Bond','Ionic','Metal_Complex','Aromatic','Hydrophobic',
                        'Carbonyl','Polar','Weak_Polar']] = peptide_feature_contacts[['Clash','Covalent','VdW Clash','VdW','Proximal','Hydrogenbond',
                        'Weak_Hydrogen_Bond','Halogen_Bond','Ionic','Metal_Complex','Aromatic','Hydrophobic',
                        'Carbonyl','Polar','Weak_Polar']].fillna(0.0).astype(int)
    peptide_feature_contacts = peptide_feature_contacts.sort_values(by=['peptide_model','peptide_name'])
    peptide_feature_contacts.to_csv('total_contacts_sum',sep='\t',index = False)
    return peptide_feature_contacts


if __name__ == "__main__":
    all_files = os.listdir('.') 
    makefile(all_files)

make_calculation(makefile(all_files))
