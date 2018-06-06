#!/usr/bin/env python

# Carry out multiple arpeggio analysis and output them to one directory for top 10 models of each kinase-ligand complex
# input: directory name of the kinase-ligand complex top10 models
# ouput: the arpeggio contacts of 10models in contacts folded and other results in all_results/


import os, sys
import arpeggio
import shutil



bin_folder = '.'
input_data_folder = '../inputdata/'
result_folder = '../results/'
refined_path = result_folder + 'refined_ouput/'

os.chdir(bin_folder)

"""creat folders to collect results"""
all_results_path = result_folder +'all_results/'
contacts_path = result_folder +'contacts/'

os.mkdir(all_results_path)
os.mkdir(contacts_path)

"""move pdb"""
def copy_pdb(foldername):
    #for forlder in foldername:
    for i in os.listdir(foldername):
        if i.find('.pdb') > 0:
            shutil.copy(foldername + '/' + i, '.')
    return

"""filter only the pdb files"""
def find_pdb(folder_content):  
    folder_pdb = []
    for i in folder_content:
        if i.find('.pdb') > 0:
            folder_pdb.append(i)
    return folder_pdb

"""excute_arpeggio"""
def do_arpeggio():  
    for pdbfile in find_pdb(os.listdir('.')):
        os.system("python arpeggio.py %s -s /B//" % pdbfile)

"""filter the .contaces files"""
def find_contacts(folder_content):  
    file_contacts = []
    for i in folder_content:
        if i.find('.contacts') > 0:
            file_contacts.append(i)
    return file_contacts

"""collect the contacts files"""
def move_contacts():  # cope all the contacts file to contacts/

    for contact in find_contacts(os.listdir('.')):
        shutil.copy(contact, contacts_path)
    return


def move_otherfile():  # move others to all_results
    for file in os.listdir('.'):
        if file.startswith('chk') == True:
            if os.path.isfile(file) == True:
                shutil.move(file, all_results_path)
            else:
                pass


"""if __name__ == "__main__":
    foldername = sys.argv[1:]
    copy_pdb(foldername)"""

copy_pdb(refined_path)
do_arpeggio()
move_contacts()
move_otherfile()

