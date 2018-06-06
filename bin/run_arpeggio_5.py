#!/usr/bin/env python

# Carry out multiple arpeggio analysis and output them to one directory for top 10 models of each kinase-ligand complex
# input: directory name of the kinase-ligand complex top10 models
# ouput: the arpeggio contacts of 10models in contacts folded and other results in all_results/


import os, sys
import arpeggio
import shutil
import glob
print(os.getcwd())

bin_folder = os.getcwd()+'/'
input_data_folder = '../inputdata/'
result_folder = '../results/'
refined_path = result_folder + 'refined_output/'
#top10_models_path = refined_path+ "top10_models/"
#os.chdir(bin_folder)

"""creat folders to collect results"""
all_results_path = refined_path +'all_results/'
contacts_path = refined_path +'contacts/'

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
        print "we clean this file", pdbfile
        os.system("python " + bin_folder+"clean_pdb.py %s -rmw" % pdbfile)   #clean PDB file first
    for pdb_clean_file in find_pdb(os.listdir('.')):
        print "we analyse this clean_file", pdb_clean_file
        os.system("python " + bin_folder+"arpeggio.py %s -s /B//" % pdb_clean_file)

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
        shutil.copy(contact, 'contacts/')
    return


def move_otherfile():  # move others to all_results
    for file in os.listdir('.'):
        #if file.startswith('chk') == True:
        if os.path.isfile(file) == True:
            shutil.move(file, 'all_results/')
        else:
            pass


"""if __name__ == "__main__":
    foldername = sys.argv[1:]
    copy_pdb(foldername)"""
os.chdir(refined_path)
for model_foldername in os.listdir('.'):
    print model_foldername
    if model_foldername[-3:] == "top":
        top10_path = model_foldername +'/' + model_foldername+'10_models/'
        copy_pdb(top10_path)
        do_arpeggio()
        move_contacts()
        move_otherfile()


