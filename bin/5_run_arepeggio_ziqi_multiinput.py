#!/usr/bin/env python

#Carry out multiple arpeggio analysis and output them to one directory for top 10 models of each kinase-ligand complex
#input: directory name of the kinase-ligand complex top10 models
#ouput: the arpeggio contacts of 10models in contacts folded and other results in all_results/


import os,sys
import arpeggio
import shutil
#copy pdb files to the excute directory(single directory)

os.mkdir('contacts')
os.mkdir('all_results')

def copy_pdb(foldername):
    for forlder in foldername:
        for i in os.listdir(forlder):
            if i.find('.pdb') >0:
                shutil.copy(forlder + '/' +i,'.')
    return 


def find_pdb(folder_content): #filter only pdb files
    folder_pdb = []
    for i in folder_content:
        if i.find('.pdb') > 0:
            folder_pdb.append(i)
    return folder_pdb 


def do_arpeggio(): #excute arpeggio
    for pdbfile in find_pdb(os.listdir('.')):
        os.system("python arpeggio.py %s -s /B//" % pdbfile)


def find_contacts(folder_content): #filter the .contacts file
    file_contacts = []
    for i in folder_content:
        if i.find('.contacts') > 0:
            file_contacts .append(i)
    return file_contacts

def move_contacts(): #cope all the contacts file to contacts/
	
    for contact in find_contacts(os.listdir('.')):
        shutil.copy(contact,'contacts/')
    return


def move_otherfile(): #move others to all_results
    for file in os.listdir('.'):
        if file.startswith('chk') == True:
            if os.path.isfile(file) == True:
                shutil.move(file,'all_results')
            else:
                pass

if __name__ == "__main__":
    foldername = sys.argv[1:]
    copy_pdb(foldername)

do_arpeggio()
move_contacts()
move_otherfile()