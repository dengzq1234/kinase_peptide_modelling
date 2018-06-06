#!/usr/bin/env python

"""
input1:kinase srtucture as model(pdb file)

input2:peptide template which will be replaced by instance(so far is 2phk)

input3:peptide which try to dock to the model(peptide_candidate)

output:a pdb file with model and peptide instance(predock_model)

"""

import sys, time, os
import subprocess
import pymol
from pymol import cmd

import __main__
#__main__.pymol_argv = ['pymol', '-c']

"""set all the necessary pathway"""
bin_folder = '.'
input_data_folder = '../inputdata/'
result_folder = '../results/'
os.chdir(bin_folder)

"""load the peptide_instance """
list_name = 'peptide_candidate'
with open(input_data_folder+list_name,'r') as f:
    instance = f.readlines()
instance = [x.strip() for x in instance]

"""convert 1 letter amino acid to 3 letter abr"""
def AAcode_1_to_3(seq):

    d = {'C': 'CYS',
         'D': 'ASP',
         'S': 'SER',
         'Q': 'GLN',
         'K': 'LYS',
         'I': 'ILE',
         'P': 'PRO',
         'T': 'THR',
         'F': 'PHE',
         'N': 'ASN',
         'G': 'GLY',
         'H': 'HIS',
         'L': 'LEU',
         'R': 'ARG',
         'W': 'TRP',
         'A': 'ALA',
         'V': 'VAL',
         'E': 'GLU',
         'Y': 'TYR',
         'M': 'MET'}
    upper_seq = seq.upper()
    triple_seq = []
    for b in upper_seq:
        triple_seq.append(d[b])
    return triple_seq

"""creat Mutagenesis extend in pymol"""
def Mutagenesis(kinase1, model, template):


    """superposition the model and template, remove template and left peptide behind """

    """replace peptide with instance peptide"""
    for pep in instance:

        cmd.delete('all')
        cmd.fetch(model)  # model_candidate(chk1,chk2model....)
        cmd.fetch(template)  # mutagenesis template
        cmd.show_as("cartoon")
        cmd.align(model, template) #superpostion of model and template
        remove_chain = "/" + template + "//A"
        cmd.remove(remove_chain)
        cmd.wizard("mutagenesis")
        peptide_position = 0

        for i in pep:
            # the peptide_position starting point depends the first position of mutagenesis peptide
            #pymol's peptide start from 1, not 0
            mutagenesis_template = '/' + template + '//B/' + str(peptide_position + 2) # because of 2phk start at 2nd of peptide
            cmd.get_wizard().do_select(mutagenesis_template)  # select peptide position of mutation
            replace_aminoacid = AAcode_1_to_3(pep)[peptide_position]
            cmd.get_wizard().set_mode(replace_aminoacid)  # select which residue want to mutate to
            cmd.get_wizard().apply()
            peptide_position += 1
        filename = result_folder + kinase1 + '_' + model + 'model_' + template + 'muta_' + pep + '.pdb' #build the canonical name
        cmd.save(filename)

    cmd.wizard(None)
    return

cmd.extend("Mutagenesis_model", Mutagenesis)

if __name__ == "__main__":
    kinase1 = sys.argv[1]
    model = sys.argv[2]
    template = sys.argv[3]

    command = 'Mutagenesis '+ sys.argv[1]+','+sys.argv[2]+','+sys.argv[3]

#pymol.finish_launching(['pymol', '-d', command])

#pymol.finish_launching(['pymol','-qc'])
#subprocess.call(['pymol', '-d', command])
#save myTraj.pdb, myMDTrajectory, state=0

#usage: $pymol 1_auto_mutagenesis_scripttemplate.py -d "Mutagenesis chk1,2ydi,2phk"
#pymol.cmd.quit()

"""call the function at command line"""
subprocess.call(['pymol', '1_auto_mutagenesis_scripttemplate.py','-qcd',command])


