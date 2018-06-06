#!/usr/bin/env python
#excute flexpepdock

#input: predock_complex
#output: 10_docked_complex

import os, sys
import subprocess
import shutil
import fnmatch

"""load the pathway"""
bin_folder = '.'
input_data_folder = '../inputdata/'
result_folder = '../results/'
os.chdir(bin_folder)


#subprocess.call(['module', 'load','Rosetta'])

"""extract all the atoms """
def get_atoms(pdb_filename):
    with open(pdb_filename, 'r') as f:
        pdb = f.readlines()
    pdb = [x.strip() for x in pdb]
    pdb_atms = []
    for line in pdb:
        if line.startswith('ATOM'):
            pdb_atms.append(line)
    pdb_atms_filename = pdb_filename[11:-4] + '_atms.pdb'
    # pdb_atms_filename = pdb_filename[:pdb_filename.find('.')]+'_atms.pdb'
    wf = open(result_folder + pdb_atms_filename, 'w')
    for line in pdb_atms:
        wf.write(line + '\n')
    wf.close()
    os.remove(pdb_filename)
    return pdb_atms


if __name__ == "__main__":
    pdb_filename = sys.argv[1]
    get_atoms(pdb_filename)

"""get the original filename. 
input =  "../results/chk1_2ydimodel_2phkmuta_LLRLCTW.pdb" 
output(core) = chk1_2ydimodel_2phkmuta_LLRLCTW
"""
core = pdb_filename[pdb_filename.rfind('/')+1:pdb_filename.rfind('.')]


"""after the the atms file, excute the prepack mode"""

"""prepare the prepack falgs"""
def prepack(pdb_atms_filename):
    with open('prepack_flags', 'w') as flag:  # prepare the prepack flags
        flag.write('-database /g/easybuild/x86_64/CentOS/7/nehalem/software/Rosetta/3.8-foss-2016b/database\n')
        flag.write('-s' + ' ' + result_folder + pdb_atms_filename + '\n')
        flag.write('-ex1\n')
        flag.write('-ex2aro\n')
        flag.write('-flexpep_prepack\n')
    subprocess.call(['FlexPepDocking.mpi.linuxgccrelease', '@prepack_flags'])


#if __name__ == "__main__":
#    pdb_atms_filename = sys.argv[1]
#    prepack(pdb_atms_filename) """

"""excute prepack to the pre_dock_complex"""
for atom_pdb in os.listdir(result_folder):
    print(atom_pdb)
    prepack(atom_pdb)

"""get the prepacked results and collect them to prepacked folder"""
prepacked_path = result_folder + 'prepacked_ouput/'
os.mkdir(prepacked_path)

# pdb_atms_prepcked_filename = pdb_atms_filename[:pdb_atms_filename.find('.')]+'_0001.pdb'
# shutil.move(pdb_atms_prepcked_filename , path+pdb_atms_prepcked_filename)
pdb_atms_prepcked_filename = core + '_atms_0001.pdb'
# os.rename(result_folder+pdb_atms_prepcked_filename,result_folder+pdb_atms_filename[:pdb_atms_filename.find('.')]+'_top.pdb' )
os.rename(pdb_atms_prepcked_filename, prepacked_path + core + '_atms_top.pdb')


# os.remove(pdb_atms_filename)


"""prepare refine step"""
def refine(pdb_atms_prepcked_filename):
    with open('refine_flags', 'w') as flag:  # prepare the refine flags
        flag.write('-database /g/easybuild/x86_64/CentOS/7/nehalem/software/Rosetta/3.8-foss-2016b/database\n')
        flag.write('-s' + ' ' + prepacked_path + pdb_atms_prepcked_filename + '\n')
        flag.write('-ex1\n')
        flag.write('-ex2aro\n')
        flag.write('-pep_refine\n')
        flag.write('-lowres_preoptimize\n') #low resolution before high resolution
        flag.write('-nstruct 2\n') #number of model generated
        flag.write('-scorefile refine.score.sc\n') #score file
        flag.write('-flexpep_score_only')
    subprocess.call(['FlexPepDocking.mpi.linuxgccrelease', '@refine_flags'])
    return


"""if __name__ == "__main__":
    pdb_atms_prepcked_filename = sys.argv[1]
    refine(pdb_atms_prepcked_filename) """

"""get the refined results and collect them to refined folder"""
for prepacked_pdb in os.listdir(prepacked_path):
    print(prepacked_pdb)
    refine(prepacked_pdb)

refined_path = result_folder + 'refined_ouput/'
os.mkdir(refined_path)

# refined_pdb_name = pdb_atms_prepcked_filename[:pdb_atms_prepcked_filename.find('.')]+'_0001.pdb'
# shutil.move(prepacked_pdb_name , path+prepacked_pdb_name)
"""move files"""
for file in os.listdir('.'):
    if fnmatch.fnmatch(file, '*_000*'):
        shutil.move(file, refined_path + file)
shutil.move('refine.score.sc', refined_path + 'refine.score.sc')
shutil.move('score.sc', prepacked_path+'score.sc')
shutil.move('refine_flags',result_folder+'refine_flags')
shutil.move('prepack_flags',result_folder+'prepack_flags')
# os.remove(pdb_atms_prepcked_filename) """
