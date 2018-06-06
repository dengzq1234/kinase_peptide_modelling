#!/usr/bin/env python
import sys, time, os
import subprocess
import shlex
import fnmatch
import pymol
from pymol import cmd
import glob
import shutil
main_folder = os.getcwd()
bin_folder = './bin/'
input_data_folder = './inputdata/'
result_folder = './results/'

#load peptide instance
def AAcode_1_to_3(seq):
    """change peptide from 1 letter to 3 letter"""

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

def ATPchange(filename):
    """change atom names of ATP to the naming as ATP params """

    ATP_naming_dic = {"C1'": 'C5 ',
     'C2': 'C9',
     "C2'": 'C4 ',
     "C3'": 'C3 ',
     'C4': 'C10',
     "C4'": 'C2 ',
     'C5': 'C7',
     "C5'": 'C1 ',
     'C6': 'C8',
     'C8': 'C6',
     'N1': 'N4',
     'N3': 'N5',
     'N6': 'N3',
     'N7': 'N2',
     'N9': 'N1',
     'O1A': 'O7 ',
     'O1B': 'O4 ',
     'O1G': 'O1 ',
     "O2'": 'O13',
     'O2A': 'O8 ',
     'O2B': 'O5 ',
     'O2G': 'O2 ',
     "O3'": 'O12',
     'O3A': 'O9 ',
     'O3B': 'O6 ',
     'O3G': 'O3 ',
     "O4'": 'O11',
     "O5'": 'O10',
     'PA': 'P3',
     'PB': 'P2',
     'PG': 'P1'}    
    
    outfile = open(filename, "r+")
    with open(filename) as f1:   
        for line in f1:
            if line.startswith('HETATM') is True:
                line2 = ' '.join(line.split())
                atom_name = line2.split(' ')[2]
                if atom_name in ATP_naming_dic:
                    if atom_name == 'C4':
                        line = line.replace(atom_name,ATP_naming_dic[atom_name])
                        #print(line)
                        newline = line[:16]+ line[17:]
                        outfile.write(newline)
                        #print(newline)
                    else:
                        outfile.write(line.replace(atom_name,ATP_naming_dic[atom_name]))
                else:
                    outfile.write(line)
                #print(atom_name)
            else:
                outfile.write(line)
    outfile.close()
    return 

def Mutagenesis(kinase1, model, template,peptide_instance):
    """superposition the model and template, remove template and leave peptide behind """

    """replace peptide with instance peptide"""
    list_name = peptide_instance
    with open(input_data_folder+list_name,'r') as f:
        instances = f.readlines()
    instance = [x.strip() for x in instances]
    print "your peptide: ", instance 
    
    for pep in instance:
        cmd.delete('all')
        cmd.fetch(model)  # model_candidate, ex.chk1,chk2...
        cmd.remove("hetatm") #remove the nonstandard residues
        cmd.fetch(template)  #mutagenesis template, ex.2phk
        peptide_template = cmd.get_fastastr( "/"+template+'//B') #get the peptide from the template and generate another one for mutagenesis 
        peptide_template = peptide_template + 'G' #peptide of 2phk is 7 amino acid long, when our input peptide is 8 aa, we need to plus one character
        
        
        for aa in peptide_template[6:].lower(): #creat template_peptide for mutagenesis

            cmd._alt(aa)
        
        
        firstaa = AAcode_1_to_3(peptide_template[6]) #translate template_peptide to 3 letter
        low_firstaa = firstaa[0].lower()

        cmd.alter(low_firstaa, 'chain = "B"') #select this template_peptide        
        cmd.show_as("cartoon")
        cmd.align(model, template) #superpostion of model and template
        cmd.align(low_firstaa,template) #superpostion of template_peptide and template

        remove_part = "("+template+" and not resn ATP"+")" 
        cmd.select("remove_part",remove_part) 

        
        cmd.remove("remove_part") #remove the template except for ATP, there are only model and template_peptide
        cmd.remove("resn hoh") #remove water
        cmd.wizard("mutagenesis")
        peptide_position = 0

        for i in pep:
            #the peptide_position starting point depends the first position of mutagenesis peptide
            #pymol's peptide start from 1, not 0
            mutagenesis_template = '/'+low_firstaa+ '///' + str(peptide_position + 2) # because of 2phk start at 2nd of peptide
            #mutagenesis_template = '/' + template + '//B/' + str(peptide_position + 2)
            cmd.get_wizard().do_select(mutagenesis_template)  # select peptide position of mutation
            replace_aminoacid = AAcode_1_to_3(pep)[peptide_position]
            cmd.get_wizard().set_mode(replace_aminoacid)  # select which residue want to mutate to
            cmd.get_wizard().apply()
            peptide_position += 1
        filename = kinase1 + '_' + model + 'model_' + template + 'muta_' + pep + '.pdb' #build the canonical name
        cmd.save(filename)
        ATPchange(filename) #change ATP naming to the format of ATP.params
    cmd.wizard(None)
    return

#execute Mutagenesis step in pymol nogui mode
cmd.extend("Mutagenesis_model", Mutagenesis)
if __name__ == "__main__":
    kinase1 = sys.argv[1]
    model = sys.argv[2]
    template = sys.argv[3]
    peptide_instance = sys.argv[4]
    command = 'Mutagenesis '+ sys.argv[1]+','+sys.argv[2]+','+sys.argv[3]+','+sys.argv[4]


#collect the models for the prepack mode
os.mkdir('results/')
subprocess.call(['pymol', sys.argv[0],'-qcd',command])
#print(glob.glob('*.pdb'))
pdb_filenames = glob.glob('*.pdb')
print pdb_filenames
core_list =[]
for pdb in pdb_filenames:
    shutil.move(pdb, result_folder + pdb)
    #get_atoms(pdb)
    core = pdb[:pdb.rfind('.')]
    core_list.append(core)
print "this is core_list:", core_list


"""get the original filename. 
input =  "../results/chk1_2ydimodel_2phkmuta_LLRLCTW.pdb" 
output(core) = chk1_2ydimodel_2phkmuta_LLRLCTW
"""




"""after the the atms file, excute the prepack mode"""

"""get the prepacked results and collect them to prepacked folder"""
prepacked_path = result_folder + 'prepacked_output/'
os.mkdir(prepacked_path)


def prepack(pdb):
	"""prepare the prepack flags"""
    with open('prepack_flags', 'w') as flag:  # prepare the prepack flags
        flag.write('-database /g/easybuild/x86_64/CentOS/7/nehalem/software/Rosetta/2017.45.59812-foss-2016b/database\n')
        flag.write('-s' + ' ' + result_folder + pdb + '\n')
        flag.write('-extra_res_fa inputdata/param_file/ATP.params\n')
        #flag.write('-constraints:cst_file /g/gibson/deng/rosetta_analysis/script_improve_test/peptite_constrain/inputdata/param_file/constraint_Psite\n')
        flag.write('-ex1\n')
        flag.write('-ex2aro\n')
        flag.write('-use_input_sc\n')
        flag.write('-flexpep_score_only\n')
        flag.write('-ignore_unrecognized_res\n')
        flag.write('-out:path:all '+prepacked_path+'\n')
        flag.write('-flexpep_prepack\n')
    subprocess.call(['FlexPepDocking.mpi.linuxgccrelease', '@prepack_flags'])


def get_allpdbfile(path):
    """get all pdb file"""
    f_list = os.listdir(path)
    prepack_pdb = []
    print f_list
    for i in f_list:
        # os.path.splitext():get the tail
        if os.path.splitext(i)[1] == '.pdb':
            prepack_pdb.append(i)
    print prepack_pdb 
    return prepack_pdb 

#prepack them
for atom_pdb in get_allpdbfile(result_folder):
    print(atom_pdb)
    prepack(atom_pdb)

#rename and output to prepack file    
for core in core_list:    
    pdb_atms_prepcked_filename = core + '_0001.pdb'
    print "we are gonna change this name:", pdb_atms_prepcked_filename
    shutil.move(prepacked_path + pdb_atms_prepcked_filename, prepacked_path + core + '_top.pdb')

#creat the directory for refine models
refined_path = result_folder + 'refined_output/'
os.mkdir(refined_path)

def low_refine(pdb_atms_prepcked_filename):
    """prepare the lowrefinement flags document"""
    pdb_atms_prepcked_filename_core = os.path.splitext(pdb_atms_prepcked_filename)[0] #core no pdb tail
    refined_model_path = refined_path + pdb_atms_prepcked_filename_core + '/'
    os.mkdir(refined_model_path)  
    with open('refine_flags_low', 'w') as flag:  # prepare the refine flags
        flag.write('-database /g/easybuild/x86_64/CentOS/7/nehalem/software/Rosetta/2017.45.59812-foss-2016b/database\n')
        flag.write('-s' + ' ' + prepacked_path + pdb_atms_prepcked_filename + '\n')
        flag.write('-extra_res_fa inputdata/param_file/ATP.params\n')
        flag.write('-constraints:cst_fa_file inputdata/param_file/constraint_super\n')
        flag.write('-constraints:cst_fa_weight 5\n')
        flag.write('-ex1\n')
        flag.write('-ex2aro\n')
        flag.write('-use_input_sc\n')
        #flag.write('-unboundrot'+ ' ' + prepacked_path + pdb_atms_prepcked_filename + '\n')
        flag.write('-lowres_preoptimize\n') #low resolution before high resolution
        flag.write('-nstruct 100\n') #number of model generated
        flag.write('-ignore_unrecognized_res\n')
        flag.write('-scorefile'+' ' + 'low_refine.score.sc\n') #score file
        flag.write('-out:path:all '+refined_model_path+'\n')
        flag.write('-pep_refine\n')
        flag.write('-flexpep_score_only')
    #subprocess.call(['/g/sbgrid/programs/x86_64-linux/rosetta/2017.36.59679/main/source/bin/FlexPepDocking.mpi.linuxgccrelease', '@refine_flags_low'])
    subprocess.call(['mpirun', '-n',"30",'FlexPepDocking.mpi.linuxgccrelease', '@refine_flags_low'])
    for low_refine_model in glob.glob(refined_model_path+'*.pdb'):
        print "this is low_refine_model name", low_refine_model
        new_name_low = os.path.splitext(low_refine_model)[0] +"_low.pdb"
        shutil.move(low_refine_model, new_name_low)
    return


"""get the refined results and collect them to refined folder"""
for prepacked_pdb in os.listdir(prepacked_path): #happend in prepack folder
    if os.path.splitext(prepacked_pdb)[1] == ".pdb":
        print(prepacked_pdb)
        low_refine(prepacked_pdb)

#refined_path_high = result_folder + 'refined_output_high/'
#os.mkdir(refined_path_high)

def high_refine(pdb_atms_prepcked_filename):
	"""prepare the highrefinement flags document"""
    pdb_atms_prepcked_filename_core = os.path.splitext(pdb_atms_prepcked_filename)[0] #with _top
    refined_model_path = refined_path + pdb_atms_prepcked_filename_core + '/'
    #os.mkdir(refined_model_path) 
    with open('refine_flags_high', 'w') as flag:  # prepare the refine flags
        flag.write('-database /g/easybuild/x86_64/CentOS/7/nehalem/software/Rosetta/2017.45.59812-foss-2016b/database\n')
        flag.write('-s' + ' ' + prepacked_path + pdb_atms_prepcked_filename + '\n')
        flag.write('-extra_res_fa inputdata/param_file/ATP.params\n')
        flag.write('-constraints:cst_fa_file inputdata/param_file/constraint_super\n')
        flag.write('-constraints:cst_fa_weight 5\n')
        #flag.write("-beta_nov15\n")
        flag.write('-ex1\n')
        flag.write('-ex2aro\n')
        flag.write('-use_input_sc\n')
        flag.write('-pep_refine\n')
        #flag.write('-unboundrot'+ ' ' + prepacked_path + pdb_atms_prepcked_filename + '\n')
        flag.write('-nstruct 100\n') #number of model generated
        flag.write('-ignore_unrecognized_res\n')
        flag.write('-out:path:all '+refined_model_path+'\n')
        flag.write('-scorefile'+' ' + 'high_refine.score.sc\n') #score file
        flag.write('-flexpep_score_only')
    subprocess.call(['mpirun', '-n','30','FlexPepDocking.mpi.linuxgccrelease', '@refine_flags_high'])
    #subprocess.call(['/g/sbgrid/programs/x86_64-linux/rosetta/2017.36.59679/main/source/bin/FlexPepDocking.mpi.linuxgccrelease', '@refine_flags_high'])
    for high_refine_model in glob.glob(refined_model_path+'*[0-9].pdb'):
        print "this is high_refine_model name", high_refine_model
        new_name_high = os.path.splitext(high_refine_model)[0] +"_high.pdb"
        shutil.move(high_refine_model, new_name_high)
    return 

#high refine models
for prepacked_pdb in os.listdir(prepacked_path):
    if os.path.splitext(prepacked_pdb)[1] == ".pdb":
        print(prepacked_pdb)
        high_refine(prepacked_pdb)

#get all the model name
os.chdir(refined_path)
refined_model_list = []
for refined_folder in os.listdir('.'):
    print refined_folder
    refined_model_list.append(refined_folder)
os.chdir(main_folder)

#move the ouput files to target folder
shutil.move('refine_flags_high',result_folder+'refine_flags_high')
shutil.move('refine_flags_low',result_folder+'refine_flags_low')
shutil.move('prepack_flags',result_folder+'prepack_flags')


def select_top10(refined_foldername):
	"""select top10 among all models based on score"""
    refined_model_path = refined_foldername+'/'
    os.chdir(refined_path+refined_foldername)
    print(os.getcwd())
    result=[]
    with open('high_refine.score.sc','r') as f1:
        for line in f1:
            line = ' '.join(line.split())
            name = line.split(' ')[-1] 
            if name.startswith('description') or name.startswith('SEQUENCE')  is True:
                pass
            else:
                name = name+'_high'
            line = line.split(' ')
            if len(line)>1 and line[1].startswith('total_score')is False:
                score = line[1]
                I_sc = line[5]
                f_atr = line[17]
                f_dun = line[18]
                fa_rep = line[22]
                fa_sol = line[23]
                hbond_sc = line[33]
                pep_sc = line[38]
                pep_sc_noref = line[39]
                rmsALL_if = line[47]
                rmsBB = line[48]
                rmsBB_if = line[49]
                startRMSbb = line[-6]
                #print(score,name)
                result.append((score,I_sc,f_atr,f_dun,fa_rep,fa_sol,hbond_sc,pep_sc,pep_sc_noref,rmsALL_if,rmsBB,rmsBB_if,startRMSbb,name))
    #print(result)
    with open('low_refine.score.sc','r') as f2:
        for line in f2:
            line = ' '.join(line.split())
            name = line.split(' ')[-1] 
            if name.startswith('description') or name.startswith('SEQUENCE')  is True:
                pass
            else:
                name = name+'_low'
            line = line.split(' ')
            if len(line)>1 and line[1].startswith('total_score')is False:
                score = line[1]
                I_sc = line[5]
                f_atr = line[17]
                f_dun = line[18]
                fa_rep = line[22]
                fa_sol = line[23]
                hbond_sc = line[33]
                pep_sc = line[38]
                pep_sc_noref = line[39]
                rmsALL_if = line[47]
                rmsBB = line[48]
                rmsBB_if = line[49]
                startRMSbb = line[-6]
                #print(score,name)
                result.append((score,I_sc,f_atr,f_dun,fa_rep,fa_sol,hbond_sc,pep_sc,pep_sc_noref,rmsALL_if,rmsBB,rmsBB_if,startRMSbb,name))
    #print(result)
    def make_string(X):
        """input: a list of float, last element is name.
           output: a list of string """
        out_list=[]
        for element in X:
            out_list.append(str(element))
        return tuple(out_list)

    result_clean= []
    for i in result:
        if i[0][0][0]=='-':
            name = i[-1]
            scores = list(map(lambda x:float(x), i[:-1]))
            combine = tuple(scores+[name])
            result_clean.append(combine)
    result_clean = sorted(result_clean)
    top10 = result_clean[0:10]
    print(top10)
    total_refine_score = "top10_total_refine.score.sc"
    os.mkdir(refined_foldername+"10_models/")

    with open(refined_foldername+"10",'w') as f3:
        for i in top10:
            f3.write(i[13]+'\n')

    with open(total_refine_score,'w') as f4:
        head= ["score","I_sc","f_atr","f_dun","fa_rep","fa_sol","hbond_sc","pep_sc","pep_sc_noref","rmsALL_if","rmsBB","rmsBB_if","startRMSbb","discription"]
        f4.write("\t".join(head)+'\n')
        for i in top10:
            f4.write("\t".join(make_string(i))+'\n')
            
    lines = [line.rstrip('\n') for line in open(refined_foldername+"10",'r')]
    print('this is line:', line)
    for model in glob.glob('*.pdb'):
        model = os.path.splitext(model)[0]
        print (model)
        if model in lines:
            print("enter top10 selection ", model)
            shutil.copyfile(model+'.pdb', refined_foldername+"10_models/"+model+'.pdb')
    
    shutil.copyfile(total_refine_score,refined_foldername+"10_models/"+total_refine_score)
    os.chdir(main_folder)
    return

for refined_foldername in os.listdir(refined_path):
    print refined_foldername
    select_top10(refined_foldername)

"""step arpeggio"""
#bin_folder = './bin/'
#input_data_folder = './inputdata/'
#result_folder = './results/'

os.chdir(bin_folder)
os.system("python run_arpeggio_5.py")

"""step integration"""
os.system("python integration_ziqi_6.py")
os.chdir(main_folder)


