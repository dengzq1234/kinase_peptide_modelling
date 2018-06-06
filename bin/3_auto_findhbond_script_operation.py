"""chimera output h-bond analysis"""


#alias chim='/Applications/Chimera.app/Contents/MacOS/chimera'
#alias pymol='/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL'
#ls -l > fileList.txt
#chim --nostatus --script "/Users/ziqideng/Desktop/replace_instance/test/auto_findhbond_test.py"

#!/usr/bin/env python
import sys
from chimera import runCommand
import os

"""change working directory to the appropriate pathway"""
os.chdir('/Users/ziqideng/Desktop/replace_instance/FindHbond/')
"""making a list of every muta_dock_model file"""
listname='/Users/ziqideng/Desktop/replace_instance/FindHbond/pdb_file/instanceList.txt'
#def findhbond():
#/Users/ziqideng/Desktop/replace_instance/chk1model_muta2phk/chk1_2ydimodel_2phkmuta_WVRANSR
os.chdir('/Users/ziqideng/Desktop/replace_instance/FindHbond/pdb_file/')

list = []
f = open(listname)
for line in f.readlines():
    list.append(line)
    f.close()


for i in list:
  cmd = 'open ' + i
  runCommand(cmd)
  runCommand('sel #all:.b')
  filename = 'hbonds lineWidth 3.0 intermodel false intraMol false intraRes true saveFile '+i[:-4]+'_hbondreport.txt'
  runCommand(filename) #command line
  """session pathway must consistent with the working directory"""
  sessionname = 'save '+'/Users/ziqideng/Desktop/replace_instance/FindHbond/pdb_file/'+i[:-4]+'_allmodles.py'
  runCommand(sessionname)
  runCommand('close all')
#runCommand('open /Users/ziqideng/Desktop/replace_instance/chk1model_muta2phk/chk1_2ydimodel_2phkmuta_WVRANSR/*.pdb')
#runCommand('sel #all:.b')
#runCommand('hbonds lineWidth 3.0 intermodel false intraMol false intraRes true saveFile /Users/ziqideng/Desktop/replace_instance/hbondtest3') #command line
os.chdir('/Users/ziqideng/Desktop/replace_instance/')
#from chimera import runCommand
#>>> runCommand('sel #all:.b')

