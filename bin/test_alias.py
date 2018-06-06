#!/usr/bin/env python
import subprocess,os
subprocess.Popen('alias pymol=/Applications/MacPyMOL.app/Contents/MacOS/MacPyMOL',shell=True)

import pymol
filename = '/Users/ziqideng/Desktop/replace_instance_pipeline/memo_script_combination/ziqi_first_project/inputdata/2ydi.pdb'
pymol.cmd.load(filename)

pymol.finish_launching()

