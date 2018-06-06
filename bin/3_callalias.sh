#!/usr/bin/sh
alias chim='/Applications/Chimera.app/Contents/MacOS/chimera'

cd pdb_file/
ls *.pdb > instanceList.txt 

chim --nostatus --script "auto_findhbond_script_template_operation.py"
