#! /usr/bin/python

import sys
import os

myFile = sys.argv[1]
myFO = "tmp.pdb"

myFrames = []
parsing = False

with open(myFile, "r") as f:
    for line in f:
        if "MODEL" in line:
            myFrames.append(line[6:9].strip())

my_identifier = "MODEL %s" % myFrames[-1]

with open(myFO, "w") as fo:

    with open(myFile, "r") as f:
        
        for line in f:
            
            if parsing:
                fo.write(line)
                parsing = True

            elif my_identifier == line[:9].strip():
                parsing = True
        
            else:
                parsing = False


cwd =  os.getcwd()
myDir =  cwd.split("/")[-1]
idx = myDir.find("_")

#cmd = "echo 0 | trjconv -f sampled.pos.pdb -s init.pdb -o tmp.pdb -dump %s"%(myFrames[-1])
#os.system(cmd)
#os.system("tail -n +4 tmp.pdb > tmp2.pdb")
os.system("sed -i.old '1s;^;CBLC >AB\\n;' tmp.pdb")
os.system("mv tmp.pdb ../my_structures_ctd/init%s.pdb" % myDir[idx:])
#os.system("rm tmp*")

