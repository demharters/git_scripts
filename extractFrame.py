#! /usr/bin/python

import sys
import os

myFile = sys.argv[1]
frame = int(sys.argv[2])

myFO = "tmp.pdb"

myFrames = []
parsing = False

with open(myFile, "r") as f:
    for line in f:
        if "MODEL" in line:
            myFrames.append(line.strip())

my_identifier = myFrames[frame-1]
print my_identifier

with open(myFO, "w") as fo:

    with open(myFile, "r") as f:
        
        for line in f:
            
            if "MODEL" in line:
                parsing = False

            if parsing:
                if "TITLE" in line or "REMARK" in line:
                    pass
                else:
                    fo.write(line)
                
                parsing = True

            elif my_identifier == line.strip():
                parsing = True
        
            else:
                parsing = False

myNewName = myFile.split(".")[0] + "_frame_" + str(frame) + ".pdb"
os.system("mv tmp.pdb %s" % myNewName)
