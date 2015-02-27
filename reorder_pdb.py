#! /usr/bin/python

#usage: reorder_pdb.py template trajectory_frame
# takes naming and ordering from template structure and applies it to trajectory frame.
# This is necessary for analysis with X3DNA

import sys,os
from MDAnalysis import *

my_templateFile = sys.argv[1]
my_trajFile = sys.argv[2]

my_index = my_trajFile.find(".pdb")
my_outFileName = my_trajFile[:my_index] + "_reordered.pdb"

try:
    os.remove(my_outFileName)
except:
    pass

my_outFile = open(my_outFileName, "a")


def reordering(my_templateFile,my_trajFile):

    
    cter = 0
    currRes = 0
    completedRes = []
    completedRes2 = []
    AtomNumber = 0
    
    with open(my_templateFile, "r") as templ:
    
        for templ_line in templ:
       
            if templ_line[0:4] == "ATOM":
    
                if cter == 0:
                    completedRes.append(int(templ_line[24:27].strip()))
                else:
                    pass
    
                if completedRes[cter] in completedRes[:-1]:
                    pass 
                
                else:
                    increment_traj = True
                    currRes += 1
    
    
                completedRes.append(int(templ_line[24:27].strip()))
                cter += 1
                
                with open(my_trajFile, "r") as traj:
                    
                    for traj_line in traj:
                        
    
                        #if increment_traj:
                        #    currRes += 1
                        #    print(currRes)
                        #else:
                        #    pass
                        
                        if traj_line[13:17].strip() == templ_line[13:17].strip() and traj_line not in completedRes2:
     
                            AtomNumber += 1
                            newLine = traj_line[0:7] + "{:>4}".format(str(AtomNumber)) + templ_line[11:27] + traj_line[27:56] + templ_line[56:] 
                            my_outFile.write(newLine)
                            completedRes2.append(traj_line)
                            break
    
                        else:
                            pass

            else:
                pass


u = Universe(my_trajFile,my_trajFile)
v = Universe(my_templateFile)

myChains = u.selectAtoms("segid A or segid B")
currFrame = 1

for ts in u.trajectory:
    my_outFile.write("MODEL " + str(currFrame) + "\n")
    currFrame += 1
    u.atoms.write("model.pdb")
    reordering(my_templateFile,"model.pdb")
    my_outFile.write("TER\n")
    my_outFile.write("ENDMDL\n")
