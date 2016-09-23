#!/usr/bin/python

import sys
from MDAnalysis import *

myFile = sys.argv[1]
myDist = sys.argv[2]
myResidue = sys.argv[3]

u = Universe(myFile,myFile)

end = myFile.find(".pdb")
fileName = myFile[0:end]
myOutName = fileName + "_waterCts_%sA.dat"%myDist
myOut = open(myOutName, "w")

for ts in u.trajectory:
    #AME = u.selectAtoms("resname AME and name N6")
    nearWaters = u.selectAtoms("byres (resname TIP and around %s (resid %s and name C5))" % (myDist, myResidue))
    nearWaterCount = nearWaters.numberOfResidues()
    myOut.write("%s\n"%nearWaterCount)

myOut.close()
