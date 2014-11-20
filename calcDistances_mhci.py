#! /usr/bin/python

from MDAnalysis import *
import os
import numpy
import math

rootdir = "/Users/sam/phd/mosaics/setup1_ampl700_shift0_period5k_uniform_0_15_mindist5"
my_distances = []

for subdirs, dirs, files in os.walk(rootdir, topdown=True):

    print subdirs
    print dirs
    splitName = subdirs.split("/")[-1]
    sequence = subdirs.split("_")[-3]
    IC50 = subdirs.split("_")[-2]
    repString = subdirs.split("_")[-1]
    
    try:
        repNumber = int(repString[3:])
    except:
        pass

    baseName = '_'.join(splitName.split("_")[:-1])
    trajName = baseName + "_sampled.pos.pdb"
    pdbName = baseName + ".pdb"
    myTrajPath = os.path.join(subdirs,trajName)
    myPDBPath = os.path.join(subdirs,pdbName)

    # rename trajectory to end in ".pdb" for MDAnalysis parser
    #newTrajName = baseName + "_sampled.pos.pdb"
    #myNewTrajPath = os.path.join(subdirs,newTrajName)
    
    #try:
    #    os.rename(myTrajPath,myNewTrajPath)
    #except:
    #    pass
   
    if "threePoint" in myTrajPath and repNumber <= 100:
        
        print subdirs
        u = Universe(myPDBPath,myTrajPath)

        a1 = u.selectAtoms("segid A and resid 26")
        a2 = u.selectAtoms("segid A and resid 97")
        a3 = u.selectAtoms("segid A and resid 117")
        c1 = u.selectAtoms("segid C and resid 1")
        c2 = u.selectAtoms("segid C and resid 5")
        c3 = u.selectAtoms("segid C and resid 9")

        for ts in u.trajectory:

            if ts.frame == 41:
                distance_Nter = numpy.linalg.norm(a1.centerOfMass() - c1.centerOfMass())
                distance_Middle = numpy.linalg.norm(a2.centerOfMass() - c2.centerOfMass())
                distance_Cter = numpy.linalg.norm(a3.centerOfMass() - c3.centerOfMass())

                avg_distance = (distance_Nter + distance_Middle + distance_Cter)/3
                
                my_distances.append(avg_distance)

