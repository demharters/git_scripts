#! /usr/bin/env python

from MDAnalysis import *
#from MDAnalysis.analysis.align import *
import numpy
import math
import sys

u = Universe("init.pdb",sys.argv[1])
v = Universe("init.pdb")

# residues
drA_glob = u.selectAtoms("segid A and resid 87:182")
dmC_glob = u.selectAtoms("segid C and resid 101:200")

drB_glob = u.selectAtoms("segid B and resid 91:190")
dmD_glob = u.selectAtoms("segid D and resid 91:193")

drA_bg = u.selectAtoms("segid A and resid 58:60")
dmC_bg = u.selectAtoms("segid C and resid 94:96")


f = open('dm_dr_dist.dat','w')

for ts in u.trajectory:
    
    distance1 = numpy.linalg.norm(drA_glob.centerOfMass() - dmC_glob.centerOfMass())
    distance2 = numpy.linalg.norm(drA_bg.centerOfMass() - dmC_bg.centerOfMass())
    distance3 = numpy.linalg.norm(drB_glob.centerOfMass() - dmD_glob.centerOfMass())
    f.write('%7.3f %7.3f %7.3f\n' % (distance1,distance2,distance3))

f.close()

