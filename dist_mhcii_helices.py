#! /usr/bin/env python

from MDAnalysis import *
#from MDAnalysis.analysis.align import *
import numpy
import math
import sys


u = Universe("init.pdb",sys.argv[1])
v = Universe("init.pdb")

# helices
a1 = u.selectAtoms("segid A and resid 46:77")
b1 = u.selectAtoms("segid B and resid 51:90")

f = open('helices_dist','w')

for ts in u.trajectory:
    distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
    f.write('%7.3f\n' % distance1)

f.close()

