#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_pdb = sys.argv[1]
my_traj = sys.argv[2]

u = Universe(my_pdb,my_traj)
v = Universe(my_pdb)

fout_name = 'pepdist123.dat' + str(my_traj)

# helices
a1 = u.selectAtoms("segid A and resid 28")
a2 = u.selectAtoms("segid A and resid 98")
a3 = u.selectAtoms("segid A and resid 116")

# peptide
c1 = u.selectAtoms("segid C and resid 1")
c2 = u.selectAtoms("segid C and resid 5")
c3 = u.selectAtoms("segid C and resid 9")


f = open(fout_name,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a1.centerOfMass() - c1.centerOfMass())
    distance2 = numpy.linalg.norm(a2.centerOfMass() - c2.centerOfMass())
    distance3 = numpy.linalg.norm(a3.centerOfMass() - c3.centerOfMass())


    f.write('%7.3f %7.3f %7.3f\n' % (distance1,distance2,distance3))

f.close()

