#! /usr/bin/env python

from MDAnalysis import *
#from MDAnalysis.analysis.align import *
import numpy
import math
import sys

my_traj = sys.argv[1]
end = my_traj.find('.pdb')

u = Universe("init.pdb",my_traj)
v = Universe("init.pdb")

# linker 1
a1 = u.selectAtoms("segid J and resid 123")
b1 = u.selectAtoms("segid N and resid 1")

# linker 2
a2 = u.selectAtoms("segid K and resid 1")
b2 = u.selectAtoms("segid M and resid 123")

fout_dist = my_traj[0:end] + '_linker_dist.dat'

f = open(fout_dist,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
    distance2 = numpy.linalg.norm(a2.centerOfMass() - b2.centerOfMass())
    f.write('%7.3f %7.3f\n' % (distance1,distance2))

f.close()

