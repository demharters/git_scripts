#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]

u = Universe("init.pdb",my_traj)
v = Universe("init.pdb")

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_dist_kink.dat'

a1 = u.selectAtoms("segid A and resid 62:78")
b1 = u.selectAtoms("segid B and resid 51:79")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
    f.write('%7.3f\n' % distance1)

f.close()
