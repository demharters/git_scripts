#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]
my_struc = sys.argv[2]

u = Universe(my_struc,my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_glob_bg_dist2.dat'

a1 = u.selectAtoms("segid A and resid 29")
b1 = u.selectAtoms("segid B and resid 151")
#a1 = u.selectAtoms("segid A and resid 21")
#a2 = u.selectAtoms("segid A and resid 137")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a1.centerOfMass() - a2.centerOfMass())
	
    f.write('%7.3f\n' % distance1)

f.close()

