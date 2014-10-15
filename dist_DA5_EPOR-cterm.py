#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]
my_struc = sys.argv[2]
u = Universe(my_struc,my_traj)
v = Universe(my_struc)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_EPOR_dist.dat'

# VHVL modules
### I changed this to match res 205 in DA10
a1 = u.selectAtoms("atom C 187 CA") 
a2 = u.selectAtoms("atom D 187 CA")

f = open(fout_name,'w')

for ts in u.trajectory:

    A1 = a1.centerOfMass()
    A2 = a2.centerOfMass()
    
    distance1 = numpy.linalg.norm(A1 - A2)
    f.write('%7.3f\n' % (distance1))

f.close()
