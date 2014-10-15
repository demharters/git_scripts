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
a1 = u.selectAtoms("atom L 207 CA") 
a2 = u.selectAtoms("atom O 205 CA")

f = open(fout_name,'w')

for ts in u.trajectory:

    A1 = a1.centerOfMass()
    
    c = numpy.array([ 2.91999817, -6.17900085, -1.69200039])
    A2_tmp = a2.centerOfMass()
    A2 = A2_tmp + c
    
    distance1 = numpy.linalg.norm(A1 - A2)
    distance2 = numpy.linalg.norm(A1 - A2_tmp)
    #    # write DA5, DA10 and DA330 distances to file
    f.write('%7.3f %7.3f\n' % (distance1,distance2))

f.close()
