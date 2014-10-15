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
fout_name = my_traj[0:end] + '_DA10_DA5EPOR_dist.dat'

a1 = u.selectAtoms("atom A 1 CA") 
a2 = u.selectAtoms("atom B 1 CA")

#for DA10
#d1 = numpy.array([-45.61500132,  -4.9470005 , -56.26099968])
#d2 = numpy.array([ 48.125, 21.04400015, -50.88700151])

#for DA5
d1 = numpy.array([-44.71899962, -14.14800072, -56.98499966])
d2 = numpy.array([ 44.71900177, -14.14800072,  56.98400116])

f = open(fout_name,'w')

for ts in u.trajectory:

    A1 = a1.centerOfMass()
    A2 = a2.centerOfMass()

    C1 = A1 + d1
    C2 = A2 + d2
    
    distance1 = numpy.linalg.norm(C1 - C2)

    #    # write DA5, DA10 and DA330 distances to file
    f.write('%7.3f\n' % (distance1))

f.close()
