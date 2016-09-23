#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]
my_struc = sys.argv[2]

u = Universe(my_struc,my_traj)

end = my_traj.find('.pdb')
fout_angle = my_traj[0:end] + '_angle.dat'

#a = u.selectAtoms("segid A and resid 78:182")
#b = u.selectAtoms("segid B and resid 91:190")
a = u.selectAtoms("segid A and resid 84:182")
b = u.selectAtoms("segid B and resid 95:190")

g = open(fout_angle,'w')

for ts in u.trajectory:

    a_1,a_2,a_3 = a.principalAxes()
    b_1,b_2,b_3 = b.principalAxes()

    angle1 = math.degrees(math.acos(numpy.dot(a_1,b_1)))
    angle2 = math.degrees(math.acos(numpy.dot(a_2,b_2)))
    angle3 = math.degrees(math.acos(numpy.dot(a_3,b_3)))

    if angle1 > 90:
	    angle1 = 180-angle1

    if angle2 > 90:
        angle2 = 180-angle2

    if angle3 > 90:
        angle3 = 180-angle3

       
    g.write('%7.3f %7.3f %7.3f\n' % (angle1,angle2,angle3))

g.close()

