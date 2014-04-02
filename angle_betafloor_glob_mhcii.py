#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]

u = Universe("init.pdb",my_traj)
v = Universe("init.pdb")

end = my_traj.find('.pdb')
fout_angle = my_traj[0:end] + '_betafloor_glob_angle.dat'

c = u.selectAtoms("(segid A and resid 4:77) or (segid B and resid 1:90)")
d = u.selectAtoms("(segid A and resid 78:182) or (segid B and resid 91:190)")

g = open(fout_angle,'w')


for ts in u.trajectory:

    
    c_1,c_2,c_3 = c.principalAxes()
    d_1,d_2,d_3 = d.principalAxes()
    
    angle1 = math.degrees(math.acos(numpy.dot(c_1,d_1)))
    angle2 = math.degrees(math.acos(numpy.dot(c_2,d_2)))
    angle3 = math.degrees(math.acos(numpy.dot(c_3,d_3)))

    if angle1 > 90:
        angle1 = 180-angle1

    if angle2 > 90:
        angle2 = 180-angle2

    if angle3 > 90:
        angle3 = 180-angle3

    g.write('%7.3f %7.3f %7.3f\n' % (angle1,angle2,angle3))

g.close()

