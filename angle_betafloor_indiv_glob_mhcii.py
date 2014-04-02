#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]

u = Universe("init.pdb",my_traj)
v = Universe("init.pdb")

end = my_traj.find('.pdb')
fout_angle1 = my_traj[0:end] + '_betafloor_glob_angle1.dat'
fout_angle2 = my_traj[0:end] + '_betafloor_glob_angle2.dat'

c = u.selectAtoms("(segid A and resid 4:77) or (segid B and resid 1:90)")
d = u.selectAtoms("segid A and resid 78:182")
e = u.selectAtoms("segid B and resid 91:190")


g = open(fout_angle1,'w')
f = open(fout_angle2,'w')

for ts in u.trajectory:

    
    c_1,c_2,c_3 = c.principalAxes()
    d_1,d_2,d_3 = d.principalAxes()
    e_1,e_2,e_3 = e.principalAxes()
    
    angle1_1 = math.degrees(math.acos(numpy.dot(c_1,d_1)))
    angle1_2 = math.degrees(math.acos(numpy.dot(c_2,d_2)))
    angle1_3 = math.degrees(math.acos(numpy.dot(c_3,d_3)))

    angle2_1 = math.degrees(math.acos(numpy.dot(c_1,e_1)))
    angle2_2 = math.degrees(math.acos(numpy.dot(c_2,e_2)))
    angle2_3 = math.degrees(math.acos(numpy.dot(c_3,e_3)))

    if angle1_1 > 90:
        angle1_1 = 180-angle1_1

    if angle1_2 > 90:
        angle1_2 = 180-angle1_2

    if angle1_3 > 90:
        angle1_3 = 180-angle1_3

    if angle2_1 > 90:
        angle2_1 = 180-angle2_1

    if angle2_2 > 90:
        angle2_2 = 180-angle2_2

    if angle2_3 > 90:
        angle2_3 = 180-angle2_3



    g.write('%7.3f %7.3f %7.3f\n' % (angle1_1,angle1_2,angle1_3))
    f.write('%7.3f %7.3f %7.3f\n' % (angle2_1,angle2_2,angle2_3))



g.close()
f.close()
