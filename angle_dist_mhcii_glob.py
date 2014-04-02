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

# globular domains 
a1 = u.selectAtoms("segid A and resid 78:182")
b1 = u.selectAtoms("segid B and resid 91:190")

fout_dist = my_traj[0:end] + '_glob_dist.dat'
fout_angle = my_traj[0:end] + '_glob_angle.dat'

f = open(fout_dist,'w')
g = open(fout_angle,'w')

#iterate through trajectory and calculate distances and angles between the two globular domains of chains A and B
for ts in u.trajectory:
        
    # get distance between globular domains
    distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())

    # get principal Axes
    a1_1,a1_2,a1_3 = a1.principalAxes()
    b1_1,b1_2,b1_3 = b1.principalAxes()
	
    # calcualte angle between the two first principle axes
    angle = math.degrees(math.acos(numpy.dot(a1_1,b1_1)))
    
    if angle > 90:
        angle = 180-angle
                
    f.write('%7.3f\n' % distance1)
    g.write('%7.3f\n' % angle)

f.close()
g.close()

