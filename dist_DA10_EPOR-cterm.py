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

# residues
a1 = u.selectAtoms("segid L and (resid 205)")
        
        #31,32,52,53,54,100,101,102,103,106,107,109,110,159,161,179,180,181,182")

b1 = u.selectAtoms("segid O and resid 205")
        
        #31,32,52,53,54,100,101,102,103,106,107,109,110,159,161,179,180,181,182")

fout_dist = my_traj[0:end] + '_EPOR_dist.dat'

f = open(fout_dist,'w')
#g = open('angle','w')

for ts in u.trajectory:
        
        distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
        #distance2 = numpy.linalg.norm(a2.centerOfMass() - b2.centerOfMass())
        #distance3 = numpy.linalg.norm(a3.centerOfMass() - b3.centerOfMass())
	#distance4 = numpy.linalg.norm(a4.centerOfMass() - b4.centerOfMass())
	
	#a4_1,a4_2,a4_3 = a4.principalAxes()
	#b4_1,b4_2,b4_3 = b4.principalAxes()
       # helix12_1,helix12_2,helix12_3 = helix12.principalAxes()
       # helix21_1,helix21_2,helix21_3 = helix21.principalAxes()
       # helix22_1,helix22_2,helix22_3 = helix22.principalAxes()

	#angle = math.degrees(math.acos(numpy.dot(a4_1,b4_1)))
       # angle2 = math.degrees(math.acos(numpy.dot(helix21_1,helix22_1)))

	#if angle > 90:
	#	angle = 180-angle
                
       # print "%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f" % (ts.frame,rmsd0,rmsd1,rmsd2,distance1,distance2,angle1,angle2)

	#f.write('%7.3f %7.3f % 7.3f % 7.3f\n' % (distance1,distance2,distance3,distance4))
	f.write('%7.3f\n' % (distance1))
	#g.write('%7.3f\n' % angle)

f.close()
#g.close()

