#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]
my_struc = sys.argv[2]

#u = Universe("init.pdb",my_traj)
#v = Universe("init.pdb")
u = Universe(my_struc,my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_dist.dat'

# residues
#a1 = u.selectAtoms("segid A and resid 46")
#b1 = u.selectAtoms("segid B and resid 90")
#a2 = u.selectAtoms("segid A and resid 67")
#b2 = u.selectAtoms("segid B and resid 64")
#a3 = u.selectAtoms("segid A and resid 77")
#b3 = u.selectAtoms("segid B and resid 51")

a1 = u.selectAtoms("segid A and resid 50:51")
b1 = u.selectAtoms("segid B and resid 85:86")
a2 = u.selectAtoms("segid A and resid 53:55")
b2 = u.selectAtoms("segid B and resid 78:83")
a3 = u.selectAtoms("segid A and resid 60:65")
b3 = u.selectAtoms("segid B and resid 65:70")
a4 = u.selectAtoms("segid A and resid 68:73")
b4 = u.selectAtoms("segid B and resid 56:61")

# helices
#a5 = u.selectAtoms("segid A and resid 50:73")
#b5 = u.selectAtoms("segid B and resid 56:86")


f = open(fout_name,'w')
#g = open('angle','w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
    distance2 = numpy.linalg.norm(a2.centerOfMass() - b2.centerOfMass())
    distance3 = numpy.linalg.norm(a3.centerOfMass() - b3.centerOfMass())
    distance4 = numpy.linalg.norm(a4.centerOfMass() - b4.centerOfMass())
    #distance5 = numpy.linalg.norm(a5.centerOfMass() - b5.centerOfMass())

    distance6 = (distance1+distance2+distance3+distance4)/4

	#a4_1,a4_2,a4_3 = a4.principalAxes()
	#b4_1,b4_2,b4_3 = b4.principalAxes()
   #     helix12_1,helix12_2,helix12_3 = helix12.principalAxes()
   #     helix21_1,helix21_2,helix21_3 = helix21.principalAxes()
   #     helix22_1,helix22_2,helix22_3 = helix22.principalAxes()

	#angle = math.degrees(math.acos(numpy.dot(a4_1,b4_1)))
       # angle2 = math.degrees(math.acos(numpy.dot(helix21_1,helix22_1)))

	#if angle > 90:
	#	angle = 180-angle
                
       # print "%6i %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f" % (ts.frame,rmsd0,rmsd1,rmsd2,distance1,distance2,angle1,angle2)

    #f.write('%7.3f %7.3f %7.3f %7.3f %7.3f\n' % (distance1,distance2,distance3,distance4,distance6))
    f.write('%7.3f\n' % (distance3))
	#f.write('%7.3f\n' % (distance2))
	#g.write('%7.3f\n' % angle)

f.close()
#g.close()

