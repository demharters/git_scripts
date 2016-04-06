#! /usr/bin/python

import sys
from MDAnalysis import *
import numpy

my_topol = sys.argv[1]
my_traj = sys.argv[2]

u = Universe(my_topol,my_traj)
#v = Universe(my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_dist.dat'
print "%s" % fout_name

myMMcpn = u.selectAtoms("segid *")
myCOM = myMMcpn.centerOfMass()

myMMcpn_tip1 = u.selectAtoms("segid A and resid 236")
myMMcpn_tip2 = u.selectAtoms("segid B and resid 236")
myMMcpn_tip3 = u.selectAtoms("segid C and resid 236")
myMMcpn_tip4 = u.selectAtoms("segid D and resid 236")
myMMcpn_tip5 = u.selectAtoms("segid D and resid 236")
myMMcpn_tip6 = u.selectAtoms("segid E and resid 236")
myMMcpn_tip7 = u.selectAtoms("segid F and resid 236")
myMMcpn_tip8 = u.selectAtoms("segid G and resid 236")
myMMcpn_tip9 = u.selectAtoms("segid H and resid 236")
myMMcpn_tip10 = u.selectAtoms("segid I and resid 236")
myMMcpn_tip11 = u.selectAtoms("segid K and resid 236")
myMMcpn_tip12 = u.selectAtoms("segid L and resid 236")
myMMcpn_tip13 = u.selectAtoms("segid M and resid 236")
myMMcpn_tip14 = u.selectAtoms("segid N and resid 236")
myMMcpn_tip15 = u.selectAtoms("segid O and resid 236")
myMMcpn_tip16 = u.selectAtoms("segid P and resid 236")


f = open(fout_name,'w')

for ts in u.trajectory:
    distance1 = numpy.linalg.norm(myCOM - myMMcpn_tip1.centerOfMass())
    distance2 = numpy.linalg.norm(myCOM - myMMcpn_tip2.centerOfMass())
    distance3 = numpy.linalg.norm(myCOM - myMMcpn_tip3.centerOfMass())
    distance4 = numpy.linalg.norm(myCOM - myMMcpn_tip4.centerOfMass())
    distance5 = numpy.linalg.norm(myCOM - myMMcpn_tip5.centerOfMass())
    distance6 = numpy.linalg.norm(myCOM - myMMcpn_tip6.centerOfMass())
    distance7 = numpy.linalg.norm(myCOM - myMMcpn_tip7.centerOfMass())
    distance8 = numpy.linalg.norm(myCOM - myMMcpn_tip8.centerOfMass())
    distance9 = numpy.linalg.norm(myCOM - myMMcpn_tip9.centerOfMass())
    distance10 = numpy.linalg.norm(myCOM - myMMcpn_tip10.centerOfMass())
    distance11 = numpy.linalg.norm(myCOM - myMMcpn_tip11.centerOfMass())
    distance12 = numpy.linalg.norm(myCOM - myMMcpn_tip12.centerOfMass())
    distance13 = numpy.linalg.norm(myCOM - myMMcpn_tip13.centerOfMass())
    distance14 = numpy.linalg.norm(myCOM - myMMcpn_tip14.centerOfMass())
    distance15 = numpy.linalg.norm(myCOM - myMMcpn_tip15.centerOfMass())
    distance16 = numpy.linalg.norm(myCOM - myMMcpn_tip16.centerOfMass())

    meanDistance1 = numpy.mean([distance1,distance2,distance3,distance4,distance5,distance6,distance7,distance8])
    meanDistance2 = numpy.mean([distance9,distance10,distance11,distance12,distance13,distance14,distance15,distance16])

    f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f%7.3f\t%7.3f\t%7.3f\t%7.3f%7.3f%7.3f\n' \
    	% (distance1,distance2,distance3,distance4,distance5,distance6,distance7,distance8,distance9,distance10,distance11,distance12,distance13,distance14,distance15,distance16,meanDistance1,meanDistance2))
