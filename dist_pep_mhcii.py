#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]

u = Universe("init.pdb",my_traj)
v = Universe("init.pdb")

end = my_traj.find('.pdb')
fout_name = 'pepdist.dat'

# helices
a = u.selectAtoms("segid A and resid 50:73")
b = u.selectAtoms("segid B and resid 56:86")

# peptide
c = u.selectAtoms("segid C")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a.centerOfMass() - c.centerOfMass())
    distance2 = numpy.linalg.norm(b.centerOfMass() - c.centerOfMass())

    f.write('%7.3f %7.3f\n' % (distance1,distance2))

f.close()

