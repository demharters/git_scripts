#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]

u = Universe("init.pdb", my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_bg_glob_dist.dat'

bg1 = u.selectAtoms("segid A and resid 50:51")
bg2 = u.selectAtoms("segid B and resid 85:86")
glob1 = u.selectAtoms("segid A and resid 84:182")
glob2 = u.selectAtoms("segid B and resid 95:190")

with open(fout_name, "w") as f:

    for ts in u.trajectory:
            
            bg_distance = numpy.linalg.norm(bg1.centerOfMass() - bg2.centerOfMass())
            glob_distance = numpy.linalg.norm(glob1.centerOfMass() - glob2.centerOfMass())
    
            f.write('%7.3f\t%7.3f\n' % (bg_distance, glob_distance))
