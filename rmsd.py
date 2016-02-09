#! /usr/bin/python

from MDAnalysis import *
import sys
from MDAnalysis.analysis.rms import *


my_traj = sys.argv[1]
u = Universe(my_traj,my_traj)
bb = u.selectAtoms('all and not resname TIP')

A = bb.coordinates()  # coordinates of first frame

my_openFile = open("rmsd.dat","w")

for ts in u.trajectory:

    B = bb.coordinates()  # coordinates of last frame
    myRmsd = rmsd(A,B)
    my_openFile.write("%s\n"%myRmsd)

my_openFile.close()

    





