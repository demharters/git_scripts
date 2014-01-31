#! /usr/bin/env python

from MDAnalysis import *
from MDAnalysis.analysis.align import rotation_matrix
import numpy
import math

u = Universe("init.pdb", "sampled.pos.pdb")
v = Universe("init.pdb")

mob1 = u.selectAtoms("segid A and resid 1:176")
ref1 = v.selectAtoms("segid A and resid 1:176")

f = open('rmsd_mda.dat','w')

for ts in u.trajectory:

        R, rmsd1 = rotation_matrix(mob1.coordinates(),ref1.coordinates())
      
	f.write('%7.3f\n' % rmsd1)

f.close() 
