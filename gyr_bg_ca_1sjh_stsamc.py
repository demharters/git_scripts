#!/usr/bin/python

from MDAnalysis import *
import numpy.linalg

u = Universe("temp.pos.pdb")
calphas = u.selectAtoms('name CA')
bg = calphas.selectAtoms('segid A and resid 4:77 or segid B and resid 1:90')

g = open('rgyr_bg','w')

for ts in u.trajectory:

	rgyr_bg = bg.radiusOfGyration()

	g.write('%f\n' % (rgyr_bg,))

g.close()
