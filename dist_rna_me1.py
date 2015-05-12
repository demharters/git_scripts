#! /usr/bin/env python

#import __main__
#__main__.pymol_argv = ['pymol','-qc']
##__main__.pymol_argv = ['pymol','']
#import sys,time,os
#import pymol
#
#pymol.finish_launching()
#pymol.cmd.feedback("disable","all","actions")
#pymol.cmd.feedback("disable","all","results")
#sys.path.append("/home/scratch/software/Pymol-script-repo-master")

import MDAnalysis
from MDAnalysis import *
from MDAnalysis.analysis.distances import *
import numpy
import math
import sys

my_traj = sys.argv[1]

#pymol.cmd.load(my_traj,my_traj[:-4])
#pymol.cmd.split_states(my_traj[:-4])
#states = pymol.cmd.get_object_list()
#
#
#for state in states:
#
#    #dih_angle = pymol.cmd.get_dihedral("%s//A/2/N1"%state, "%s//A/2/C6"%state, "%s//A/2/N6"%state, "%s//A/2/H61"%state) 
#    dih_angle = pymol.cmd.get_dihedral("%s//A/2/N1"%state, "%s//A/2/C6"%state, "%s//A/2/N6"%state, "%s//A/2/H61"%state) 
#    print dih_angle

u = Universe("init.pdb",my_traj)
v = Universe("init.pdb")

u = Universe(my_traj,my_traj)
v = Universe(my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_dist.dat'

a1 = u.selectAtoms("segid A and resid 2 and name H51")
b1 = u.selectAtoms("segid A and resid 1 and name O6")
#da1 = u.selectAtoms("segid A and resid 2")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
    f.write('%7.3f\n' % (distance1))

f.close()

