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

#u = Universe("init.pdb",my_traj)
#v = Universe("init.pdb")

u = Universe(my_traj,my_traj)
v = Universe(my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_dist.dat'

o6_1 = u.selectAtoms("segid A and resid 9 and name O6")
chm_1 = u.selectAtoms("segid A and resid 8 and name HO")
o6_2 = u.selectAtoms("segid B and resid 19 and name O6")
chm_2 = u.selectAtoms("segid B and resid 18 and name HO")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    distance1 = numpy.linalg.norm(o6_1.centerOfMass() - chm_1.centerOfMass())
    distance2 = numpy.linalg.norm(o6_2.centerOfMass() - chm_2.centerOfMass())
    
    f.write('%7.3f\t%7.3f\n' % (distance1,distance2))
f.close()

