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

#o6_1 = u.selectAtoms("segid A and resid 9 and name O6")
#chm_1 = u.selectAtoms("segid A and resid 8 and name HO")
#o6_2 = u.selectAtoms("segid B and resid 19 and name O6")
#chm_2 = u.selectAtoms("segid B and resid 18 and name HO")

A9_1 = u.selectAtoms("segid A and resid 9 and name O4'")
A8_1 = u.selectAtoms("segid A and resid 8 and name O4'")
A9_2 = u.selectAtoms("segid A and resid 9 and name C1'")
A8_2 = u.selectAtoms("segid A and resid 8 and name C1'")
A9_3 = u.selectAtoms("segid A and resid 9 and name C2'")
A8_3 = u.selectAtoms("segid A and resid 8 and name C2'")
A9_4 = u.selectAtoms("segid A and resid 9 and name C3'")
A8_4 = u.selectAtoms("segid A and resid 8 and name C3'")
A9_5 = u.selectAtoms("segid A and resid 9 and name C4'")
A8_5 = u.selectAtoms("segid A and resid 8 and name C4'")

A9_ring = u.selectAtoms("segid A and resid 9 and (name O4' or name C1' or name C2' or name C3' or name C4')")
A8_ring = u.selectAtoms("segid A and resid 8 and (name O4' or name C1' or name C2' or name C3' or name C4')")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    #distance1 = numpy.linalg.norm(o6_1.centerOfMass() - chm_1.centerOfMass())
    #distance2 = numpy.linalg.norm(o6_2.centerOfMass() - chm_2.centerOfMass())
    
    distance1 = numpy.linalg.norm(A9_1.centerOfMass() - A8_1.centerOfMass())
    distance2 = numpy.linalg.norm(A9_2.centerOfMass() - A8_2.centerOfMass())
    distance3 = numpy.linalg.norm(A9_3.centerOfMass() - A8_3.centerOfMass())
    distance4 = numpy.linalg.norm(A9_4.centerOfMass() - A8_4.centerOfMass())
    distance5 = numpy.linalg.norm(A9_5.centerOfMass() - A8_5.centerOfMass())
    distance6 = numpy.linalg.norm(A9_ring.centerOfMass() - A8_ring.centerOfMass())
    
    f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n' % (distance1,distance2,distance3,distance4,distance5,distance6))
f.close()

