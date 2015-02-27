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

#ahm_h1 = u.selectAtoms("segid A and resid 3 and name HO")
#ahm_h2 = u.selectAtoms("segid B and resid 9 and name HO")
##acc_1 = u.selectAtoms("segid A and resid 5 and name P")
##acc_2 = u.selectAtoms("segid A and resid 4 and name N7")
##acc_3 = u.selectAtoms("segid A and resid 4 and name N1")
#acc_1 = u.selectAtoms("segid A and resid 4 and name O4")
#acc_2 = u.selectAtoms("segid B and resid 10 and name O4")
#acc_5 = u.selectAtoms("segid B and resid 17 and name N3")
#acc_7 = u.selectAtoms("segid B and resid 15 and name N1")
bp_base11 = u.selectAtoms("segid A and resid 4 and name O4")
bp_base12 = u.selectAtoms("segid B and resid 9 and name H62")
bp_base13 = u.selectAtoms("segid A and resid 4 and name O2")
bp_base14 = u.selectAtoms("segid B and resid 9 and name H2")
bp_base21 = u.selectAtoms("segid B and resid 10 and name O4")
bp_base22 = u.selectAtoms("segid A and resid 3 and name H62")
bp_base23 = u.selectAtoms("segid B and resid 10 and name O2")
bp_base24 = u.selectAtoms("segid A and resid 3 and name H2")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    #distance1 = numpy.linalg.norm(ahm_h1.centerOfMass() - bp_base11.centerOfMass())
    #distance2 = numpy.linalg.norm(ahm_h2.centerOfMass() - bp_base21.centerOfMass())
    distance3 = numpy.linalg.norm(bp_base11.centerOfMass() - bp_base12.centerOfMass())
    distance4 = numpy.linalg.norm(bp_base13.centerOfMass() - bp_base14.centerOfMass())
    distance5 = numpy.linalg.norm(bp_base21.centerOfMass() - bp_base22.centerOfMass())
    distance6 = numpy.linalg.norm(bp_base23.centerOfMass() - bp_base24.centerOfMass())
    f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\n' % (distance3,distance4,distance5,distance6))
    #f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\n' % (distance3,distance4,distance1,distance2))
    #f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n' % (distance1,distance2,distance3,distance4,distance5,distance6,distance7,bp_distance))
f.close()

