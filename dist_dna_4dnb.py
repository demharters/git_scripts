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

#a1 = u.selectAtoms("segid A and resid 3 and name HO")
#b1 = u.selectAtoms("segid A and resid 4 and name O4")
#
#a2 = u.selectAtoms("segid B and resid 9 and name HO")
#b2 = u.selectAtoms("segid B and resid 10 and name O4")

#a1 = u.selectAtoms("segid A and resid 6 and name HO")
#b1 = u.selectAtoms("segid A and resid 7 and name O4")
#
#a2 = u.selectAtoms("segid B and resid 18 and name HO")
#b2 = u.selectAtoms("segid B and resid 19 and name O4")

#a1 = u.selectAtoms("segid A and resid 3 and name C5M")
#b1 = u.selectAtoms("segid A and resid 4 and name O4")
#
#a2 = u.selectAtoms("segid B and resid 9 and name C5M")
#b2 = u.selectAtoms("segid B and resid 10 and name O4")

#a1 = u.selectAtoms("bynum 88")
#b1 = u.selectAtoms("bynum 123")
#a1 = u.selectAtoms("segid A and resid 3 and name HO")
#b1 = u.selectAtoms("segid A and resid 4 and name O4")
#a2 = u.selectAtoms("segid A and resid 3 and name H62")
#b2 = u.selectAtoms("segid B and resid 10 and name O4")
#
##a3 = u.selectAtoms("segid B and resid 9 and atom.number 283")
#a3 = u.selectAtoms("bynum 283")
#b3 = u.selectAtoms("segid B and resid 10 and name O4")
#a4 = u.selectAtoms("segid B and resid 9 and name H62")
#b4 = u.selectAtoms("segid A and resid 4 and name O4")
#da1 = u.selectAtoms("segid A and resid 2")

ahm_h = u.selectAtoms("segid A and resid 5 and name HO")
#acc_1 = u.selectAtoms("segid A and resid 5 and name P")
acc_2 = u.selectAtoms("segid A and resid 4 and name N7")
#acc_3 = u.selectAtoms("segid A and resid 4 and name N1")
acc_4 = u.selectAtoms("segid A and resid 6 and name O4")
#acc_5 = u.selectAtoms("segid B and resid 17 and name N3")
acc_6 = u.selectAtoms("segid B and resid 16 and name N3")
#acc_7 = u.selectAtoms("segid B and resid 15 and name N1")
bp_base1 = u.selectAtoms("segid A and resid 5 and name N1")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    #distance1 = numpy.linalg.norm(a1.centerOfMass() - b1.centerOfMass())
    #distance2 = numpy.linalg.norm(a2.centerOfMass() - b2.centerOfMass())
    #distance3 = numpy.linalg.norm(a3.centerOfMass() - b3.centerOfMass())
    #distance4 = numpy.linalg.norm(a4.centerOfMass() - b4.centerOfMass())
    #f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\n' % (distance1,distance2,distance3,distance4))

    #distance1 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_1.centerOfMass())
    distance2 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_2.centerOfMass())
    #distance3 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_3.centerOfMass())
    distance4 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_4.centerOfMass())
    #distance5 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_5.centerOfMass())
    #distance6 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_6.centerOfMass())
    #distance7 = numpy.linalg.norm(ahm_h.centerOfMass() - acc_7.centerOfMass())
    bp_distance = numpy.linalg.norm(acc_6.centerOfMass() - bp_base1.centerOfMass())
    f.write('%7.3f\t%7.3f\t%7.3f\n' % (distance2,distance4,bp_distance))
    #f.write('%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\t%7.3f\n' % (distance1,distance2,distance3,distance4,distance5,distance6,distance7,bp_distance))
f.close()

