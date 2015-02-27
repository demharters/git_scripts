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

#import MDAnalysis
from MDAnalysis import *
import MDAnalysis.analysis.contacts
import numpy
import math
import sys

my_traj = sys.argv[1]

#u = Universe("init.pdb",my_traj)
#v = Universe("init.pdb")

u = Universe(my_traj,my_traj)
v = Universe(my_traj)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_dist.dat'

sel_epimod1 = "segid A and resid 3 and name HO"
sel_epimod1env = "segid A and (resid 4 or resid 2)"

a1 = u.selectAtoms(sel_epimod1)
b1 = u.selectAtoms(sel_epimod1env)

f = open(fout_name,'w')

#for ts in u.trajectory:
#       
#    CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(u, selection=(sel_epimod1, sel_epimod1env), refgroup=(a1,b1), radius=6.0, outfile="CA_epimod1.dat")
#    f.write('%7.3f\n' % (CA1))

#f.close()


CA1 = MDAnalysis.analysis.contacts.ContactAnalysis1(u, selection=(sel_epimod1, sel_epimod1env), refgroup=(a1,b1), radius=6.0, outfile="CA_epimod1.dat")

CA1.run(force=True)

CA1.plot(filename="CA_epimod1.pdf",linewidth=3, color="blue")
