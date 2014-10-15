#! /usr/bin/env python

from MDAnalysis import *
import numpy
import math
import sys

my_traj = sys.argv[1]
my_struc = sys.argv[2]
u = Universe(my_struc,my_traj)
v = Universe(my_struc)

end = my_traj.find('.pdb')
fout_name = my_traj[0:end] + '_EPOR_dist.dat'

# VHVL modules
a1 = v.selectAtoms("segid A and resid 1:123 or segid B and resid 129:237")
a2 = v.selectAtoms("segid B and resid 1:123 or segid A and resid 129:237")

f = open(fout_name,'w')

for ts in u.trajectory:
        
    A1 = a1.centerOfMass()
    A2 = a2.centerOfMass()
    #B1 = numpy.array([43.820398,-16.425929,-24.652051])
    #B2 = numpy.array([80.901129,-11.980774,-22.144420])
   
    # vhvl1 COM in DA10 according to pymol script centerOfMass.py
    #44.603714 -15.694621 -23.663141


    # difference between DA10 and DA5 COM
    #d1 = A1-B1
    #d2 = A2-B2
    d1 = numpy.array([2.92885319,-11.80629801,-3.65362573])
    d2 = numpy.array([-3.3932897,11.09212982,-1.43605291])
    
    # COM adjusted to fit DA5 COM
    A1_new = A1 + d1
    A2_new = A2 + d2
    
    # EPO-R C-termini
    C1_structure = v.selectAtoms("segid L and resid 215")
    C2_structure = v.selectAtoms("segid O and resid 215")
    
    # calculate vectors connecting EPO-R C-termini and COM of VHVL modules
    c1 = C1_structure.centerOfMass() - A1
    c2 = C2_structure.centerOfMass() - A2
    
    # DA10 vector
    c1_1 = numpy.array([-44.77196217,23.96018064,8.06541591])
    c1_2 = numpy.array([ 44.12922434,-25.01678919,8.39097504])
    
    # DA5 vector
    c2_1 = numpy.array([-39.46969247,-33.32577937,-42.38484518])
    c2_2 = numpy.array([44.19245885,29.43997329,-35.14102179])
   

    # calculate postion of DA10 EPO-R C-termini by adding c vectors to COM of VHVL modules
    C1_1_norm = A1 + c1_1
    C1_2_norm = A2 + c1_2
    # calculate postion of DA5 EPO-R C-termini by adding c vectors to COM of VHVL modules
    C2_1_norm = A1 + c2_1
    C2_2_norm = A2 + c2_2
    
    # calculate distance between the two EPO-R C-termini
    #distance = numpy.linalg.norm(C1_norm - C2_norm)
    
    distance1 = numpy.linalg.norm(C1_1_norm - C1_2_norm)
    distance2 = numpy.linalg.norm(C2_1_norm - C2_2_norm)
    
    #    # write DA5, DA10 and DA330 distances to file
    #f.write('%7.3f %7.3f\n' % (distance1,distance2))
    print C2_1_norm
    print C2_2_norm
    print c1,c2,c1_1,c1_2
    #

f.close()
