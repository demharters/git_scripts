#! /usr/bin/python

import __main__
__main__.pymol_argv = ['pymol','-qc']
#__main__.pymol_argv = ['pymol','']
import sys,time,os
import pymol

pymol.finish_launching()
pymol.cmd.feedback("disable","all","actions")
pymol.cmd.feedback("disable","all","results")

sys.path.append("/home/scratch/software/Pymol-script-repo-master")
import my_elbow_angle_tcr

my_file = sys.argv[1]
chain1 = sys.argv[2]
chain2 = sys.argv[3]

end = my_file.find(".pdb")
file_name = my_file[0:end]

pymol.cmd.load(my_file,file_name)

#my_out_name = file_name + "_elbows.dat"
my_out_name = "elbows.dat"
my_out = open(my_out_name, "a") ## be careful to check whether file already existed!

angle = my_elbow_angle_tcr.elbow_angle(file_name,chain1,chain2)
my_out.write("%s\t%i\n" %(file_name,angle))

my_out.close()
pymol.cmd.quit()
