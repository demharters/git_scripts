#! /usr/bin/python

import __main__
__main__.pymol_argv = ['pymol','-qc']
#__main__.pymol_argv = ['pymol','']
import sys,time,os
import pymol
pymol.finish_launching()

sys.path.append("/home/scratch/software/Pymol-script-repo-master")
import my_elbow_angle_ab

my_file = sys.argv[1]
chain1 = sys.argv[2]
chain2 = sys.argv[3]
limit_h = sys.argv[4]
limit_l = sys.argv[5]

end = my_file.find(".pdb")
file_name = my_file[0:end]

pymol.cmd.load(my_file,file_name)
pymol.cmd.split_states(file_name)

states = pymol.cmd.get_object_list()

# skip first entry (first state is identical)
iterstates = iter(states)
next(iterstates)

my_out_name = file_name + "_elbows.dat"
my_out = open(my_out_name, "w")

for state in iterstates:
    angle = my_elbow_angle_ab.elbow_angle(state,chain1,chain2,limit_h,limit_l)
    my_out.write("%i\n" % angle)

my_out.close()
pymol.cmd.quit()
