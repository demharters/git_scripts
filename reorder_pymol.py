#! /usr/bin/python

import __main__
__main__.pymol_argv = ['pymol','-qc']
#__main__.pymol_argv = ['pymol','']
import sys,time,os
import pymol
pymol.finish_launching()

sys.path.append("/home/scratch/software/Pymol-script-repo-master")
import my_elbow_angle_tcr

my_file = sys.argv[1]

end = my_file.find(".pdb")
file_name = my_file[0:end]

pymol.cmd.load(my_file,file_name)

pymol.cmd.save("init_reordered.pdb",file_name)

pymol.cmd.quit()
