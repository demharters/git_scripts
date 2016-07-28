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
fileName = my_file[0:end]
newFileName = fileName + "_reordered.pdb"

pymol.cmd.load(my_file,fileName)

pymol.cmd.save(newFileName,fileName,state=0)

pymol.cmd.quit()
