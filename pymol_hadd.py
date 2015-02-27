#! /usr/bin/python

import __main__
__main__.pymol_argv = ['pymol','-qc']
#__main__.pymol_argv = ['pymol','']
import sys,time,os
sys.path.append("/home/scratch/software/pymol/")
import pymol
pymol.finish_launching()

my_file = sys.argv[1]

end = my_file.find(".pdb")
file_name = my_file[0:end]

pymol.cmd.load(my_file,file_name)
pymol.cmd.h_add(file_name)

my_out_name = file_name + "_hadd.pdb"

pymol.cmd.save(my_out_name, file_name)
pymol.cmd.quit()
