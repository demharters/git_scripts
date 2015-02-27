#! /usr/bin/python

import __main__
__main__.pymol_argv = ['pymol','-qc']
#__main__.pymol_argv = ['pymol','']
import sys,time,os
import pymol
pymol.finish_launching()

my_file = sys.argv[1]

end = my_file.find(".pdb")
file_name = my_file[0:end]

pymol.cmd.load(my_file,file_name)
#pymol.cmd.split_states(file_name)
#pymol.cmd.select("res 18 & name N1 | res 18 & name N3 | res 18 & name C4")
pymol.cmd.intra_fit("res 18 & name N1 | res 18 & name N3 | res 18 & name C4")

my_out_name = file_name + "_aligned_OH.pdb"

pymol.cmd.save(my_out_name, file_name, state=0)
pymol.cmd.quit()
