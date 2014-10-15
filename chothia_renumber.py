#! /usr/bin/python
# usage: chothia_renumber.py structure_file reference_file 
import warnings
from Bio import BiopythonWarning
from Bio.PDB import *
import sys,subprocess

my_file = sys.argv[1]
my_chothia_file = sys.argv[2]

parser = PDBParser(QUIET=True)
warnings.simplefilter('ignore', BiopythonWarning)

# get trajectory
my_traj = parser.get_structure('traj', my_file)

# get reference file
my_chothia_structure = parser.get_structure('chothia', my_chothia_file)

# get chain information
my_chothia_chain_H = my_chothia_structure[0]['H']
my_chothia_chain_L = my_chothia_structure[0]['L']

# get chothia residue information for both chains (H and L)
reslist_chothia_H = Selection.unfold_entities(my_chothia_chain_H, 'R')
reslist_chothia_L = Selection.unfold_entities(my_chothia_chain_L, 'R')

for my_model in my_traj:

    
    print "model"
    my_chain_H = my_model['H']
    my_chain_L = my_model['L']

    # get all residues
    reslist_H = Selection.unfold_entities(my_chain_H, 'R')
    reslist_L = Selection.unfold_entities(my_chain_L, 'R')
    
    # replace old numbering with chothia numbering from chothia pdb file

    print "L"
    for res_no in range(0, len(reslist_chothia_L)):
        
        reslist_L[res_no].id = reslist_chothia_L[res_no].id
    
    print "H"
    for res_no in range(0, len(reslist_chothia_H)):
        
        reslist_H[res_no].id = reslist_chothia_H[res_no].id


#for my_model in range(0,len(my_traj)):
#
#    my_chain_H = my_traj[my_model]['H']
#    my_chain_L = my_traj[my_model]['L']
#
#    # get all residues
#    reslist_H = Selection.unfold_entities(my_chain_H, 'R')
#    
#    reslist_L = Selection.unfold_entities(my_chain_L, 'R')
#    
#    # replace old numbering with chothia numbering from chothia pdb file
#    for res_no in range(0, len(reslist_chothia_L)):
#        
#        reslist_L[res_no].id = reslist_chothia_L[res_no].id
#    
#    for res_no in range(0, len(reslist_chothia_H)):
#        
#        reslist_H[res_no].id = reslist_chothia_H[res_no].id
        
# Save PDB file
end = my_file.find(".pdb")
my_out = my_file[0:end] + "_renum.pdb"
io = PDBIO()
io.set_structure(my_traj)
io.save(my_out)
#subprocess.call("fix_pdb_traj.py %s" % my_out, shell = True)

