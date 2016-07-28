#!/usr/bin/python
import glob
import math, os, sys, re
import numpy
import time
from PDBUtils import PDBchain
from os.path import isfile, join
import matplotlib.pyplot as plt
import pickle

#Define your atomic cutoffs -- in Angstrom
O_CUTOFF = 6
H_CUTOFF = 6

#Define the pairs of atoms you want to have water between
def get_pairs():
	pairs = []
	#Example pair
	pairs.append((("A",4,'C5'),('A',5,'O2P')))

	return pairs

#Parse trajectory
def parse_trajectory(t_file):
	
	models = dict()
	curr_mod = -1
	mod_lines = []
	a_lines = []
	index = 0
	for line in open(t_file,'r'):
		line = line.strip()
		if 'ENDMDL' in  line:
			
			models[curr_mod] = numpy.array(zip(*sorted(models[curr_mod][0]))[1])
			#Condition on model number -- skip adding the model to the models object if the number has some property.
			#if curr_mod>10:
			#	continue
			
			
		if 'MODEL' in line:
			curr_mod = int(line.replace('MODEL',''))
			mod_lines = []
			a_lines = []
			models[curr_mod] = (mod_lines,a_lines)
			index = 0
		else:
			if 'ATOM' == line[0:4]:
				a_name = line[11:16].strip()
				
		        	x = float(line[30:38].strip()); y = float(line[38:46].strip()); z = float(line[46:54].strip())
				mod_lines.append((index,line))
				a_lines.append((index,a_name))
				index+=1
	
	return models

#Euclidean distance given differences in coordinates, x, y, z
def euclid(x,y,z):
	result = math.sqrt(x*x+y*y+z*z)
		
	return result

#Check water distance from the specified water -- oxygen -- NB: this returns only the first water under the cutoff!
def check_water_oxygen(h2o,water_atoms,residue,candidate_atom):
	#First Residue
	for a1 in residue.atoms:
		if (a1.typ ==candidate_atom):
			for a2 in water_atoms:
				if a2.typ != 'OH2':
					continue
				distance =  euclid(a1.coords[0]-a2.coords[0],a1.coords[1]-a2.coords[1],a1.coords[2]-a2.coords[2])
				if distance<O_CUTOFF:
					return [h2o],distance
	#If we didn't find anything below the specified cutoff, too bad.
	return [],-1

#Check water distance from the specified water -- only hydrogen -- 
def check_water_hydrogen(h2o,water_atoms,residue,candidate_atom):
	#First Residue
	for a1 in residue.atoms:
		if (a1.typ ==candidate_atom):
			for a2 in water_atoms:
				if a2.typ == 'OH2':
					continue
				distance =  euclid(a1.coords[0]-a2.coords[0],a1.coords[1]-a2.coords[1],a1.coords[2]-a2.coords[2])
				if distance<H_CUTOFF:
					return [h2o],distance
	#If we didn't find anything below the specified cutoff, too bad.
	return [],-1

#Calculate waters that are between the specified pairs
def hb_distances(residue_1,residue_2,waters,pair):
	
	#Hydrogen
	candidates_H = []
	#Oxygen
	candidates_O = []
	#For each water molecule in the structure, check its distance from the given pair
	for h2o in waters.residues:
			#Only consider waters, no counterions etc.
			if waters.residues[h2o].res_type !='TIP':
				continue
			#Hydrogen style
			candidate_1,dist_1 = check_water_hydrogen(h2o,waters.residues[h2o].atoms,residue_1,pair[0][2])
			candidate_2,dist_2 = check_water_hydrogen(h2o,waters.residues[h2o].atoms,residue_2,pair[1][2])
			if len(candidate_1)!=0 and len(candidate_2)!=0:
				candidates_H.append((h2o,dist_1,dist_2))
			#Oxygen style
			candidate_1,dist_1 = check_water_oxygen(h2o,waters.residues[h2o].atoms,residue_1,pair[0][2])
			candidate_2,dist_2 = check_water_oxygen(h2o,waters.residues[h2o].atoms,residue_2,pair[1][2])
			if len(candidate_1)!=0 and len(candidate_2)!=0:
				candidates_O.append((h2o,dist_1,dist_2))

	return candidates_H,candidates_O
	


#For a given run (source), given pairs of atoms (pairs) generate pickled file of water occupancy statistics
def compile_data(source,pickle_name,pairs):
	
	#Parse the trajectory file
	models = parse_trajectory(source)
	
	#Do we have an oxygen between a given atom pair?
	o_water_data = []
	#Do we have hydrogen between a given atom pair?
	full_water_data = []
	for model in models:
		
		cont_pair = ""
		chains = dict()
		
		waters = PDBchain(models[model],' ')
		#NB: assmuing there is only one chain in the structure.		
		nt_chain_A = PDBchain(models[model],'A')
		chains = {'A':nt_chain_A}
		for pair in pairs:
			
			a_1_chain = pair[0][0]
			a_2_chain = pair[1][0]

			a_1_residue = pair[0][1]
			a_2_residue = pair[1][1]
			
			dist_H,dist_O = hb_distances(chains[a_1_chain].residues[(a_1_residue,'')],chains[a_2_chain].residues[(a_2_residue,'')],waters,pair)
			o_water_data.append((model,pair,dist_O))
			full_water_data.append((model,pair,dist_H))
	plotdata = dict()
	#Are there any waters if we only look at oxygens?
	plotdata['o_waters'] = o_water_data
	#Are there any waters if we only look at hydrogens?
	plotdata['full_waters'] = full_water_data
	pickle.dump(plotdata,open(pickle_name,'wb'))

#Create pickled data for faster visualization
def pickle_data(input_trajectory,output_file):

	#get definitions of the pairs.
	pairs = get_pairs()
	#You need to have a pickles directory...
	if not os.path.exists('pickles'):
		os.mkdir('pickles')
	output_pickle = join('pickles',output_file)
	compile_data(input_trajectory,output_pickle,pairs)
	
	#Calculate statistics for non-methylated and methylated

#Main function
if __name__ == '__main__':
	import sys

	#Usage:
	#?>python WaterInteractions.py input_trajectory_file pickle_output_name
	#Result: results of the calculation will be stored in pickles/pickle_output_name

	#Get input trajectory
	trajectory = sys.argv[1]
	#How to call the output pickle file?
	output_file = sys.argv[2]
	pickle_data(trajectory,output_file)
	
	
			
	
    





