#!check reader
import sys
from pdb_reader import PDB_reader
import math
import numpy as np
import numpy.linalg
import itertools
import random
from matplotlib.widgets import Slider 
from matplotlib.ticker import LogFormatter 

import matplotlib.pyplot as plt
import matplotlib.colors
from scipy.interpolate import spline
#import seaborn as sns

#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.tick_params(labelsize=20, pad = 10, length=8, width=2)

file = sys.argv
pdb_input=[]

with open(file[1], "r") as myfile:
	for i, line in enumerate(myfile):
		pdb_input.append(line.strip('\n'))
		
print(myfile.name)
		
atom_name, residue_name, residue_num, chainID, atom_pos, segid, dimensions = PDB_reader.read_pdb(pdb_input)

#carbon_list=['C','CA','CA5','CB','CG','CD1','CD','CD2','CE1','CE2','CZ','CE','CH','CT','CI','CK','CG1','CG2','CZ1','C03']
carbon_list=['CA','C']
#nitrogen_list=['N','NG','NZ']
nitrogen_list=['N']
oxygen_list=['O']
hydrogen_list=['H']
nitrogen_segments=[]
oxygen_segments=[]
carbons=[]
nitrogens=[]
oxygens=[]
hydrogens=[]
i=0
distances=[]
for at in atom_name:
	if at in carbon_list:
		carbons.append(atom_pos[i])
	elif at in nitrogen_list:
		nitrogens.append(atom_pos[i])
		nitrogen_segments.append(segid[i])
	elif at in oxygen_list:
		oxygens.append(atom_pos[i])
		oxygen_segments.append(segid[i])
	elif at in hydrogen_list:
		hydrogens.append(atom_pos[i])
	i+=1

#print ('Carbons: ' + str(len(carbons)), 'Oxygens: ' + str(len(oxygens)), 'Nitrogens: ' + str(len(nitrogens)))	

#carbons_x = [x[0] for x in carbons]
#carbons_y = [y[1] for y in carbons]
#carbons_z = [z[2] for z in carbons]

nitrogens_x = [x[0] for x in nitrogens]
nitrogens_y = [y[1] for y in nitrogens]
nitrogens_z = [z[2] for z in nitrogens]


oxygens_x = [x[0] for x in oxygens]
oxygens_y = [y[1] for y in oxygens]
oxygens_z = [z[2] for z in oxygens]

def correct_distance(u1, u2, ind):
	dm = float(dimensions[ind])
	dis = u1-u2
	if dis > dm / 2:
		dis = dis-dm
	elif dis <= -dm /2:
		dis=dis+dm
	return dis	

def calculate_distance(x1,x2,y1,y2,z1,z2):
	dx=correct_distance(x1,x2,0)
	dy=correct_distance(y1,y2,1)
	dz=correct_distance(z1,z2,2)
	distance = math.sqrt(math.pow(dx,2) + math.pow(dy,2) + math.pow(dz,2))
	return distance
m=0
for N1 in nitrogens:
	distance_min = 1000
	n=0
	for O2 in oxygens:
		if nitrogen_segments[m] != oxygen_segments[n]:
			pass
		else:
			distance = calculate_distance(N1[0],O2[0],N1[1],O2[1],N1[2],O2[2])
			if distance < distance_min:
				distance_min = distance
			else:
				pass
		n+=1
	#print(m,nitrogen_segments[m],distance_min)
	m+=1
	distances.append(distance_min)

data = open(str(myfile.name) + '_backbone_distance_short_NO_intra.txt', 'w')
i=0
while i < len(distances):
		data.write(str(distances[i])+'\n')
		i+=1
	
data.close()
