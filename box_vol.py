import os
import sys
import MDAnalysis as md
#import my_plotter as myplot
import math
import numpy as np
import scipy as sp
import numpy.linalg
import itertools
import matplotlib.pyplot as plt
import seaborn as sns

##########################################################################################################################################

#md.core.flags['use_pbc'] = True

###########################################################################################################################################

color_bndr = ['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']

for file in os.listdir("./../../cis_flat/12_8_v04/"):
    if file.endswith(".psf"):
        PSF_file = './../../cis_flat/12_8_v04/'+file

sns.set_style("white", {'xtick.direction': u'in','ytick.direction': u'in'})
sns.set_context("poster")

#traj_list = []
#for file in os.listdir("../data"):
#    if file.endswith(".dcd"):
#        traj_list.append('../data/'+file)

traj_list = ['./../../cis_flat/12_8_v04/T298/restart07/Ndc_Nte9cis_flat.dcd']
print traj_list
print PSF_file, traj_list

u            = md.Universe(PSF_file, traj_list)

time         = np.array([ts.time for ts in u.trajectory])
time         = np.multiply ( time, 0.001 )
box_vol      = np.array([ts.dimensions[0:3] for ts in u.trajectory])

#########################################################################################################################################

fig, ax1 = plt.subplots(figsize=(6,5))
plt.interactive(False)
ax1.set_ylabel(r'$L$', fontsize=24)
ax1.set_xlabel(r'$t (ns)$', fontsize=24)
ax1.plot(time, box_vol[:,0],'o', label = r'$L$', color = color_bndr[0], markersize = 6.0)
ax2 = ax1.twinx()
ax2.set_ylabel(r'$B,H$', fontsize=24)
ax2.plot(time, box_vol[:,1],'o', label = r'$B$', color = color_bndr[1], markersize = 6.0)
ax2.plot(time, box_vol[:,2],'o', label = r'$H$', color = color_bndr[2], markersize = 6.0)
ax1.autoscale()
ax1.legend(loc='best')
ax2.legend(loc='best')
plt.tight_layout()
plt.savefig('bx_volT298.pdf')

nylayer = 4
nzlayer = 12

fig, ax1 = plt.subplots(figsize=(6,5))
plt.interactive(False)
ax1.set_ylabel(r'$c$', fontsize=24)
ax1.set_xlabel(r'$t (ns)$', fontsize=24)
ax1.plot(time, np.divide(box_vol[:,1], nylayer) ,'o', label = r'$c$', color = color_bndr[0], markersize = 6.0)
ax2 = ax1.twinx()
ax2.set_ylabel(r'$a$', fontsize=24)
ax2.plot(time, np.divide(box_vol[:,2], nzlayer),'o', label = r'$a$', color = color_bndr[1], markersize = 6.0)
ax1.autoscale()
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
plt.tight_layout()
plt.savefig('lat_paramT298.pdf')

#########################################################################################################################################
