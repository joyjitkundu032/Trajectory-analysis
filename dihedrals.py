import os
import sys
import MDAnalysis as md
#import my_plotter as myplot
import matplotlib.pyplot as plt
import math
import numpy as np
import scipy as sp
import numpy.linalg
import itertools
#import backbone_dihedrals as bb_dihedrals

import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(labelsize=20, pad = 10, length=8, width=2)

def dihedral(p0, p1, p2, p3):
    
    "Return angle [-pi,pi] formed by vertices p0-p1-p2-p3."

    b1  = np.subtract(p1, p0)  
    b2  = np.subtract(p2, p1)
    b3  = np.subtract(p3, p2)

    yy  = np.dot(np.cross(np.cross(b1,b2), np.cross(b2,b3)), np.divide(b2,np.linalg.norm(b2)))
    xx  = np.dot(np.cross(b1,b2), np.cross(b2,b3))

    ang = math.atan2 (yy, xx)

    return ang


def backbone_dihedral_angles (no_residues, bb_Nitro, bb_Calph, bb_CCarb):

    phi = np.zeros(no_residues-1); psi = np.zeros(no_residues-1); ome = np.zeros(no_residues-1);

    for ires in np.arange(no_residues-1):

        phi[ires]     = dihedral (bb_CCarb[ires,:],   bb_Nitro[ires+1,:], bb_Calph[ires+1,:], bb_CCarb[ires+1,:])
        psi[ires]     = dihedral (bb_Nitro[ires+1,:], bb_Calph[ires+1,:], bb_CCarb[ires+1,:], bb_Nitro[ires+2,:])
        ome[ires]     = dihedral (bb_Calph[ires,:],   bb_CCarb[ires,:],   bb_Nitro[ires+1,:], bb_Calph[ires+1,:])

    return phi, psi, ome


def calpha_distances (no_residues, bb_Nitro, bb_Calph, bb_CCarb):

    calpha_dist = np.zeros(no_residues-1);

    for ires in np.arange(no_residues-1):

        p0  = bb_Calph[ires]; p1 = bb_Calph[ires+1]
        b1  = np.subtract(p1, p0)  
        calpha_dist[ires] = np.linalg.norm(b1)
        
    return calpha_dist
##########################################################################################################################################
#md.core.flags['use_pbc'] = True
##########################################################################################################################################

for file in os.listdir("cis_flat/12_8_v04/"):
    if file.endswith(".psf"):
        PSF_file = 'cis_flat/12_8_v04/'+file

sns.set_style("white", {'xtick.direction': u'in','ytick.direction': u'in'})
sns.set_context("poster")

#traj_list = []
#for file in os.listdir("../data"):
#    if file.endswith(".dcd"):
#        traj_list.append('../data/'+file)

traj_list = ['cis_flat/12_8_v04/T298/restart08/Ndc_Nte9cis_flat.dcd','cis_flat/12_8_v04/T298/restart09/Ndc_Nte9cis_flat.dcd','cis_flat/12_8_v04/T298/restart10/Ndc_Nte9cis_flat.dcd','cis_flat/12_8_v04/T298/restart11/Ndc_Nte9cis_flat.dcd']
print traj_list
print PSF_file, traj_list

#traj_list    = ['../data/lam_r6.dcd','../data/lam_r7.dcd']
u            = md.Universe(PSF_file, traj_list)

nresidue     = 18
nchains      = 96

bt_Nitro = []; bt_Calph = []; bt_Ccarb = [];

traj_length  = len(u.trajectory)
calpha_dist  = np.zeros((nresidue-1, nchains, traj_length))

psi_traj_U   = np.zeros((nresidue-1, nchains, traj_length)); 
phi_traj_U   = np.zeros((nresidue-1, nchains, traj_length)); 
ome_traj_U   = np.zeros((nresidue-1, nchains, traj_length)); 

#phi_U        = np.zeros(0); psi_U        = np.zeros(0); ome_U        = np.zeros(0);

for ichain in np.arange(nchains):   
    segid = 'U'+str(ichain+1)
        
    bt_Nitro.append ( u.select_atoms ( "segid "+segid+" and (name N or name NT)" )) ;
    bt_Calph.append ( u.select_atoms ( "segid "+segid+" and name CA" )) ;
    bt_Ccarb.append ( u.select_atoms ( "segid "+segid+" and name C"  )) ;

index = 0

print traj_length

for ts in u.trajectory:     # iterate through all frames  

    for ichain in np.arange(nchains):                           

        bb_Nitro = bt_Nitro[ichain]; bb_Calph = bt_Calph[ichain]; bb_Ccarb = bt_Ccarb[ichain];
        phi_curr_traj_U, psi_curr_traj_U, ome_curr_traj_U = backbone_dihedral_angles (nresidue, bb_Nitro.positions, bb_Calph.positions, bb_Ccarb.positions)
        
        phi_traj_U[:, ichain, index] = phi_curr_traj_U; 
        psi_traj_U[:, ichain, index] = psi_curr_traj_U; 
        ome_traj_U[:, ichain, index] = ome_curr_traj_U; 

        calpha_dist_curr = calpha_distances(nresidue, bb_Nitro.positions, bb_Calph.positions, bb_Ccarb.positions)
        calpha_dist[:,ichain, index] = calpha_dist_curr

    index += 1

#myplot.make_histo ( calpha_dist, r'$r$', r'$C_{\alpha} - C_{\alpha}$', 'calpha_hist.pdf')

plt.figure(figsize=(6,5))
plt.interactive(False)
plt.ylabel( r'$\psi$', fontsize=24)
plt.xlabel( r'$\phi$', fontsize=24)
plt.xlim([-math.pi, math.pi])
plt.ylim([-math.pi, math.pi])
lbls = [-math.pi, -math.pi/2, 0.0, math.pi/2, math.pi]
plt.xticks(lbls,[r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"])
plt.yticks(lbls,[r"$-\pi$", r"$-\pi/2$", "$0$", r"$\pi/2$", r"$\pi$"])
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=4, mode="expand", borderaxespad=0.)
cmap = sns.cubehelix_palette(light = 1, as_cmap=True)
plt.hist2d ( phi_traj_U.flatten(), psi_traj_U.flatten(), bins = 180, cmap=cmap)
plt.colorbar()
plt.tight_layout()
plt.savefig('hist_ramachandran.pdf', bbox_inches="tight")

#########################################################################################################################################


