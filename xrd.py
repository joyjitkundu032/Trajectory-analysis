import os
import sys
import MDAnalysis as md
import math
import numpy as np
import numpy.linalg
import itertools

from matplotlib.ticker import LogFormatter , FormatStrFormatter

import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(labelsize=20, pad = 10, length=8, width=2)

sns.set_style("white", {'xtick.direction': u'in','ytick.direction': u'in'})
sns.set_context("poster")

##########################################################################################################################################

#md.core.flags['use_pbc'] = True

##########################################################################################################################################

def make_image(image, xlbl, ylbl, filename):
    
    fig, ax1 = plt.subplots(figsize=(6,5))
    plt.interactive(False)

    plt.ylabel ( ylbl, fontsize=24 )
    plt.xlabel ( xlbl, fontsize=24 )

    #img = plt.imshow(np.log(image), interpolation='None', origin='lower', cmap = "RdBu")

    img       = plt.imshow(np.transpose(np.log(image)), interpolation='None', origin='lower', cmap = "RdBu", norm=matplotlib.colors.LogNorm(), aspect='auto')
    formatter = LogFormatter(10, labelOnlyBase=False) 
    cb        = plt.colorbar()
    ax1.autoscale()

    npoints_x = image.shape[0]
    npoints_y = image.shape[1]

    xlbls = [0, npoints_x/4, npoints_x/2, (3*npoints_x)/4, npoints_x]
    ylbls = [0, npoints_y/4, npoints_y/2, (3*npoints_y)/4, npoints_y]

    tickvals_x = np.around ( np.array ([qx[nresx-nres_p_x], qx [ nresx -nres_p_x/2], qx[nresx], qx[nresx + nres_p_x/2], qx[nres_p_x + nresx]]), decimals = 3)
    tickvals_y = np.around ( np.array ([qy[nresy-nres_p_y], qy [ nresy -nres_p_y/2], qy[nresy], qy[nresy + nres_p_y/2], qy[nres_p_y + nresy]]), decimals = 3)

    plt.xticks ( xlbls,tickvals_x ) 
    plt.yticks ( ylbls,tickvals_y ) 

    plt.tight_layout()
    plt.savefig(filename)
    
    return

def make_plot(x,y, xlbl, ylbl, filename):
    
    fig, ax1 = plt.subplots(figsize=(6,5))
    plt.interactive(False)
    plt.ylabel(ylbl, fontsize=24)
    plt.xlabel(xlbl, fontsize=24)
    plt.semilogy(x, y, linewidth = 1)
    ax1.autoscale()
    plt.tight_layout()
    plt.savefig(filename)
    #plt.show()
    return
                      

PSF_file     = 'solid_planar_6_16_3_PA02.psf'
traj_list    = ['Ndc_Nte9cis.dcd']
nres         = 80
u            = md.Universe(PSF_file, traj_list)

nresx        = nres
nresy        = nres
nresz	     = 3*nres

carbons      = u.select_atoms ( "(type C or type CA or type CAY or type CY or type CA5 or type CB or type CG or type CD or type CE or type CZ or type CH or type CT or type CI or type CK )" ) ;
hydrogens    = u.select_atoms ( "(type HA1 or type HA2 or type HA3 or type HA4 or type HB1 or type HB2 or type HG1 or type HG2 or type HD1 or type HD2 or type HE1 or type HE2 or type HZ1 or type HZ2 or type HH1 or type HH2 or type HT3 or type HT4 or type HI1 or type HI2 or type HK1 or type HK2 or type HK3 or type HY1 or type HY2 or type HY3 or type HN)" ) ;
oxygens      = u.select_atoms ( "(type O or type OY or type OG or type OZ or type OI )" ) ;
nitrogens    = u.select_atoms ( "(type N)" ) ;

n_qx         = np.linspace ( -nresx, nresx, num = (2*nresx) + 1, endpoint = True, dtype=np.float64 )
n_qy         = np.linspace ( -nresy, nresy, num = (2*nresy) + 1, endpoint = True, dtype=np.float64 )
n_qz         = np.linspace ( -nresz, nresz, num = (2*nresz) + 1, endpoint = True, dtype=np.float64 )

I_save       = np.zeros((len(n_qx), len(n_qy),len(n_qz))) 
Ig           = (len(nitrogens)*49) + (len(carbons)*36) + (len(oxygens)*64) + (len(hydrogens))

all_atoms   = u.select_atoms ( "(type C or type CT2 or type CA or type CC or type H or type HB or type HA or type HC or type HP or type O or type OC or type NH1 or type NH2 or type NH3)" ) ;

atom_numbr  = np.multiply ( all_atoms.masses, 0.50)
atom_numbr[ atom_numbr==0.50395] = 1.0  # Make change for hydrogen

#########################################################################################################################################

ntrj  = 0
index = 0

for ts in u.trajectory:     # iterate through all frames

    curr_time = ts.time

    if (index > 8):
        
        print curr_time

        bx = ts.dimensions[0:3]

        I  = np.zeros((len(n_qx), len(n_qy),len(n_qz)), dtype=np.complex128)
        qx = np.multiply ( np.copy ( n_qx ), (2.0*math.pi)/bx[0] )
        qy = np.multiply ( np.copy ( n_qy ), (2.0*math.pi)/bx[1] )
	qz = np.multiply ( np.copy ( n_qz ), (2.0*math.pi)/bx[2] )

        pos_x   = all_atoms.positions[:,0]
        unmdx   = np.copy(pos_x)

        #pos_x [ unmdx > bx[0]* 0.5   ] = np.subtract ( unmdx [ unmdx > bx[0]* 0.5  ], bx[0] )
        #pos_x [ unmdx < bx[0]*(-0.5) ] = np.add      ( unmdx [ unmdx < bx[0]*(-0.5)], bx[0] )
        
        pos_y   = all_atoms.positions[:,1]
        unmdy   = np.copy(pos_y)

        #pos_y [ unmdy > bx[1]*0.5    ] = np.subtract ( unmdy [ unmdy > bx[1]*0.5   ], bx[1] )
        #pos_y [ unmdy < bx[1]*(-0.5) ] = np.add      ( unmdy [ unmdy < bx[1]*(-0.5)], bx[1] )

	pos_z   = all_atoms.positions[:,2]
        unmdy   = np.copy(pos_z)
        
        kx = 0
        for iqx in qx:
            ky = 0
            for iqy in qy:
		kz=0
            	for iqz in qz:
    
                    fac  = np.add ( np.multiply (pos_x, iqx), np.multiply (pos_y, iqy),np.multiply (pos_z, iqz))
                    fac1 = ( np.sum( np.multiply ( np.cos(fac), atom_numbr)) ) + 1j * ( np.sum( np.multiply ( np.sin(fac), atom_numbr)))

                    I[kx,ky,kz] += fac1 #
                    kz +=1
                ky += 1
            kx += 1
    
        ntrj  += 1        
        I_save = np.add ( np.copy(I_save), np.square(np.absolute(I)) )

    index += 1

I_save = np.divide ( np.copy(I_save), (Ig*ntrj) )

nres_p_x  = int((math.floor((bx[0]*3.0)/(math.pi*2.0))))
nres_p_y  = int((math.floor((bx[1]*3.0)/(math.pi*2.0))))
nres_p_z  = int((math.floor((bx[2]*3.0)/(math.pi*2.0))))

#make_image ( I_save [nresx - nres_p_x: nresx + 1 + nres_p_x, nresy - nres_p_y: nresy + 1 + nres_p_y] , r'$q_x  (A^{\circ -1})$',r'$q_y (A^{\circ -1}$)', '2d_scatter.pdf')

rad_res = 0.05
nrad    = int( round ((3.0/rad_res)))

q_xyz     = np.linspace (0,3.0,num=nrad+1, endpoint=True)

I_xyz     = np.zeros_like (q_xyz)
n_xyz     = np.zeros_like (q_xyz)
d_q_xyz_c = q_xyz[1] - q_xyz[0]

kx = 0
for iqx in qx:
    ky = 0
    for iqy in qy:
	kz=0;
	for iqz in qz:
            qxyz_c = math.sqrt((iqx**2) + (iqy**2)+ (iqz**2))

            if ( qxyz_c <= 3.0):
                
                bin_qxyz = int(qxyz_c / d_q_xyz_c)
                I_xyz [ bin_qxyz ] += (I_save [kx , ky, kz])
                n_xyz [ bin_qxyz ] += 1.0 
            kz += 1
        ky += 1
    kx += 1

#print n_xy

I_xyz = np.divide ( np.copy (I_xyz) , n_xyz)
#make_plot  ( q_xy, I_xy , r'$q_{xy} A^{\circ -1}$',r'$I(q_{xy})$', 'xrd_2d.pdf')

np.savetxt ( 'I_qr.dat',  zip(q_xyz, I_xyz), fmt="%16.8f", delimiter=' ')
#np.savetxt ( 'I_q.dat', I_save,          fmt="%16.8f", delimiter=' ')

#########################################################################################################################################
