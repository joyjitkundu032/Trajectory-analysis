import os
import sys
import MDAnalysis as md
import my_plotter as myplot
import backbone_dihedrals as bb_dihedrals
import math
import numpy as np
import scipy as sp
import numpy.linalg
import itertools

from scipy import signal

#md.core.flags['use_pbc'] = True

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx], idx

for file in os.listdir("../input"):
    if file.endswith(".psf"):
        PSF_file = '../input/'+file

traj_list = []
for file in os.listdir("../data"):
    if file.endswith(".dcd"):
        traj_list.append('../data/'+file)

print traj_list.sort()
print PSF_file, traj_list

#traj_list = ['../data/spep_r3.dcd']

nresidue     = int(sys.argv[1])
u            = md.Universe(PSF_file, traj_list)

traj_length  = len(u.trajectory)
prin_moments = np.zeros((3,traj_length))
time         = np.zeros(traj_length)

psi_traj_U   = np.zeros((nresidue-1, traj_length)); phi_traj_U   = np.zeros((nresidue-1, traj_length)); ome_traj_U   = np.zeros((nresidue-1, traj_length)); 
phi_U        = np.zeros(0); psi_U        = np.zeros(0); ome_U        = np.zeros(0); 

calpha_dist  = np.zeros((nresidue-1, traj_length))
rho_bins     = 500
water_dens   = np.zeros(rho_bins) ; anion_dens   = np.zeros(rho_bins); ction_dens   = np.zeros(rho_bins)

water_O      = u.select_atoms ( "(segid WT1 or segid WT2 or segid WT3 or segid WT4) and name OH2" ) ;
anion_O      = u.select_atoms ( "segid U and name OE2" ) ;
ction_O      = u.select_atoms ( "segid U and name NZ" ) ;

peptoid      = u.select_atoms ( "segid U" ) ;

bb_Nitro_U   = u.select_atoms ( "segid U and (name N or name NR)" )  ;
bb_Calph_U   = u.select_atoms ( "segid U and name CA" ) ;
bb_Ccarb_U   = u.select_atoms ( "segid U and name C" )  ; 

nterm_U      = u.select_atoms ( "atom U 1 N" ) ;   
cterm_U      = u.select_atoms ( "atom U "+str(nresidue)+" NR" ) ; 
bb_atoms_U   = u.select_atoms ( "segid U and (name N or name CA or name C)" ) ;

COM_ring     = np.zeros((nresidue/2, traj_length, 3)); 
COM_anin     = np.zeros((nresidue/4, traj_length, 3))
COM_ctin     = np.zeros((nresidue/4, traj_length, 3))

solv_t       = [] ; n_water_solv = [] ; leg_res      = [] ; index        = 0 ; nequib       = 0 ; interface_loc_z = [];

r2g_U        = np.array([(bb_atoms_U.radius_of_gyration()) for ts in u.trajectory])
end_to_end_U = np.array([(np.linalg.norm(cterm_U.positions - nterm_U.positions)) for ts in u.trajectory])
COM_pep      = np.array([(peptoid.center_of_mass()) for ts in u.trajectory])

for ires in range(1,nresidue/2,2):
    ctin                = u.select_atoms("resid " + str(ires) + " and name NZ")
    COM_ctin[index,:,:] = np.array([(ctin.positions[0,:]) for ts in u.trajectory])
    index              += 1

index = 0
for ires in range(nresidue/2,nresidue,2):
    anin                = u.select_atoms("resid " + str(ires+1) + " and (name OE2 or name OE1)")
    COM_anin[index,:,:] = np.array([(anin.center_of_mass()) for ts in u.trajectory])
    index              += 1

for ires in range(1,nresidue,2):
    COM_phenyl               = u.select_atoms("resid " + str(ires+1) + " and (name CG or name CD1 or name CD2 or name CE1 or name CE2 or name CZ)") ;
    COM_ring[(ires-1)/2,:,:] = np.array([(COM_phenyl.center_of_mass()) for ts in u.trajectory])

########################################################################################################################################################

intr_index_0 = rho_bins / 4; intr_index_f = ((3 * rho_bins) / 4)

anin_pos     = np.zeros((nresidue/4, 3))

index = 0
for ts in u.trajectory:     # iterate through all frames

    curr_time = ts.time
    print 'Time = ', curr_time
    
    if ((index+1)%40 == 0):

        water_solv   = u.select_atoms("resname TIP3 and (around 6.0 (segid U))")
        n_water_solv = np.append( n_water_solv, water_solv.n_atoms)
        solv_t       = np.append( solv_t, curr_time)

    A,B,C,alpha, beta, gamma = ts.dimensions

    bx                  = np.array((A, B, C))
    time      [index]   = curr_time
    #prin_moments[:,index], principal_axes  = np.linalg.eig( peptoid.moment_of_inertia() )
    prin_moments[:,index], principal_axes  = np.linalg.eig( bb_atoms_U.moment_of_inertia() )
    
    phi_curr_traj_U, psi_curr_traj_U, ome_curr_traj_U = bb_dihedrals.backbone_dihedral_angles (nresidue, bb_Nitro_U.positions, \
                                                                                               bb_Calph_U.positions, bb_Ccarb_U.positions)
    phi_traj_U[:,index] = phi_curr_traj_U; psi_traj_U[:, index] = psi_curr_traj_U; ome_traj_U[:, index] = ome_curr_traj_U;    

    calpha_dist_curr = bb_dihedrals.calpha_distances(nresidue, bb_Nitro_U.positions, bb_Calph_U.positions, bb_Ccarb_U.positions)
    calpha_dist[:, index] = calpha_dist_curr

    index            += 1

    water_histo, w_edges = np.histogram (water_O.positions[:,2], bins = rho_bins, range=(-C/2.0, C/2.0))
    anion_histo, a_edges = np.histogram (COM_anin [:,index-1,2], bins = rho_bins, range=(-C/2.0, C/2.0))
    ction_histo, c_edges = np.histogram (ction_O.positions[:,2], bins = rho_bins, range=(-C/2.0, C/2.0))
    
    bin_width            = w_edges[1] - w_edges[0]
    water_dens_curr      = np.divide ( water_histo, (A * B * bin_width) )

    #rho_inter, indx      = find_nearest( water_dens_curr[rho_bins/4:(3*rho_bins)/4], np.amax (water_dens_curr [rho_bins/4:(3*rho_bins)/4])/2.0) 

    rho_inter, indx      = find_nearest ( water_dens_curr[intr_index_0:intr_index_f], 0.0165 ) 
    interface_loc_z      = np.append ( interface_loc_z, w_edges[intr_index_0 + indx] + (0.50 * bin_width))

    indx        = intr_index_0 + indx;  intr_index_0 = indx - (rho_bins/4);    intr_index_f = indx + (rho_bins/4)
    print intr_index_0, intr_index_f
    if (intr_index_f > rho_bins):
        intr_index_f = rho_bins
        intr_index_0 = rho_bins/2
        
    if (intr_index_0 < 0):
        intr_index_0 = 0
        intr_index_f = rho_bins/2

    if (curr_time >= 100000):

        nequib              += 1

        water_dens           = np.add(water_dens, water_dens_curr)
        anion_dens           = np.add(anion_dens, np.divide(anion_histo, (A * B * (a_edges[1] - a_edges[0]))))
        ction_dens           = np.add(ction_dens, np.divide(ction_histo, (A * B * (c_edges[1] - c_edges[0]))))

        phi_U = np.append(phi_U, phi_curr_traj_U); psi_U = np.append(psi_U, psi_curr_traj_U); ome_U = np.append(ome_U, ome_curr_traj_U);
        
########################################################################################################################################################

water_dens = np.divide (water_dens, nequib); anion_dens = np.divide (anion_dens, nequib); ction_dens = np.divide (ction_dens, nequib)
time       = np.multiply(time , 0.001  ) ;   solv_t     = np.multiply(solv_t, 0.001 ) ;

fe_file = open('r2g.dat','w')
fe_file.write("%16s %16s\n" %( "time", "R_g"));
for a,b in itertools.izip(time, r2g_U):
    fe_file.write("%16.8f %16.8f\n" %( a, b));
fe_file.close();

fe_file = open('e2e.dat','w')
fe_file.write("%16s %16s\n" %( "time", "R_E"));
for a,b in itertools.izip(time, end_to_end_U):
    fe_file.write("%16.8f %16.8f\n" %( a, b));
fe_file.close();

fe_file = open('nw.dat','w')
fe_file.write("%16s %16s\n" %( "time", "N_w"));
for a,b in itertools.izip(solv_t,n_water_solv):
    fe_file.write("%16.8f %16.8f\n" %( a, b));
fe_file.close();

fe_file = open('dens_prof.dat','w')
fe_file.write("%16s %16s %16s %16s %16s %16s\n" %( "z_w", "rho_w", "z_a", "rho_a", "z_c", "rho_c" ));
for a,b,c,d,e,f in itertools.izip( w_edges[:-1], water_dens, a_edges[:-1], anion_dens, c_edges[:-1], ction_dens):
    fe_file.write("%16.8f %16.8f %16.8f %16.8f %16.8f %16.8f\n" %( a, b, c, d, e, f));
fe_file.close();

Pep_phenyl_z = np.zeros((nresidue/2, traj_length)); Pep_anion_z = np.zeros((nresidue/4, traj_length)); Pep_cation_z = np.zeros((nresidue/4, traj_length));

index = 0
for ires in range(1,nresidue,2):
    Pep_phenyl_z [index, :] = np.subtract ( COM_ring [ index , : , 2], interface_loc_z ) 
    index += 1

index = 0
for ires in range(nresidue/2,nresidue,2):
    Pep_anion_z [index, :] = np.subtract ( COM_anin[index,:,2], interface_loc_z )
    index += 1

index = 0
for ires in range(1, nresidue/2, 2):
    Pep_cation_z [index, :] = np.subtract ( COM_ctin[index,:,2], interface_loc_z )
    index += 1

Pep_z = np.subtract ( COM_pep[:,2], interface_loc_z )

myplot.make_plot  ( time,   end_to_end_U,      r'$t (ns)$', r'$R_{EE}$', 'e2e.pdf'   ) 
myplot.make_plot  ( time,   interface_loc_z,   r'$t (ns)$', r'$z_{interface}$', 'interface_z.pdf'   ) 
myplot.make_plot  ( time,   r2g_U,             r'$t (ns)$', r'$R_g$',    'r2g.pdf'   )
myplot.make_rama  ( phi_U,  psi_U,             r'$\phi$',   r'$\psi$',   'ramachandran_U.pdf', nresidue )

print np.shape(phi_U)

#myplot.make_rama  ( phi_U[:,0:4],  psi_U[:,0:4],   r'$\phi$',   r'$\psi$',   'ramachandran_ex1.pdf' )
#myplot.make_rama  ( phi_U[:,4:15], psi_U[:,4:15],  r'$\phi$',   r'$\psi$',   'ramachandran_int.pdf' )
#myplot.make_rama  ( phi_U[:,15:19],psi_U[:,15:19], r'$\phi$',   r'$\psi$',   'ramachandran_ex2.pdf' )

myplot.make_plot  ( time, Pep_phenyl_z[(nresidue/2)-1,:], r'$t (ns)$', r'$z_{PT}$', 'COM_terminal_z.pdf' )

myplot.make_plot2 ( time, Pep_z, np.zeros(traj_length),          r'$t (ns)$', r'$r_z$',    'COM_z.pdf' )
myplot.make_plotn ( time, Pep_phenyl_z, nresidue/2, r'$t (ns)$', r'$z_P$', 'COM_phenyl_z.pdf', hbar = False )
myplot.make_plotn ( time, Pep_anion_z,  nresidue/4, r'$t (ns)$', r'$z_A$', 'COM_anion_z.pdf',  hbar = False )
myplot.make_plotn ( time, Pep_cation_z, nresidue/4, r'$t (ns)$', r'$z_C$', 'COM_cation_z.pdf', hbar = False )

myplot.make_plotn ( time, prin_moments, 3, r'$t (ns)$',r'$\mathbf{I}$','moments.pdf' , leg=['x','y','z'] )

myplot.make_plot  ( solv_t, n_water_solv, r'$t (ns)$', r'$N_w$', 'water_number.pdf' )
myplot.make_plot3 ( w_edges[:-1], water_dens, a_edges[:-1], anion_dens, c_edges[:-1], ction_dens,r'$z$', r'$\rho(z)$','dens_prof_ions.pdf', xlim=(-C/2, 5.0) , leg=['Water','Anion','Cation'])

myplot.make_image ( psi_traj_U, r'$\psi_U$',   'Frame (Time)', 'Residue', 'psi_traj_U.pdf' )
myplot.make_image ( phi_traj_U, r'$\phi_U$',   'Frame (Time)', 'Residue', 'phi_traj_U.pdf' )
myplot.make_image ( ome_traj_U, r'$\omega_U$', 'Frame (Time)', 'Residue', 'ome_traj_U.pdf' )
myplot.make_histo ( calpha_dist, r'$r$', r'$C_{\alpha} - C_{\alpha}$', 'calpha_hist.pdf')

#tyme, nonbonded, total = np.loadtxt('pep_ener.dat', usecols=(0,9,10), skiprows=1,unpack=True); tyme = np.multiply(tyme , 0.01)
#myplot.make_plots ( time, nonbonded, total, r'$t (ns)$',r'$E$','energy.pdf' ) 

