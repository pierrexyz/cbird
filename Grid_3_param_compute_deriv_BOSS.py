#!/usr/bin/env python
from __future__ import print_function
import time
import numpy as np
import os
import pandas as pd
import tools_7_param
import wrapper_7_param
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

t0 = time.time()

"""The goal of this script is to build the two-point finite difference derivatives
with respect to three cosmological parameters : Omega_matter, h and ln(10^10 As).
It returns the power spectrum at fiducial cosmology, the first derivative and 
the hessian. To do so, it computes a grid of cosmologies that will be used for 
the finite differences."""

###########################################
#  Globals ##############################
###########################################

# PATH GLOBALS & tests to check they are as expected

# Save this path (expect data and other needed executables present in relative paths!!!)
THIS_PATH = os.path.dirname(__file__)

# Output path: it will be in scratch
OUTPATH = os.path.abspath(os.path.join(THIS_PATH,'output'))
if not os.path.isdir(OUTPATH):
    raise Exception(OUTPATH + ' not there!')

# The CLASS Boltzmann code in_stallation location
# Expects the CLASS code compiled here
CLASSPATH = os.path.abspath(os.path.join(THIS_PATH, 'class_public'))

# How to find the EFT model files and tools_7_param
# Expects the EFT code compiled here
EFT_PATH = os.path.abspath(os.path.join(THIS_PATH, 'cbird-precise'))

# COSMOLOGICAL GLOBALS: fiducial model (should match input sim data!)
dfcosmo = pd.read_csv(os.path.join(THIS_PATH, 'input',
                                   'BOSS_3param.csv'), index_col=0)

simtype = 'BOSS'
series_cosmo = dfcosmo.loc[simtype]
gridname = 'GridSDSS'

# Cosmological parameters : BOSS
    #Nb : o_x = O_x * h**2
ocfid = series_cosmo.loc['Omega_m'] * \
    series_cosmo.loc['h']**2 - series_cosmo.loc['omega_b']
obfid = series_cosmo.loc['omega_b']
fb = obfid/ocfid
Omega_mfid = series_cosmo['Omega_m']
h_fid = series_cosmo['h']
lnAs_fid = series_cosmo['lnAs']
ns_fid = series_cosmo['ns']
Omega_kfid = series_cosmo['Omega_k']
w_fid = series_cosmo['w']
z_pk = series_cosmo['z_pk']

#one-sigma departure from Planck 2018
sigma_ob = series_cosmo['sigma_ob']
sigma_Omega_m = series_cosmo['sigma_Omega_m']
sigma_h = series_cosmo['sigma_h']
sigma_lnAs = series_cosmo['sigma_lnAs']
sigma_ns = series_cosmo['sigma_ns']
sigma_Omega_k = series_cosmo['sigma_Omega_k']
sigma_w = series_cosmo['sigma_w']

fid = np.array([obfid, Omega_mfid, h_fid, lnAs_fid, ns_fid, Omega_kfid, w_fid])
sigma = np.array([sigma_ob, sigma_Omega_m, sigma_h, sigma_lnAs, sigma_ns, sigma_Omega_k, sigma_w])

np.save("precomputed/BOSS/3param/theta_fid.npy", fid)
np.save("precomputed/BOSS/3param/sigma.npy", sigma)
np.save("precomputed/BOSS/3param/fb.npy", fb)




###########################################
#  Grid ##############################
###########################################

"""the grid here aims just to construct finite difference for Om, As and h """

#number of cosmology we need to compute to get the hessian: 
#n_cosmo = 1 for fid + 3 * 2 for direct derivatives + 3 * 4 for mixed derivatives 
n_cosmo = 19

theta_tab_deriv = np.zeros((n_cosmo, len(fid)))

for i in range(n_cosmo): 
    theta_tab_deriv[i] = fid
    
#fiducial cosmo, then 1 param variations, then 2 param variation_s (++,+-,-+,--)
    #it is creating the grid with all the required cosmologies for doing 
    #the finite difference
cnt = 6      #filling count
index_param = [1, 2, 3]    #the code is made for 7 param : here just look at Om, h and lnAs which are respectively labelled 1, 2 and 3
for i in range(3):
    theta_tab_deriv[i*2 + 1, index_param[i]] += sigma[index_param[i]]
    theta_tab_deriv[i*2 + 2, index_param[i]] -= sigma[index_param[i]]
    for j in range(3):
        if j > i :
            theta_tab_deriv[cnt + 1, index_param[i]] += sigma[index_param[i]] #++
            theta_tab_deriv[cnt + 1, index_param[j]] += sigma[index_param[j]]
            theta_tab_deriv[cnt + 2, index_param[i]] += sigma[index_param[i]] #+-
            theta_tab_deriv[cnt + 2, index_param[j]] -= sigma[index_param[j]]
            theta_tab_deriv[cnt + 3, index_param[i]] -= sigma[index_param[i]] #-+
            theta_tab_deriv[cnt + 3, index_param[j]] += sigma[index_param[j]]
            theta_tab_deriv[cnt + 4, index_param[i]] -= sigma[index_param[i]] #--
            theta_tab_deriv[cnt + 4, index_param[j]] -= sigma[index_param[j]]
            cnt += 4    



###########################################
#  Function_s  ###########################
###########################################


def get_EFT_pars(ob, Om, h, lnAs, n_s, Omega_k, w, config_file, nrank):
    omega_m = Om*h**2
    pars = tools_7_param.get_config(bigconfig_file=config_file, cat=True)
    pars['omega_cdm'] = omega_m / (1. + fb)
    pars['omega_b'] = ob
    pars['omega_m'] = omega_m
    pars['h'] = h
    pars['ln10^{10}A_s'] = lnAs
    pars['n_s'] = n_s
    pars['Omega_k'] = Omega_k
    pars['w0_fld'] = w    
    return pars


def CompPterms(theta, config_file, nrank, run=0):
    """Run the power spectra specified parameters.
    Assume PS are given in a file concatenated by multipoles.
    Inputs:
    - theta: Set of parameters for a single power spectrum
    - nrank: Rank of the MPI process
    - run: Index of the run
    Outputs:
    Plin, Ploop"""
    ob, Om, h, lnAs, n_s, Omega_k, w = theta
    # Get the EFT parameters
    pars = get_EFT_pars(ob, Om, h, lnAs, n_s, Omega_k, w, config_file, nrank)
    # Run_s Pierre's code and save output to folder in path.
    path = wrapper_7_param.run_zbEFT(pars)
    # Get the k-values
    # For l=0 
    Plin = np.loadtxt(os.path.abspath(
        os.path.join(path, 'PowerSpectraLinear.dat')))
    P1loop = np.loadtxt(os.path.abspath(
        os.path.join(path, 'PowerSpectra1loop.dat')))
    
    return Plin, P1loop

###########################################
#  Data  ###########################
###########################################


ncores = size
Ntot = len(theta_tab_deriv)
sizered = Ntot/ncores


def testcs(arraytheta, nrank):
    arrayred = arraytheta[nrank*sizered:(nrank+1)*sizered]
    allPlin = []
    allP1loop = []

    for i, theta in enumerate(arrayred):
        print(i)
        idx = nrank * sizered + i
        Plin, P1loop = CompPterms(theta, 'config_3_param_deriv_BOSS.ini', nrank)

        idxcol = np.full([Plin.shape[0], 1], idx)
        
        allPlin.append(np.hstack([Plin, idxcol]))
        allP1loop.append(np.hstack([P1loop, idxcol]))
        
    np.save("output/BOSS/3param/Plin_rank%s.npy" % str(nrank), np.array(allPlin))
    np.save("output/BOSS/3param/P1loop_rank%s.npy" % str(nrank), np.array(allP1loop))
    return

testcs(theta_tab_deriv, rank)

print(time.time() - t0)



###########################################
#  Derivatives  ###########################
###########################################

#load output from c_bird
P1loop = np.load('output/BOSS/3param/P1loop_rank0.npy')

#separate monopole, quadrupole and hexadecapole
P1loop_l0 = P1loop[:,:50,:]
P1loop_l2 = P1loop[:,50:100,:]
P1loop_l4 = P1loop[:,100:,:]

#number of parameters
n_pars = 3

kvals = P1loop_l0[0,:,0]
P1loop_l0_fid = P1loop_l0[0,:,1:-1]
P1loop_l2_fid = P1loop_l2[0,:,1:-1]
P1loop_l4_fid = P1loop_l4[0,:,1:-1]

#fd = first derivative
#hess = hessian matrix

fd_l0 = np.zeros((n_pars, len(P1loop_l0[0,:,0]), len(P1loop_l0[0,0,1:-1])))
fd_l2 = np.zeros((n_pars, len(P1loop_l0[0,:,0]), len(P1loop_l0[0,0,1:-1])))
fd_l4 = np.zeros((n_pars, len(P1loop_l0[0,:,0]), len(P1loop_l0[0,0,1:-1])))

hess_l0 = np.zeros((n_pars, n_pars, len(P1loop_l0[0,:,0]), len(P1loop_l0[0,0,1:-1])))
hess_l2 = np.zeros((n_pars, n_pars, len(P1loop_l0[0,:,0]), len(P1loop_l0[0,0,1:-1])))
hess_l4 = np.zeros((n_pars, n_pars, len(P1loop_l0[0,:,0]), len(P1loop_l0[0,0,1:-1])))

#the counting to fill the hessian is stupid but I didn't find better. Sorry if
#it's unclear. The idea is just to use the computed cosmologies to get the first 
#derivative and the hessian by first order finite difference. 
cnt = 0
for i in range(len(P1loop_l0)):  
    if i < 2*n_pars : 
        if i % 2 == 1 :
            print('i = ', i )
            print('int(i/2) = ', int(i/2) )
            fd_l0[int(i/2)] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1])/(2*sigma[int(i/2)+1])#int(1/2) = 0, 1, 2 and sigma_interest = 1, 2, 3 
            fd_l2[int(i/2)] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1])/(2*sigma[int(i/2)+1])
            fd_l4[int(i/2)] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1])/(2*sigma[int(i/2)+1])
            
            hess_l0[int(i/2), int(i/2)] = (P1loop_l0[i,:,1:-1] - 2*P1loop_l0_fid + P1loop_l0[i+1,:,1:-1])/(sigma[int(i/2)+1]**2)
            hess_l2[int(i/2), int(i/2)] = (P1loop_l2[i,:,1:-1] - 2*P1loop_l2_fid + P1loop_l2[i+1,:,1:-1])/(sigma[int(i/2)+1]**2)
            hess_l4[int(i/2), int(i/2)] = (P1loop_l4[i,:,1:-1] - 2*P1loop_l4_fid + P1loop_l4[i+1,:,1:-1])/(sigma[int(i/2)+1]**2)

    else : 
        if i == 7 :
            hess_l0[0, 1] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1] - P1loop_l0[i+2,:,1:-1] + P1loop_l0[i+3,:,1:-1])/(4*sigma[1]*sigma[2])
            hess_l2[0, 1] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1] - P1loop_l2[i+2,:,1:-1] + P1loop_l2[i+3,:,1:-1])/(4*sigma[1]*sigma[2])                
            hess_l4[0, 1] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1] - P1loop_l4[i+2,:,1:-1] + P1loop_l4[i+3,:,1:-1])/(4*sigma[1]*sigma[2])
                            
            hess_l0[1, 0] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1] - P1loop_l0[i+2,:,1:-1] + P1loop_l0[i+3,:,1:-1])/(4*sigma[1]*sigma[2])
            hess_l2[1, 0] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1] - P1loop_l2[i+2,:,1:-1] + P1loop_l2[i+3,:,1:-1])/(4*sigma[1]*sigma[2])                
            hess_l4[1, 0] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1] - P1loop_l4[i+2,:,1:-1] + P1loop_l4[i+3,:,1:-1])/(4*sigma[1]*sigma[2])
  

        if i == 11 :
            hess_l0[0, 2] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1] - P1loop_l0[i+2,:,1:-1] + P1loop_l0[i+3,:,1:-1])/(4*sigma[1]*sigma[3])
            hess_l2[0, 2] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1] - P1loop_l2[i+2,:,1:-1] + P1loop_l2[i+3,:,1:-1])/(4*sigma[1]*sigma[3])                
            hess_l4[0, 2] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1] - P1loop_l4[i+2,:,1:-1] + P1loop_l4[i+3,:,1:-1])/(4*sigma[1]*sigma[3])
                            
            hess_l0[2, 0] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1] - P1loop_l0[i+2,:,1:-1] + P1loop_l0[i+3,:,1:-1])/(4*sigma[1]*sigma[3])
            hess_l2[2, 0] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1] - P1loop_l2[i+2,:,1:-1] + P1loop_l2[i+3,:,1:-1])/(4*sigma[1]*sigma[3])                
            hess_l4[2, 0] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1] - P1loop_l4[i+2,:,1:-1] + P1loop_l4[i+3,:,1:-1])/(4*sigma[1]*sigma[3])
     

        if i == 15 :
            hess_l0[1, 2] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1] - P1loop_l0[i+2,:,1:-1] + P1loop_l0[i+3,:,1:-1])/(4*sigma[3]*sigma[2])
            hess_l2[1, 2] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1] - P1loop_l2[i+2,:,1:-1] + P1loop_l2[i+3,:,1:-1])/(4*sigma[3]*sigma[2])                
            hess_l4[1, 2] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1] - P1loop_l4[i+2,:,1:-1] + P1loop_l4[i+3,:,1:-1])/(4*sigma[3]*sigma[2])
                            
            hess_l0[1, 2] = (P1loop_l0[i,:,1:-1] - P1loop_l0[i+1,:,1:-1] - P1loop_l0[i+2,:,1:-1] + P1loop_l0[i+3,:,1:-1])/(4*sigma[3]*sigma[2])
            hess_l2[1, 2] = (P1loop_l2[i,:,1:-1] - P1loop_l2[i+1,:,1:-1] - P1loop_l2[i+2,:,1:-1] + P1loop_l2[i+3,:,1:-1])/(4*sigma[3]*sigma[2])                
            hess_l4[1, 2] = (P1loop_l4[i,:,1:-1] - P1loop_l4[i+1,:,1:-1] - P1loop_l4[i+2,:,1:-1] + P1loop_l4[i+3,:,1:-1])/(4*sigma[3]*sigma[2])
                                           

#save derivatives
np.save("precomputed/BOSS/3param/kvals.npy", kvals)
np.save("precomputed/BOSS/3param/P1loop_l0_fid.npy", P1loop_l0_fid)
np.save("precomputed/BOSS/3param/P1loop_l2_fid.npy", P1loop_l2_fid)
np.save("precomputed/BOSS/3param/P1loop_l4_fid.npy", P1loop_l4_fid)
np.save("precomputed/BOSS/3param/fd_l0.npy", fd_l0)
np.save("precomputed/BOSS/3param/fd_l2.npy", fd_l2)
np.save("precomputed/BOSS/3param/fd_l4.npy", fd_l4)
np.save("precomputed/BOSS/3param/hess_l0.npy", hess_l0)
np.save("precomputed/BOSS/3param/hess_l2.npy", hess_l2)
np.save("precomputed/BOSS/3param/hess_l4.npy", hess_l4)


