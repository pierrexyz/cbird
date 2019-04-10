#!/usr/bin/env python
from __future__ import print_function
import time
import numpy as np
# import scipy.stats
import os
import pandas as pd
import tools_7_param
import wrapper_7_param
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

t0 = time.time()


"""The goal of this script is to computed the direct power spectrum for 
a range of test cosmologies and the Taylor expansion at second order around 
a defined fiducial cosmology. It returns the direct power spectra, the taylor ones,
the ratio taylor/direct and the difference"""


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



#Fiducial cosmology 
theta_fid = np.load('precomputed/BOSS/theta_fid_3param.npy')
sigma = np.load('precomputed/BOSS/sigma_3param.npy')
fb = np.load('precomputed/BOSS/fb_3param.npy')



###########################################
#  Grid ##############################
###########################################
#

"""GRID TEST HYPERSPHERE"""

#for test purposes
#theta_fid = np.array([1., 1., 1., 1., 1., 1., 1.])
#sigma = np.array([1., 1., 1., 1., 1., 1., 1.])

N = 100     #number of desired test cosmologies
dim = 3     #here, we focus on 3 parameters
n_sigma = 3     #number of sigma for the radius of the hypersphere

norm = np.random.normal
normal_deviates = norm(size=(dim, N))

#build a sphere of radius one
radius = np.sqrt((normal_deviates**2).sum(axis=0))
points = np.transpose(normal_deviates/radius)

theta_tab_taylor = np.zeros((N+1,7))
index_param = [1, 2, 3]    #the code is built to investigate 7 parameters variations
#here, we just use Om, h and lnAs which are labelled 1, 2, 3

for i in range(len(theta_tab_taylor)):
    if i == 0 : 
        theta_tab_taylor[i] = theta_fid    #the first cosmo is always the fiducial
    else : 
        theta_tab_taylor[i] = theta_fid 
        for j in index_param :
            theta_tab_taylor[i, j] += points[i-1, j-1] * n_sigma * sigma[j]

np.save("input/theta_tab_taylor_hypersphere_BOSS_3param.npy", theta_tab_taylor)
#obtain N cosmo distributed on a sphere of radius 3 sigma


###########################################
#  Function_s  ###########################
###########################################


def get_EFT_pars(ob, Om, h, lnAs, n_s, Omega_k, w, config_file, nrank):
    omega_m = Om*h**2
    #omega_k = Omega_k*h**2
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
Ntot = len(theta_tab_taylor)
sizered = Ntot/ncores


def testcs(arraytheta, nrank):
    arrayred = arraytheta[nrank*sizered:(nrank+1)*sizered]
    allPlin = []
    allP1loop = []

    for i, theta in enumerate(arrayred):
        print(i)
        idx = nrank * sizered + i
        Plin, P1loop = CompPterms(theta, 'config_3_param_taylor_crash_BOSS.ini', nrank)

        idxcol = np.full([Plin.shape[0], 1], idx)
        
        allPlin.append(np.hstack([Plin, idxcol]))
        allP1loop.append(np.hstack([P1loop, idxcol]))
        
    np.save("results/BOSS/3sigma3param/Plin_direct_rank%s.npy" % str(nrank), np.array(allPlin))
    np.save("results/BOSS/3sigma3param/P1loop_direct_rank%s.npy" % str(nrank), np.array(allP1loop))
    return

testcs(theta_tab_taylor, rank)

print(time.time() - t0)



###########################################
#  Taylor expansion part ##################
###########################################


"""LOAD PRECOMPUTED PARAM"""

theta_fid = np.load("precomputed/BOSS/3param/theta_fid.npy")
sigma = np.load("precomputed/BOSS/3param/sigma.npy")

P1loop_l0_fid = np.load("precomputed/BOSS/3param/P1loop_l0_fid.npy")
P1loop_l2_fid = np.load("precomputed/BOSS/3param/P1loop_l2_fid.npy")
P1loop_l4_fid = np.load("precomputed/BOSS/3param/P1loop_l4_fid.npy")

fd_l0 = np.load("precomputed/BOSS/3param/fd_l0.npy")
fd_l2 = np.load("precomputed/BOSS/3param/fd_l2.npy")
fd_l4 = np.load("precomputed/BOSS/3param/fd_l4.npy")

hess_l0 = np.load("precomputed/BOSS/3param/hess_l0.npy")
hess_l2 = np.load("precomputed/BOSS/3param/hess_l2.npy")
hess_l4 = np.load("precomputed/BOSS/3param/hess_l4.npy")


"""Taylor EFT power spectrum"""

def TaylorEFT_3param(theta, theta_fid, P1loop_fid, fd, hess):
    dtheta = theta - theta_fid
    dtheta_interest = dtheta[1:4]
    P1loop = P1loop_fid + np.tensordot(dtheta_interest, fd, axes=(0,0)) + 0.5 * np.tensordot(dtheta_interest, np.tensordot(hess, dtheta_interest, axes=(0,0)), axes=(0,0))
    
    return P1loop


"""TEST"""

theta_tab = np.load('input/theta_tab_taylor_hypersphere_BOSS_3param.npy')

P1loop = np.load('results/BOSS/3sigma3param/P1loop_direct_rank0.npy')

P1loop_l0 = P1loop[:,:50,:]
P1loop_l2 = P1loop[:,50:100,:]
P1loop_l4 = P1loop[:,100:,:]

kvals = P1loop_l0[0,:,0]
P_direct_l0 = P1loop_l0[:,:,1:-1]
P_direct_l2 = P1loop_l2[:,:,1:-1]
P_direct_l4 = P1loop_l4[:,:,1:-1]

#compute Taylor :
P_taylor_l0 = np.zeros(np.shape(P_direct_l0))
P_taylor_l2 = np.zeros(np.shape(P_direct_l0))
P_taylor_l4 = np.zeros(np.shape(P_direct_l0))

for i in range(len(theta_tab)):
    P_taylor_l0[i] = TaylorEFT_3param(theta_tab[i],theta_fid, P1loop_l0_fid, fd_l0, hess_l0)
    P_taylor_l2[i] = TaylorEFT_3param(theta_tab[i],theta_fid, P1loop_l2_fid, fd_l2, hess_l2)
    P_taylor_l4[i] = TaylorEFT_3param(theta_tab[i],theta_fid, P1loop_l4_fid, fd_l4, hess_l4)
    
#comparison : 
ratio_l0 = P_taylor_l0/P_direct_l0
ratio_l2 = P_taylor_l2/P_direct_l2
ratio_l4 = P_taylor_l4/P_direct_l4

diff_l0 = P_taylor_l0 - P_direct_l0
diff_l2 = P_taylor_l2 - P_direct_l2
diff_l4 = P_taylor_l4 - P_direct_l4


np.save("results/BOSS/3sigma3param/ratio_l0.npy", ratio_l0)
np.save("results/BOSS/3sigma3param/ratio_l2.npy", ratio_l2)
np.save("results/BOSS/3sigma3param/ratio_l4.npy", ratio_l4)

np.save("results/BOSS/3sigma3param/diff_l0.npy", diff_l0)
np.save("results/BOSS/3sigma3param/diff_l2.npy", diff_l2)
np.save("results/BOSS/3sigma3param/diff_l4.npy", diff_l4)

np.save("results/BOSS/3sigma3param/P_taylor_l0.npy", P_taylor_l0)
np.save("results/BOSS/3sigma3param/P_taylor_l2.npy", P_taylor_l2)
np.save("results/BOSS/3sigma3param/P_taylor_l4.npy", P_taylor_l4)

np.save("results/BOSS/3sigma3param/P_direct_l0.npy", P_direct_l0)
np.save("results/BOSS/3sigma3param/P_direct_l2.npy", P_direct_l2)
np.save("results/BOSS/3sigma3param/P_direct_l4.npy", P_direct_l4)