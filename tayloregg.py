#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import os
import os.path as opa
from scipy.integrate import quad
from scipy.interpolate import interp1d

THIS_PATH = opa.dirname(opa.abspath(__file__))
OUTPATH = opa.abspath(opa.join(THIS_PATH, 'taylornest/'))
if not opa.isdir(OUTPATH): raise Exception(OUTPATH + ' not there!')
if not opa.isdir(opa.join(OUTPATH, 'intermediary/')): os.makedirs(opa.join(OUTPATH, 'intermediary/'))
INPATH = '/exports/pierre/EFTofBOSS/input/'
if not opa.isdir(INPATH): INPATH = opa.abspath(opa.join(THIS_PATH, 'input/'))

CLASSPATH = opa.abspath(opa.join(THIS_PATH, 'class/')) 
EFT_PATH = opa.abspath(THIS_PATH) 

if opa.isdir(EFT_PATH):
    import sys
    if EFT_PATH not in sys.path: 
        sys.path.append(EFT_PATH)
    import zbeft
else:
    raise Exception('Module not found at ' + EFT_PATH)

import APpowerspectraNkmu as AP

################################################################################################################

def Hubble(Om, z):
    return ((Om)*(1+z)**3.+(1-Om))**0.5

def DA(Om, z):
    r = quad(lambda x:1./Hubble(Om,x), 0, z)[0]
    return r/(1+z)

def get_AP_param(z_pk, Om, Om_fid):
    qperp  =  DA(Om, z_pk) / DA(Om_fid, z_pk)
    qpar  =  Hubble(Om_fid, z_pk) / Hubble(Om, z_pk)
    return qperp, qpar

def import_cosmoref_from_DataFrameCosmosim(simtype):
    import pandas as pd
    dfcosmo = pd.read_csv(opa.join(INPATH, 'DataFrameCosmosims.csv'),index_col=0)
    series_cosmo = dfcosmo.loc[simtype]
    omega_bfid = series_cosmo.loc['omega_b']
    omega_cfid =  series_cosmo.loc['Omega_m']*series_cosmo.loc['h']**2-series_cosmo.loc['omega_b']
    fb = omega_bfid/omega_cfid
    Omega_mfid = dfcosmo.loc[simtype,'Omega_m']
    hfid = dfcosmo.loc[simtype,'h']
    lnAsfid = dfcosmo.loc[simtype,'lnAs']
    z_pk = dfcosmo.loc[simtype,'z_pk']
    nsfid = dfcosmo.loc[simtype,'ns']

    # Reference cosmology and error-step differences
    cosmoref_full = np.array( [lnAsfid, Omega_mfid, hfid, omega_bfid, nsfid] )
    #cosmodelta_full = np.array([0.014, 0.0073, 0.0054, 0.000015, 0.0042])     # Planck2018
    #cosmodelta_full = np.array([0.03, 0.003, 0.006, 0.00002, 0.005]) # ~ 1/3 stat. errors from full challenge
    cosmodelta_full = np.array([0.015, 0.007, 0.006, 0.00002, 0.005]) # Planckish

    cosmoref_shift = np.array( [3.16, 0.295, 0.648, 0.295*0.648**2 * fb/(1.+fb), nsfid] )
    cosmoref = cosmoref_shift#cosmoref_full
    cosmodelta = cosmodelta_full
    #print (fb)

    return z_pk, cosmoref, cosmodelta

def get_EFT_pars(z_pk, theta, nrank = 0):
    fb = 0.1861897575130901 # Omega_b / Omega_c ChallengeA
    lnAs, Om, h, omega_b, ns = theta
    omega_m = Om*h**2
    omega_b = omega_m*fb/(1.+fb)
    Outpath=opa.abspath(opa.join(OUTPATH, 'intermediary/'))
    keys = ['outpath', 'zbEFT_path', 'CLASS_path', 'ln10^{10}A_s', 'n_s', 'h', 'omega_b', 'omega_cdm', 'z_pk', 'pid']
    values = [Outpath, EFT_PATH, CLASSPATH, lnAs, ns, h, omega_b, omega_m-omega_b, z_pk, '%s'%nrank]
    return dict(zip(keys, values))

def CompPterms(z_pk, theta, kout, Om_fid = 0, TableNkmu = None, nrank = 0):
    
    pars = get_EFT_pars(z_pk, theta, nrank)
    #path = "/home/pchonje/Documents/EFTofLSS/cbird/cbird-dev/taylorbird-nest/intermediary/zbEFT_103089189_time.st_956fd03-c5d9b6d-5d0281e_0/"
    path = zbeft.run_zbEFT(pars)
    kPlin = np.loadtxt(opa.abspath(opa.join(path, 'PowerSpectraLinear.dat')))
    Ploop = np.loadtxt(opa.abspath(opa.join(path, 'PowerSpectra1loop.dat')))[:,1:]
    Plin = kPlin[:,1:]
    k = kPlin[:,0]

    # remove spike-bugs from cbird
    bugs = [0.15, 0.3]
    for bug in bugs:
        bugplace = np.where(k==bug)[0]
        #print ("the bug %s is at %s"%(bug, bugplace))
        Plin = np.delete(Plin, bugplace, 0)
        Ploop = np.delete(Ploop, bugplace, 0)
        k = np.delete(k, bugplace, 0)

    k = k[:len(k)/3]

    Plin = np.swapaxes(Plin.reshape(3, len(k), Plin.shape[1] ), axis1 = 1, axis2 = 2)
    Ploop = np.swapaxes(Ploop.reshape(3, len(k), Ploop.shape[1] ), axis1 = 1, axis2 = 2)

    if Om_fid == 0:
        kred = k
    else: # Apply AP effect
        kred = k[k <= kout.max() + 0.1]
        qperp, qpar = get_AP_param(z_pk, theta[1], Om_fid)
        if TableNkmu is None: # No wedges
            Plin = AP.changetoAPnobinning(Plin, k, kred, qperp, qpar, nbinsmu=100)
            Ploop = AP.changetoAPnobinning(Ploop, k, kred, qperp, qpar, nbinsmu=100)
        else:
            Plin = AP.changetoAPbinning(Plin, k, kred, qperp, qpar, TableNkmu)
            Ploop = AP.changetoAPbinning(Ploop, k, kred, qperp, qpar, TableNkmu)
    
    Plin = interp1d(kred, Plin)(kout)
    Ploop = interp1d(kred, Ploop)(kout)
    return Plin, Ploop

def testcs(z_pk, arraytheta, kout, Om_fid = 0, TableNkmu = None, nrank = 0):
    allPlin = []
    allPloop = []
    for i, theta in enumerate(arraytheta):
        Plin, Ploop = CompPterms(z_pk, theta, kout, Om_fid, TableNkmu, nrank)
        allPlin.append(Plin)
        allPloop.append(Ploop)
    return np.array(allPlin), np.array(allPloop)

# Produce power spectrum at cosmologies on a compact stencil for numerical evaluation of first and second derivatives around reference cosmology
def egging(nest, z_pk, cosmoref, cosmodelta, kout, Om_fid = 0, TableNkmu = None, ncosmo = 0, nrank = 0):
    thetatab = []
    thetaname = []
    sign = [-1, 1]

    if ncosmo == 0:
        n = len(cosmoref)
    else:
        n = ncosmo

    if (rank == 0):
        print ("Taylorbird with %s cosmological parameters"%n)
    
    # compact stencil of cosmologies
    for i, c in enumerate(cosmoref[:n]):
        for s in sign:
            theta = np.copy(cosmoref)
            theta[i] += s * cosmodelta[i]
            thetatab.append(theta)
            thetaname.append('%s'%(s*(i+1)))
            for i2, c2 in enumerate([x for j, x in enumerate(cosmoref[:n]) if j > i]):
                for s2 in sign:
                    theta2 = np.copy(theta)
                    theta2[i2+1+i] += s2 * cosmodelta[i2+1+i]
                    thetatab.append(theta2)
                    thetaname.append('%s%s'%(s*(i+1),s2*(i2+1+i+1)))

    if (rank == 0):
        print ("Evaluating %s cosmologies..."%len(thetatab))
    # Evaluation of the power spectrum at all needed nearby cosmologies
    sizered=len(thetatab)/size
    arrayred = thetatab[nrank*sizered:(nrank+1)*sizered]

    Plin_send, Ploop_send = testcs(z_pk, arrayred, kout, Om_fid, TableNkmu, nrank)

    Plintab = None
    Plooptab = None

    if nrank == 0:
        # Number of: nearby cosmologies to evaluate, power spectra, multipoles, ks
        Plintab = np.empty((size * Plin_send.shape[0], Plin_send.shape[1], Plin_send.shape[2], Plin_send.shape[3]), dtype = Plin_send.dtype)
        Plooptab = np.empty((size * Ploop_send.shape[0], Ploop_send.shape[1], Ploop_send.shape[2], Ploop_send.shape[3]), dtype = Plin_send.dtype)
        if sizered * size != len(thetatab):
            arrayextra = thetatab[size*sizered:]
            Plin_extra, Ploop_extra = testcs(z_pk, arrayextra, kout, Om_fid, TableNkmu, size)
        # Evaluation of the power spectrum at reference cosmology
        print ("Evaluating reference cosmology...")
        Plin_ref, Ploop_ref = CompPterms(z_pk, cosmoref, kout, Om_fid, TableNkmu)
        np.save(opa.join(nest, "Plin_ref.npy"), Plin_ref)
        np.save(opa.join(nest, "Ploop_ref.npy"), Ploop_ref)

    comm.Gather(Plin_send, Plintab, root = 0)
    comm.Gather(Ploop_send, Plooptab, root = 0)

    if nrank == 0:
        if sizered * size != len(thetatab):
            Plintab = np.append(Plintab, Plin_extra, axis = 0)
            Plooptab = np.append(Plooptab, Ploop_extra, axis = 0)

        # first and second derivatives around reference cosmology
        fdlin = np.zeros((len(cosmoref), Plin_send.shape[1], Plin_send.shape[2], Plin_send.shape[3]))
        sdlin = np.zeros((len(cosmoref), len(cosmoref), Plin_send.shape[1], Plin_send.shape[2], Plin_send.shape[3]))
        fdloop = np.zeros((len(cosmoref), Ploop_send.shape[1], Ploop_send.shape[2], Ploop_send.shape[3]))
        sdloop = np.zeros((len(cosmoref), len(cosmoref), Ploop_send.shape[1], Ploop_send.shape[2], Ploop_send.shape[3]))
        #fdname = []
        #sdname = []

        for i, c in enumerate(cosmoref[:n]):
            fdlin[i] = 0.5 * ( Plintab[thetaname.index('%s'%(i+1))] - Plintab[thetaname.index('%s'%(-i-1))] ) / cosmodelta[i] 
            fdloop[i] = 0.5 * ( Plooptab[thetaname.index('%s'%(i+1))] - Plooptab[thetaname.index('%s'%(-i-1))] ) / cosmodelta[i]
            #fdname.append('%s'%(i+1))
            sdlin[i,i] = ( Plintab[thetaname.index('%s'%(i+1))] -2*Plin_ref + Plintab[thetaname.index('%s'%(-i-1))] ) / cosmodelta[i]**2
            sdloop[i,i] = ( Plooptab[thetaname.index('%s'%(i+1))] -2*Ploop_ref + Plooptab[thetaname.index('%s'%(-i-1))] ) / cosmodelta[i]**2
            #sdname.append('%s%s'%(i+1,i+1))
            for i2, c2 in enumerate([x for j, x in enumerate(cosmoref[:n]) if j > i]):
                k = i+1
                l = i2+1+i+1
                sdlin[k-1,l-1] = ( Plintab[thetaname.index('%s%s'%(k,l))] - Plintab[thetaname.index('%s%s'%(k,-l))] - Plintab[thetaname.index('%s%s'%(-k,l))] + Plintab[thetaname.index('%s%s'%(-k,-l))] ) / (4.*cosmodelta[k-1]*cosmodelta[l-1]) 
                sdloop[k-1,l-1] = ( Plooptab[thetaname.index('%s%s'%(k,l))] - Plooptab[thetaname.index('%s%s'%(k,-l))] - Plooptab[thetaname.index('%s%s'%(-k,l))] + Plooptab[thetaname.index('%s%s'%(-k,-l))] ) / (4.*cosmodelta[k-1]*cosmodelta[l-1])
                sdlin[l-1,k-1] = sdlin[k-1,l-1]
                sdloop[l-1,k-1] = sdloop[k-1,l-1]
                #sdname.append('%s%s'%(k,l))

        np.save(opa.join(nest,"fdlin.npy"), fdlin)
        np.save(opa.join(nest,"sdlin.npy"), sdlin)
        np.save(opa.join(nest,"fdloop.npy"), fdloop)
        np.save(opa.join(nest,"sdloop.npy"), sdloop)
    return

########################################################################################################################################################

if __name__ ==  "__main__":

    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()

    if (rank == 0):
        print('size:', size, ', rank:', rank)

    boxnumber = sys.argv[1]
    kmin = 0.01
    kmax = float(sys.argv[2])
    simtype = sys.argv[3]
    ZONE = 'NGC'
    if "Challenge" in simtype:
        ZONE = ''

    ##################################################
    ##### Loading covariance and Nkmu binning ###
    if 'Challenge' in simtype:
        simname = 'Challenge%s'%boxnumber
        TableNkmu = np.loadtxt(opa.join(INPATH,'Binning/Nkmu%s%s.dat' % (simtype, boxnumber))).T
        if 'Quarter' not in simtype:
            simtype_false = 'ChallengeQuarter%s'%boxnumber
            #print('Using quarter covariance divided by 4.25 instead of full')
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype_false,ZONE))) / 4.25
        else:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype,ZONE)))
    else:
        simname = simtype
        TableNkmu = None
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%sdata.dat'%(simtype,ZONE)))

    ##################################################
    ##### Loading power spectrum data ###
    kPS, PSdata, _ = np.loadtxt(opa.join(INPATH, 'DataSims/ps1D_%s%s_%s.dat' % (simtype, ZONE, boxnumber)), unpack=True)

    indexkred = np.argwhere((kPS < kmax) & (kPS > kmin))[:, 0]
    xdata = kPS[indexkred]
    ydata = PSdata[indexkred]

    kmask = np.array([False] * len(kPS))
    kmask[indexkred] = True
    kmask = kmask[:int(len(kmask) / 3)]

    Covred = Full_Cov[indexkred.reshape((len(indexkred), 1)), indexkred]
    Cinv = np.linalg.inv(Covred)

    masktriangle = []

    if 'Challenge' in simtype:  # periodic box, no window function
        kpred = xdata[:len(xdata)/3]
        Cinvw = Cinv
        Cinvww = Cinv
    else:                       # apply window function
        import WindowFunctionFourier
        kpred, Cinvw, Cinvww = WindowFunctionFourier.apply_window_covariance(ZONE, Cinv, xdata, thin=2, indexkred=kmask, masktriangle=masktriangle)

    ##################################################
    ##### Creating Taylor-egg ##############

    #nest = opa.join(OUTPATH, simname)
    nest = opa.join(OUTPATH, '%s-shifted'%simname)
    if not opa.isdir(nest): os.makedirs(nest)
    if rank ==0:
        print ('taylornest: %s'%nest)
        np.save(opa.join(nest,'k.npy'), kpred)

    z_pk, cosmoref, cosmodelta = import_cosmoref_from_DataFrameCosmosim(simname)
    Om_fid = 0.307115 #cosmoref[1]
    egging(nest, z_pk, cosmoref, cosmodelta, kpred, Om_fid, TableNkmu, 3, rank)
    
