#!/usr/bin/env python
from __future__ import print_function
import sys
from scipy import special, integrate
from scipy.interpolate import interp1d
import numpy as np
import os
import scipy
import sys


INPATH = 'input'
ZONE = sys.argv[1]


def check_if_multipoles_k_array(setk):
    return setk[int(len(setk) / 3)] == setk[0]


def Hubble(Om, z):
    return ((Om) * (1 + z)**3. + (1 - Om))**0.5


def DA(Om, z):
    r = integrate.quad(lambda x: 1. / Hubble(Om, x), 0, z)[0]
    return r / (1 + z)


def get_AP_param(z_pk, Om, Om_fid):
    qperp = DA(Om, z_pk) / DA(Om_fid, z_pk)
    qpar = Hubble(Om_fid, z_pk) / Hubble(Om, z_pk)
    return qperp, qpar


def W2D(x):
    return (2. * special.j1(x)) / x


def Hllp(l, lp, x):
    if l == 2 and lp == 0:
        return x ** 2 - 1.
    if l == 4 and lp == 0:
        return 1.75 * x**4 - 2.5 * x**2 + 0.75
    if l == 4 and lp == 2:
        return x**4 - x**2
    if l == 6 and lp == 0:
        return 4.125 * x**6 - 7.875 * x**4 + 4.375 * x**2 - 0.625
    if l == 6 and lp == 2:
        return 2.75 * x**6 - 4.5 * x**4 + 7. / 4. * x**2  # PZ: why 7/4
    if l == 6 and lp == 4:
        return x**6 - x**4
    else:
        return x * 0.


def fllp_IR(l, lp, k, q, Dfc):
    # IR q < k
    # q is an array, k is a scalar
    if l == lp:
        return (q / k) * W2D(q * Dfc) * (q / k)**l
    else:
        return (q / k) * W2D(q * Dfc) * (2. * l + 1.) / 2. * Hllp(max(l, lp), min(l, lp), q / k)


def fllp_UV(l, lp, k, q, Dfc):
    # UV q > k
    # q is an array, k is a scalar
    if l == lp:
        return W2D(q * Dfc) * (k / q)**l
    else:
        return W2D(q * Dfc) * (2. * l + 1.) / 2. * Hllp(max(l, lp), min(l, lp), k / q)


def dPuncorr(kout, fs=0.6, Dfc=0.43 / 0.6777):  # PZ: change kPS to kout
    """
    Compute the uncorrelated contribution of fiber collisions

    kPS : a cbird wavenumber output, typically a (39,) np array
    fs : fraction of the survey affected by fiber collisions
    Dfc : angular distance of the fiber channel Dfc(z = 0.55) = 0.43Mpc
    """
    dPunc = np.zeros((3, len(kout)))
    for l in [0, 2, 4]:
        dPunc[int(l / 2)] = - fs * np.pi * Dfc**2. * (2. * np.pi / kout) * (2. * l + 1.) / 2. * \
            special.legendre(l)(0) * (1. - (kout * Dfc)**2 / 8.) 
    return dPunc


def dPcorr_trust(kout, kPS, PS, ktrust=0.25, fs=0.6, Dfc=0.43 / 0.6777):  # PZ: added kout
    """
    Compute the correlated contribution of fiber collisions

    kPS : a cbird wavenumber output, typically a (39,) np array
    PS : a cbird power spectrum output, typically a (3, 39) np array
    ktrust : a UV cutoff
    fs : fraction of the survey affected by fiber collisions
    Dfc : angular distance of the fiber channel Dfc(z = 0.55) = 0.43Mpc
    """
    # import a very precise q array
    # q_ref = np.loadtxt('input_egg/Window_functions/k_LightConeHectorNGC.dat')

    q_ref = np.geomspace(min(kPS), ktrust, num=1024)  # PZ: construct the array within the interpolation range of PS
    # create log bin
    dq_ref = q_ref[1:] - q_ref[:-1]
    dq_ref = np.concatenate([[0], dq_ref])

    PS_interp = scipy.interpolate.interp1d(kPS, PS, axis=-1, bounds_error=False, fill_value='extrapolate')(q_ref)

    dPcorr = np.zeros(shape=(PS.shape[0], PS.shape[1], len(kout)))  # PZ: transformed PS to array of PS
    for j in range(PS.shape[1]):
        for l in [0, 2, 4]:
            for lp in [0, 2, 4]:
                # for i in range(len(kPS)):
                #    k=kPS[i]
                for i, k in enumerate(kout):
                    if lp <= l:
                        maskIR = (q_ref < k)
                        dPcorr[int(l / 2.), j, i] += - 0.5 * fs * Dfc**2 * np.einsum('q,q,q,q->', q_ref[maskIR],
                                                                                     dq_ref[maskIR], PS_interp[int(lp / 2.), j, maskIR], fllp_IR(l, lp, k, q_ref[maskIR], Dfc))
                    if lp >= l:
                        maskUV = ((q_ref > k) & (q_ref < ktrust))
                        dPcorr[int(l / 2.), j, i] += - 0.5 * fs * Dfc**2 * np.einsum('q,q,q,q->', q_ref[maskUV],
                                                                                     dq_ref[maskUV], PS_interp[int(lp / 2.), j, maskUV], fllp_UV(l, lp, k, q_ref[maskUV], Dfc))
    return dPcorr


def Pk_fibercollided(kout, kPS, PS, ktrust=0.25, fs=0.6, Dfc=0.43):
    """
    Apply fiber collision filter to the power spectrum and return fiber collided power spectrum

    kPS : a cbird wavenumber output, typically a (39,) np array
    PS : a cbird power spectrum output, typically a (3, 39) np array
    ktrust : a UV cutoff
    fs : fraction of the survey affected by fiber collisions. fs = 0.6 : BOSS DR12 value 
    Dfc : angular distance of the fiber channel Dfc(z = 0.55) = 0.43Mpc
    """
    PS_corr = dPcorr_trust(kout, kPS, PS, ktrust, fs, Dfc)
    return PS_corr

def load_window(simtype, zone, setkout, withmask=True, windowk=0.1):
    """
    Pre-load the window function to apply to the power spectrum by convolution in Fourier space

    Qll is an array of shape l,l',k',k where so that P_l(k) = \int dk' \sum_l' Q_{l l'}(k',k) P_l'(k')

    The original k on which Qll is evaluated is given by setk_or and k' by setkp_or.

    Inputs:
    ------
    setkout: the array of k on which the results should be evaluated.
    withmask: whether to only do the convolution over a small range around k
    windowk: the size of said range

    Output:
    ------
    PStransformed: the multipoles of power spectrum evaluted on setkout with the window function applied

    """
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]

    # Load window matrices
    Qll = np.load(opa.join(INPATH,'Window_functions/Qll_%s%s.npy'%(simtype, zone)))  
    setk_or = np.loadtxt(opa.join(INPATH,'Window_functions/kp_%s%s.txt'%(simtype, zone)))  
    setkp_or = np.loadtxt(opa.join(INPATH,'Window_functions/k_%s%s.dat'%(simtype, zone)))

    # Apply masking centered around the value of k
    if withmask:
        kpgrid,kgrid = np.meshgrid(setkp_or,setk_or,indexing='ij')
        mask = (kpgrid<kgrid+windowk)&(kpgrid>kgrid-windowk)
        Qll = np.einsum('lpkn,kn->lpkn',Qll,mask)
    
    # the spacing (needed to do the convolution as a sum)
    deltak = setkp_or[1:] - setkp_or[:-1]
    deltak = np.concatenate([[0],deltak])
    Qll_weighted = np.einsum('lpkn,k->lpkn',Qll,deltak)
    
    # Only keep value of setkp_or in the relevant range
    #maskred = ((setkp_or>setkout.min()-0.1*windowk)&(setkp_or<setkout.max()+windowk))
    maskred = ((setkp_or>setkout.min())&(setkp_or<setkout.max()+windowk))
    kpred = setkp_or[maskred]
    
    Qll_weighted_red = Qll_weighted[:,:,maskred,:]
    
    # Interpolate Qll(k) on setkout
    Qll_data = scipy.interpolate.interp1d(setk_or,Qll_weighted_red,axis=-1)(setkout)

    return kpred, Qll_data

def apply_window_PS(kpredwin, Qll_data, setPS, PS):
    """
    Apply the window function to the power spectrum by doing a convolution directly in fourier space
    
    Inputs:
    ------
    kpredwin: thin array of k' to perform the convolution integral
    Qll_data: pre-loaded window functions (evaluated at chosen k points) Q_{l l'}(k',k) x integral measure dk': see load_window(...)
    setPS: the array of k on which PS is evaluated
    PS: the multipoles of the PS (non concatenated), shape (3,len(setPS))

    Output:
    ------
    PStransformed: the multipoles of power spectrum evaluted on setkout with the window function applied

    """
    PS_red = scipy.interpolate.interp1d(setPS,PS,axis=-1,bounds_error=False,fill_value='extrapolate')(kpredwin)
    
    # (multipole l, multipole ' p, k, k' m) , (multipole ', power pectra s, k' m)
    PStransformed = np.einsum('lpkm,psk->lsm',Qll_data,PS_red)
    
    return PStransformed


def changetoAP(Pk, setkin, setkout, qperp, qpar, nbinsmu=500):

    muacc = np.linspace(0., 1., nbinsmu)

    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin): setkin = setkin[:len(setkin) / 3]
    if check_if_multipoles_k_array(setkout): setkout = setkout[:len(setkout) / 3]

    # Interpolate the multipoles
    Pkint = interp1d(setkin, Pk, axis=-1, kind='cubic', bounds_error=False, fill_value='extrapolate')

    # Define the grid with the right kmax and kmin and reshape into (k,mu)
    kgrid, mugrid = np.meshgrid(setkout, muacc, indexing='ij')

    # AP factors
    F = float(qpar / qperp)
    k = kgrid / qperp * (1 + mugrid**2 * (F**-2 - 1))**0.5
    mup = mugrid / F * (1 + mugrid**2 * (F**-2 - 1))**-0.5

    # Goes from the multipoles back to P(k,mu) and apply AP
    arrayLegendremup = np.array([special.legendre(0)(mup),
                                 special.legendre(2)(mup),
                                 special.legendre(4)(mup)])

    arrayLegendremugrid = np.array([2 * (2 * 0 + 1.) / (2 * qperp**2 * qpar) * special.legendre(0)(mugrid),
                                    2 * (2 * 2. + 1.) / (2 * qperp**2 * qpar) * special.legendre(2)(mugrid),
                                    2 * (2 * 4. + 1.) / (2 * qperp**2 * qpar) * special.legendre(4)(mugrid)])

    # Pkint(k).shape: (multipoles, power spectra, ks, mus): (lpkm)
    Pkmu = np.einsum('lpkm,lkm->pkm', Pkint(k), arrayLegendremup)

    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    # Pkmu.shape: (power spectra, ks, mus): (pkm)
    Integrandmu = np.einsum('pkm,lkm->lpkm', Pkmu, arrayLegendremugrid)

    Pk_AP = np.trapz(Integrandmu, x=mugrid, axis=-1)

    return Pk_AP


if __name__ == "__main__":

    simname = 'LightConeHector'
    ZONE = 'NGC'

    nrun = int(sys.argv[2])
    runs = int(sys.argv[3])

    z_pk, Om_fid = 0.32, 0.3153
    kmin = 0.01
    kmax = 0.3

    griddir = os.path.join(os.environ['GROUP_HOME'], 'GridsEFT') # path of the unprocessed grid
    outgrid = os.path.join(os.environ['GROUP_HOME'], 'modgrid')  # path for the processed grid
    gridname = "z0p32-A_s-h-omega_cdm-omega_b-n_s-Sum_mnu"       # gridname
    pspath = os.path.join('input', 'DataSims')                   # path of (k, power spectrum) data
    kPS, PSdata, _ = np.loadtxt(os.path.join(pspath, 'ps1D_%%s_data.dat') % (simname,ZONE), unpack=True)
    indexkred = np.where((kPS <= kmax) & (kPS >= kmin))[0]
    xdata = kPS[indexkred]
    kpred = xdata[:int(len(xdata) / 3)]

    thetatab = np.load(os.path.join(griddir, "Tablecoord_%s.npy" % gridname))  # Saved non-flattened, as (npar, sizegrid, ..., sizegrid)
    thetatab = np.transpose(thetatab.reshape((thetatab.shape[0], -1)))  # We flatten it, so there is correspondence to the PS

    # Load the PS grids directly in their format, flattened and concatenated for multipoles: (sizegrid**npar, len(k) * 3, columns)
    plingrid = np.load(os.path.join(griddir, "TablePlin_%s.npy" % gridname))
    ploopgrid = np.load(os.path.join(griddir, "TablePloop_%s.npy" % gridname))
    kfull = plingrid[0, :, 0]

    if check_if_multipoles_k_array(kfull): kfull = kfull[:int(len(kfull) / 3)]

    # Load window
    if 'GC' in ZONE:
        kpredwin, Qllwin = load_window(simname, ZONE, kpred, withmask=True, windowk=0.05)
    
    lenrun = int(len(thetatab) / runs)
    thetarun = thetatab[nrun * lenrun:(nrun + 1) * lenrun]
    sizered = len(thetarun)

    arrayred = thetarun[:sizered]

    allPlin = []
    allPloop = []
    for i, theta in enumerate(arrayred):
        idx = nrun * lenrun + i
        h = theta[1]
        omc = theta[2]
        omb = theta[3]
        qperp, qpar = get_AP_param(z_pk, (omc + omb) / (h*h), Om_fid)
        # Now put the PS in the order needed by the AP function.
        # Notice that we refer to the right parameters by idx
        Plin = np.swapaxes(plingrid[idx].reshape(3, len(kfull), plingrid.shape[-1]), axis1=1, axis2=2)[:, 1:, :]
        Ploop = np.swapaxes(ploopgrid[idx].reshape(3, len(kfull), ploopgrid.shape[-1]), axis1=1, axis2=2)[:, 1:, :]
        
        # AP effect
        if TableNkmu is not None:
            Plin = changetoAPbinning(Plin, kfull, kfull, qperp, qpar, TableNkmu)
            Ploop = changetoAPbinning(Ploop, kfull, kfull, qperp, qpar, TableNkmu)
        else:
            Plin = changetoAPnobinning(Plin, kfull, kfull, qperp, qpar, nbinsmu=100)
            Ploop = changetoAPnobinning(Ploop, kfull, kfull, qperp, qpar, nbinsmu=100)

        # Window function
        if 'GC' in ZONE:
            Plin = apply_window_PS(kpredwin, Qllwin, kfull, Plin)
            Ploop = apply_window_PS(kpredwin, Qllwin, kfull, Ploop)

        # fiber collisions
        ktrust = 0.2
        dPlin = dPcorr_trust(kpred, kfull, Plin, ktrust=ktrust)
        dPloop = dPcorr_trust(kpred, kfull, Ploop, ktrust=ktrust)

        allk = np.hstack([kpred, kpred, kpred])

        Plinconc = np.vstack([allk, np.concatenate(PlinW + dPlin, axis=-1)])
        Ploopconc = np.vstack([allk, np.concatenate(PloopW + dPloop, axis=-1)])
    
        idxcol = np.full([Plinconc.shape[0], 1], idx)
        allPlin.append(Plinconc.T)
        allPloop.append(Ploopconc.T)
        if ((i + 1) % 200 == 0):
            print("theta check: ", thetatab[idx], theta)
    np.save(os.path.join(outgrid, "Plin_run%s.npy" % str(nrun)), np.array(allPlin))
    np.save(os.path.join(outgrid, "Ploop_run%s.npy" % str(nrun)), np.array(allPloop))
