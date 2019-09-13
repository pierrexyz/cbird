# MCMC sampler for arXiv:1909.XXXXX
# the EFT power spectrum is evaluated with TBiRd.

# import sys
from scipy import stats, optimize, special
import time
import pandas as pd
# from scipy.interpolate import interp1d
# from multiprocessing import cpu_count
import numpy as np
import emcee
import Grid
import tbird
import os
import sys
###########################################
#  Globals ##############################
###########################################
THIS_PATH = os.path.dirname(__file__)

# Data paths
INPATH = os.path.abspath(os.path.join(THIS_PATH, 'input'))

GRIDPATH = os.path.abspath(THIS_PATH, 'GridsEFT')
OUTPATH = THIS_PATH

CHAINPATH = os.path.abspath(os.path.join(OUTPATH, 'chains'))

if not os.path.isdir(OUTPATH):
    raise Exception(OUTPATH + ' not there!')


###########################################
#  Functions  ###########################
###########################################

# rescale number density
nd = 4.e-4 
km = 0.7 
knl = 0.7 
k0p7 = 0.7
shotnoiseprior = 100. 
k2prior = 2. / 4.

C = 299792.458 # speed of light [km/s]
OG = 2.47282e-5 # omega_gamma = Omega_gamma h^2: normalized physical photon density today (T_cmb = 2.7255 (CLASS))
NUR = 3.046 #Nur: Number of ultra-relativistic species (CLASS):
ORAD = (1.+ NUR*7./8.*(4./11.)**(4./3.))*OG # omega_radiation
ZD = 1059.94 # Baryon-photon decoupling redshift (PLANCK 2018 TT,TE,EE+lowP+lensing (Table 4))
SIGMA_ZD = 0.30 # rd(zd+sigma)-dr(zd-sigma) < 0.2 sigma_rd: we take zd to be a delta function
RD = 147.09 # Sound horizon at decoupling [Mpc] (PLANCK 2018 TT,TE,EE+lowP+lensing (Table 4)):
SIGMA_RD = 0.26


def dPuncorr(kout, fs=0.6, Dfc=0.43 / 0.6777):
    """
    Compute the uncorrelated contribution of fiber collisions
    kPS : a cbird wavenumber output, typically a (39,) np array
    fs : fraction of the survey affected by fiber collisions
    Dfc : angular distance of the fiber channel Dfc(z = 0.55) = 0.43Mpc
    """
    dPunc = np.zeros((3, len(kout)))
    for l in [0, 2, 4]:
        dPunc[int(l / 2)] = (- fs * np.pi * Dfc**2. * (2. * np.pi / kout) * (2. * l + 1.) / 2. *
                             special.legendre(l)(0) * (1. - (kout * Dfc)**2 / 8.))
    return dPunc


def gelman_rubin_convergence(withinchainvar, meanchain, n, Nchains, ndim):
    """ Calculate Gelman & Rubin diagnostic
     1. Remove the first half of the current chains
     2. Calculate the within chain and between chain variances
     3. estimate your variance from the within chain and between chain variance
     4. Calculate the potential scale reduction parameter

    Inputs
    ------
        withinchainvar : array of the variances within each chains
        meanchain : array of the means within each chains
        n : length of the chains
        Nchains : number oc chains
        ndim : number of varied parameters
    Outputs
    ------
        The gelman rubin criteria
    """

    meanall = np.mean(meanchain, axis=0)
    W = np.mean(withinchainvar, axis=0)
    B = np.arange(ndim, dtype=np.float)
    for jj in range(0, ndim):
        B[jj] = 0.
    for jj in range(0, Nchains):
        B = B + n * (meanall - meanchain[jj])**2 / (Nchains - 1.)
    estvar = (1. - 1. / n) * W + B / n
    scalereduction = np.sqrt(estvar / W)

    return scalereduction


def get_Pi_for_marg(Ploop, kfull, kmin, kmax, b1, model = 1, bisp=None, masktriangle=None, withhex = True):

    nk = len(kfull)
    Pi = np.array([ 1.*(Ploop[:,16,:]+b1*Ploop[:,13,:]) / km**2 ])

    #if withhex: 
    Pi = np.concatenate(( Pi, np.array([ 1.*(Ploop[:,17,:]+b1*Ploop[:,14,:]) / km**2 ]) ))

    if model == 1: 
        Onel0 = np.array([np.ones(nk),np.zeros(nk),np.zeros(nk)]) # shot-noise mono
        kl2 = np.array([np.zeros(nk), kfull, np.zeros(nk)]) # k^2 quad
        Pi = np.concatenate(( Pi, 
            np.array([
                0.5* (Ploop[:,3,:]+b1*Ploop[:,7,:]),
                0.5* (Ploop[:,15,:]+b1*Ploop[:,12,:]) / knl**2,
                Onel0 * shotnoiseprior,
                0.25*kl2**2 / nd / km**2
                    ]) ))

    if model == 2:
        kl2 = np.array([np.zeros(nk), kfull, np.zeros(nk)]) # k^2 quad
        Pi = np.concatenate(( Pi, 
            np.array([
                0.5* (Ploop[:,3,:]+b1*Ploop[:,7,:]),
                0.5* (Ploop[:,15,:]+b1*Ploop[:,12,:]) / knl**2,
                0.25*kl2**2 / nd / km**2
                    ]) ))

    kmask = np.where((kfull > kmin) & (kfull < kmax+0.0001))[0]

    if withhex: Pi = Pi[:,:, kmask]
    else:       Pi = Pi[:,:2, kmask]
    Pi = Pi.reshape( (Pi.shape[0], -1) )

    if bisp is not None: # if with bisp, we can marginalized over the bisp shot noise P11 * P11 and b8^2/nd^2
        nparams = Pi.shape[0]
        nkpred = Pi.shape[1]
        nkbisp = sum(masktriangle)
        #newPi = np.zeros( shape=(nparams+1, nkpred+nkbisp) ) 
        newPi = np.zeros( shape=(nparams, nkpred+nkbisp) )
        newPi[:nparams, :nkpred] = Pi
        #newPi[-1, nkpred:] = bisp[4][masktriangle] * 0.00952 * shotnoiseprior * 10.
        #newPi[-1, nkpred:] = np.ones(nkbisp) * shotnoiseprior**2 * 4. * 100.
        Pi = 1.*newPi
    
    return Pi


def get_Covbi_for_marg(Pi_data, Cinv, sigma=200):
    Covbi = np.dot(Pi_data, np.dot(Cinv, Pi_data.T)) + 1. / sigma**2 * np.identity(Pi_data.shape[0])
    return Covbi


def import_simspec_from_DataFrameCosmosim(simtype):
    dfcosmo = pd.read_csv(os.path.join(INPATH, 'DataFrameCosmosims.csv'), index_col=0)
    series_cosmo = dfcosmo.loc[simtype]
    omega_bfid = series_cosmo.loc['omega_b']
    omega_cfid = series_cosmo.loc['Omega_m'] * series_cosmo.loc['h']**2 - series_cosmo.loc['omega_b']
    fb = omega_bfid / omega_cfid
    Omega_mfid = dfcosmo.loc[simtype, 'Omega_m']
    hfid = dfcosmo.loc[simtype, 'h']
    lnAsfid = dfcosmo.loc[simtype, 'lnAs']
    z_pk = dfcosmo.loc[simtype, 'z_pk']
    nsfid = dfcosmo.loc[simtype, 'ns']

    Om_AP = Omega_mfid  # Change DataFrameCosmosim for suppport Om_AP

    cosmoref = np.array([lnAsfid, Omega_mfid, hfid, omega_bfid, nsfid])
    return z_pk, cosmoref, Om_AP, fb


def import_data(simtype, boxnumber, kmin, kmax, kminbisp=0, kmaxbisp=0, ZONE=''):

    ##################################################
    # Loading covariance and Nkmu binning ##############
    if ('Challenge' in simtype) and ('Japan' not in simtype):
        if 'Quarter' in simtype:  # challenge quarter box
            Full_Cov = np.loadtxt(os.path.join(INPATH, 'Covariance/Cov%s%s.dat' % (simtype, boxnumber)))
        elif 'Hybrid' in simtype:  # challenge quarter covariance full power spectrum
            simtype = 'Challenge'
            TableNkmu = np.loadtxt(os.path.join(INPATH, 'Binning/Nkmu%s%s.dat' % (simtype, boxnumber))).T
            Full_Cov = np.loadtxt(os.path.join(INPATH, 'Covariance/Cov%s%s.dat' % ('ChallengeQuarter', boxnumber)))
        else:  # challenge full box
            # print('Using quarter covariance divided by 4.25 instead of full')
            Full_Cov = np.loadtxt(os.path.join(INPATH,
                                               'Covariance/Cov%s%s.dat' %
                                               ('ChallengeQuarter', boxnumber))) / 4.25
        TableNkmu = np.loadtxt(os.path.join(INPATH, 'Binning/Nkmu%s%s.dat' % (simtype, boxnumber))).T
    elif 'Japan' in simtype:
        TableNkmu = None
        Full_Cov = np.loadtxt(os.path.join(INPATH, 'Covariance/CovFull%s_%s.dat' % (simtype, boxnumber)))
    else:
        TableNkmu = None
        if "mean" in boxnumber:
            print('mean')
            Full_Cov = np.loadtxt(os.path.join(INPATH, 'Covariance/Cov%s%sdata.dat' % (simtype, ZONE))) / 16.
        else:
            Full_Cov = np.loadtxt(os.path.join(INPATH, 'Covariance/Cov%s%sdata.dat' % (simtype, ZONE)))

    #################################################
    # Loading power spectrum data ####
    kPS, PSdata, _ = np.loadtxt(os.path.join(INPATH, 'DataSims/ps1D_%s%s_%s.dat' %
                                             (simtype, ZONE, boxnumber)), unpack=True)

    indexkred = np.argwhere((kPS <= kmax) & (kPS >= kmin))[:, 0]
    xdata = kPS[indexkred]
    ydata = PSdata[indexkred]

    kmask = np.array([False] * len(kPS))
    kmask[indexkred] = True
    kmask = kmask[:int(len(kmask) / 3)]

    Covred = Full_Cov[indexkred.reshape((len(indexkred), 1)), indexkred]
    kpred = xdata[:len(xdata) / 3]
    N = int(2 * len(ydata) / 3)
    Covred = Full_Cov[indexkred.reshape((len(indexkred), 1)), indexkred][:N, :N]
    Cinv = np.linalg.inv(Covred)
    xdata = xdata[:N]
    ydata = ydata[:N]
    chi2data = np.dot(ydata, np.dot(Cinv, ydata))
    Cinvdata = np.dot(ydata, Cinv)

    return kpred, chi2data, Cinvdata, Cinv, TableNkmu, xdata, ydata, np.sqrt(np.diag(Covred))


def rs(Om, h, f_fid):
    om = Om * h**2
    ob = om * f_fid / (f_fid + 1.)  # f_fid: fiducial ratio omega_b/omega_c
    R = 0.75 * ob / OG
    result = ((2. * C / 100. / np.sqrt(3. * R * om)) *
              np.log((np.sqrt(1. + ZD + R) + np.sqrt((1. + ZD) * R * ORAD / om + R)) / np.sqrt(1. + ZD) /
                     (1. + np.sqrt(R * ORAD / om))))
    return result


def Hubble(Om, z):
    return ((Om) * (1 + z)**3. + (1 - Om))**0.5


def DA(Om, z):
    r = scipy.integrate.quad(lambda x: 1. / Hubble(Om, x), 0, z)[0]
    return r / (1 + z)


def check_if_multipoles_k_array(setk):
    """Check if we have 3 identical sets of k in the same file"""
    return setk[int(len(setk) / 2)] == setk[0]


def computePS(cvals, plin, ploop, setkin, kmin, kmax, sigsq=0, Puncorr=0):
    plin0, plin2 = plin
    ploop0, ploop2 = ploop[:, :18, :]
    b1, c2, b3, c4, b5, b6, b7, b8, b9, b10, b11, e1, e2, e3 = cvals

    b2 = (c2 + c4) / np.sqrt(2.)
    b4 = (c2 - c4) / np.sqrt(2.)

    # the columns of the Ploop data files.
    cvals = np.array([1, b1, b2, b3, b4, b1 * b1, b1 * b2, b1 * b3, b1 * b4, b2 * b2, b2 * b4, b4 * b4,
                      b1 * b5 / knl**2, b1 * b6 / km**2, b1 * b7 / km**2, b5 / knl**2, b6 / km**2, b7 / km**2])

    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin) / 3]

    kmask = np.where((setkin >= kmin) & (setkin <= kmax))[0]

    P0 = (np.dot(cvals, ploop0) +
          plin0[0] + b1 * plin0[1] + b1 * b1 * plin0[2]
          - 2 * (-b1 + b2 + b4)**2 * sigsq
          + b8 / nd + b9 / nd / km**2 * setkin**2)
    P2 = (np.dot(cvals, ploop2) +
          plin2[0] + b1 * plin2[1] + b1 * b1 * plin2[2]
          + b10 / nd / km**2 * setkin**2)
    # P4 = (np.dot(cvals, ploop4) +
    #       plin4[0] + b1 * plin4[1] + b1 * b1 * plin4[2]
    #       + e3 * (b1 * ploop4e1b1 + ploop4e1))

    if Puncorr is not 0:
        P0 += Puncorr[0]
        P2 += Puncorr[1]

    return np.array([P0[kmask], P2[kmask]])


def match_para(theta, free_para, fix_para):

    value_array = np.arange(len(free_para), dtype=np.float)
    counter = 0
    for i in range(len(free_para)):
        if free_para[i] is True:
            value_array[i] = theta[counter]
            counter += 1
        else:
            value_array[i] = fix_para[i]

    return value_array


def lnprior(theta, free_para, fix_para, bounds, withPlanck, cosmoref):
    value_array = match_para(theta, free_para, fix_para)

    lnAs, Om, h, omb, ns, Summnu, b1, c2, b3, c4, b5, b6, b7, b8, b9, b10, b11, e1, e2, e3 = value_array
    ombtrue = cosmoref[3]
    nstrue = cosmoref[4]

    withinprior = True
    for i in range(6 + 2):
        withinprior = (withinprior) and (bounds[i][0] <= value_array[i] <= bounds[i][1])  # cosmo + b1, c2
    if withinprior:
        if withPlanck:
    	    prior = - 0.5 * ((b3 / 2.)**2 + (c4 / 2.)**2 + (b5 / 2.)**2 + (b6 / 8.)**2 + (b7 / 4.)**2 +
    		    (b8 / 400.)**2 + (b10 / 2.)**2  + ((omb - ombtrue) / 0.00015)**2 + ((ns - nstrue)/0.0044)**2)
        else:
    	    prior = - 0.5 * ((b3 / 2.)**2 + (c4 / 2.)**2 + (b5 / 2.)**2 + (b6 / 8.)**2 + (b7 / 4.)**2 +
    		    (b8 / 400.)**2 + (b10 / 2.)**2 + ((omb - ombtrue) / 0.0007425)**2)
        return prior
    else:
        return -np.inf


def lnlike(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, cosmoref,
           linder, loopder, sigma_prior=100, marg_gaussian=False, model=2):

    lnAs, Om, h, omb, ns, Summnu, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, e1, e2, e3 = match_para(
        theta, free_para, fix_para)
    # omega_c = Om * h * h / (1. + f_fid)
    # omega_b = f_fid * omega_c
    omega_c = Om * h * h - omb
    omega_b = omb
    dtheta = np.array((np.exp(lnAs) * 1e-10 - Grid.valueref[0], h - Grid.valueref[1],
                       omega_c - Grid.valueref[2], omega_b - Grid.valueref[3],
                       ns - Grid.valueref[4], Summnu - Grid.valueref[5]))

    Plin = tbird.get_PSTaylor(dtheta, linder)
    Ploop = tbird.get_PSTaylor(dtheta, loopder)

    kfull = Plin[0, :, 0]
    if check_if_multipoles_k_array(kfull):
        kfull = kfull[:int(len(kfull) / 2)]
    Puncorr = dPuncorr(kfull)

    Plin = np.swapaxes(Plin.reshape(2, len(kfull), Plin.shape[-1]), axis1=1, axis2=2)[:, 1:, :]
    Ploop = np.swapaxes(Ploop.reshape(2, len(kfull), Ploop.shape[-1]), axis1=1, axis2=2)[:, 1:, :]

    bs = np.array([b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, b11, e1, e2, e3])

    modelX = computePS(bs, Plin, Ploop, kfull, kmin, kmax, Puncorr=Puncorr).reshape(-1)

    if marg_gaussian:
        Pi = get_Pi_for_marg(Ploop, kfull, kmin, kmax, b1, model)

        Pi = Pi.reshape((Pi.shape[0], -1))

        Covbi = get_Covbi_for_marg(Pi, Cinv, sigma=sigma_prior)
        Cinvbi = np.linalg.inv(Covbi)
        vectorbi = np.dot(modelX, np.dot(Cinv, Pi.T)) - np.dot(Cinvdata, Pi.T)
        chi2nomar = np.dot(modelX, np.dot(Cinv, modelX)) - 2. * np.dot(Cinvdata, modelX) + chi2data
        chi2mar = -np.dot(vectorbi, np.dot(Cinvbi, vectorbi)) + np.log(np.linalg.det(Covbi))
        chi2 = chi2mar + chi2nomar - 5 * np.log(2. * np.pi)
    else:
        chi2 = np.dot(modelX - ydata, np.dot(Cinv, modelX - ydata))

    return -0.5 * chi2


def lnprob(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, cosmoref,
           linder, loopder, sigma_prior=100, marg_gaussian=False, model=1):

    lp = lnprior(theta, free_para, fix_para, bounds, withPlanck, cosmoref)

    if not np.isfinite(lp):
        dummy = -np.inf
    else:
        dummy = lp + lnlike(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, cosmoref,
                            linder, loopder, sigma_prior=sigma_prior, marg_gaussian=marg_gaussian, model=model)
    return dummy


###########################################
#  Main program  ########################
###########################################

if __name__ == "__main__":

    # print("started")
    simtype = str(sys.argv[1])
    boxnumber = sys.argv[2]
    kmin = float(sys.argv[3])
    kmax = float(sys.argv[4])
    marg_gaussian = int(sys.argv[5])  # 0 or 1 (True of False)
    planckchain = int(sys.argv[6])  # Boolean, whether to use or not the Planck fiducial
    if planckchain == 1:
        withPlanck = True
    else:
        withPlanck = False
    model = int(sys.argv[7])
    ZONE = str(sys.argv[8])
    gridname = str(sys.argv[9])
    if "Challenge" in simtype:
        ZONE = ''
    if "Nseries" in simtype:
        ZONE = ''
    runtype = simtype + ZONE
    if withPlanck:
        runtype += 'wPlanck'
    if marg_gaussian:
        runtype += 'gaussMarg'

    simname = simtype
    if 'LightConeHector' in simtype:
        simname = 'LightConeHector'

    withBisp = False
    kminbisp = 0.04
    kmaxbisp = 0  # float(sys.argv[6]) ########################### no current suppport for bisp
    if kmaxbisp > 0:
        withBisp = True
        print('kmax bispectrum is bigger than zero, withBisp is %s' % withBisp)

    ##################################################
    # Loading power spectrum data, covariance and Nkmu binning | loading simulation specification
    if "Challenge" in simtype:
        if 'D' in boxnumber:
            simname = 'ChallengeD'
            # RD = 147.253
        elif 'Japan' in simtype:
            simname = simtype
        else:
            simname = 'ChallengeA'
    print(simtype)
    print(boxnumber)
    kpred, chi2data, Cinvdata, Cinv, TableNkmu, xdata, ydata, yerror = import_data(
        simtype, boxnumber, kmin, kmax, kminbisp, kmaxbisp, ZONE)

    simspec = import_simspec_from_DataFrameCosmosim(simname)
    z_pk, cosmosim, Om_fid, fb = simspec

    ##################################################
    # Parameters to minimize ######
    if marg_gaussian:
        free_para = [True, True, True, True, True, True,
                     True, True, False, withBisp,
                     False, False, False,
                     withBisp, False, False, False,
                     False, False, False]
        a = 2.0
    else:
        a = 1.3
        if model == 1:
            free_model = [False, True, withBisp, False, True, False]
        elif model == 2:
            free_model = [False, True, withBisp, False, False, False]
        else:
            free_model = [False, False, withBisp, False, False, False]

        free_para = [True, True, True, True, True, True,
                     True, True, True,
                     False, True, True, False, True] + free_model

    runtype += '-model%s' % model
    ndim = sum(free_para)
    fix_para = np.concatenate((list(cosmosim) + [0.2], [2.] + 13 * [0])
                              )  # if free_para is false read the value in fix_para

    lnAsmin, lnAsmax, Ommin, Ommax, hmin, hmax = (1.5, 4., 0.2, 0.4, 0.5, 0.85)
    ombmin, ombmax, nsmin, nsmax, Summnumin, Summnumax = (0.01, 0.03, 0.87, 1.07, 0.059, 1.)  # DES priors
    print("got grid!")

    ##################################################
    # Uniform prior on the As, h, Om and b_i ######
    priorsup = 10.
    priorgauss = 4.
    bfmintab = np.concatenate(([lnAsmin, Ommin, hmin, ombmin, nsmin, Summnumin], [0., -4.] +
                               12 * [-priorsup]))  # We require b_1 > 0
    bfmaxtab = np.concatenate(([lnAsmax, Ommax, hmax, ombmax, nsmax, Summnumax], [4., 4.] + 12 * [priorsup]))
    bounds = zip(bfmintab, bfmaxtab)

    # Guess for the \sigma, to help with the initial position of walkers #####
    onesigma = np.array([0.1, 0.03, 0.03, 0.003, 0.01, 0.0001, 0.5, 0.5, 0.5,
                         0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 0.5, 0.5, 0.5])

    # find max likelihood
    tryini = fix_para   # initial guess

    linder = np.load(os.path.abspath(
        os.path.join(GRIDPATH, 'DerPlin%s.npy' % gridname)), allow_pickle=True)
    loopder = np.load(os.path.abspath(
        os.path.join(GRIDPATH, 'DerPloop%s.npy' % gridname)), allow_pickle=True)

    def lnprobloaded(theta):
        return lnprob(theta, chi2data, Cinvdata, Cinv, kmin, kmax,
                      free_para, fix_para, bounds, cosmosim,
                      linder, loopder,
                      sigma_prior=priorgauss, marg_gaussian=marg_gaussian,
                      model=model)
    print("Loaded probability")

    # def chi2(theta):
    #     return -2. * lnprobloaded(theta)

    # result = optimize.minimize(chi2, tryini, method='SLSQP', bounds=bounds, options={'maxiter': 1000})

    # bestfit = match_para(result["x"], free_para, fix_para)
    # minchi2 = result["fun"]
    # dof = len(xdata)-sum(free_para)
    # pvalue = 1. - stats.chi2.cdf(minchi2, dof)

    # print('minchi2 = ' + str(minchi2), dof)
    # np.savetxt(os.path.join(OUTPATH, "minchi2%sbox_%skmax_%s.txt") %
    #            (runtype, boxnumber, kmax), np.concatenate([bestfit, [minchi2, dof, pvalue]]))

    bestfit = match_para(tryini, free_para, fix_para)

    ###################################
    # run MCMC #######################
    ###################################

    # Start MCMC
    t0 = time.time()
    Nchains = 4
    nwalkers = 12 * ndim
    # nthreads = cpu_count() / Nchains
    nthreads = 1
    # Set up the sampler.
    pos = []
    sampler = []
    # rstate = np.random.get_state()

    t_try = time.time()
    for jj in range(0, Nchains):
        initialpos = []
        for ii in xrange(nwalkers):
            accepted = False
            t_try = time.time()
            while (not accepted):
                trialfiducial = np.random.normal(loc=bestfit[free_para], scale=onesigma[free_para])
                accepted = np.isfinite(lnprior(trialfiducial, free_para, fix_para, bounds, withPlanck))
            if accepted:
                initialpos.append(trialfiducial)
        pos.append(initialpos)
        print("Found initial pos in %s seconds!" % str(time.time() - t_try))
        sampler.append(emcee.EnsembleSampler(nwalkers, ndim, lnprobloaded, a=a, threads=nthreads))
        time.sleep(5)

    ##################################################
    epsilon = 0.005
    print("Defined epsilon")
    withinchainvar = np.zeros((Nchains, ndim))
    meanchain = np.zeros((Nchains, ndim))
    scalereduction = np.arange(ndim, dtype=np.float)
    for jj in range(0, ndim):
        scalereduction[jj] = 2.
    minlength = 5000
    chainstep = minlength
    itercounter = 0
    loopcriteria = 1
    ichaincheck = 50
    ithin = 1
    t1 = time.time()

    while loopcriteria:
        itercounter = itercounter + chainstep
        print("chain length  = ", itercounter, " minlength  = ", minlength)
        samplesJG = []
        for jj in range(0, Nchains):
            c = 0
            for result in sampler[jj].sample(pos[jj], iterations=chainstep, storechain=True, thin=ithin):
                pos[jj] = result[0]
                chainchi2 = -2. * result[1]
                # rstate = result[2]
            # we do the convergence test on the second half of the current chain (itercounter/2)
            chainsamples = sampler[jj].chain[:, itercounter / 2:, :].reshape((-1, ndim))
            withinchainvar[jj] = np.var(chainsamples, axis=0)
            meanchain[jj] = np.mean(chainsamples, axis=0)
            samplesJG.append(chainsamples)
        scalereduction = gelman_rubin_convergence(withinchainvar, meanchain, itercounter / 2, Nchains, ndim)
        print("scalereduction  =  ", scalereduction)
        loopcriteria = 0
        for jj in range(0, ndim):
            if np.absolute(1 - scalereduction[jj]) > epsilon:
                loopcriteria = 1

        chainstep = ichaincheck

    print("Done.")

    ################################
    stringname = 'omb'
    ##################################################
    # Saving results #######################
    samplerchain = sampler[0].chain.reshape((-1, ndim))
    lnprobchain = sampler[0].lnprobability.reshape(-1)
    np.save(os.path.join(CHAINPATH, "samplerchain%sbox_%skmin_%skmax_%s_%s") %
            (runtype, boxnumber, kmin, kmax, stringname), samplerchain)
    np.save(os.path.join(CHAINPATH, "lnlikechain%sbox_%skmin_%skmax_%s_%s") %
            (runtype, boxnumber, kmin, kmax, stringname), lnprobchain)

    bestfit = match_para(samplerchain[lnprobchain.argmax(), :], free_para, fix_para)
    minchi2 = -2 * lnprobchain.max()
    dof = len(xdata) - sum(free_para)
    pvalue = 1. - stats.chi2.cdf(minchi2, dof)
    np.savetxt(os.path.join(CHAINPATH, "bestfit-mcmc%sbox_%skmin_%skmax_%s_%s.txt") %
               (runtype, boxnumber, kmin, kmax, stringname), np.concatenate([bestfit, [minchi2, dof, pvalue]]), fmt='%.3f')

    # Compute the quantiles ; columns: expectation values, +sigma, -sigma, meansigma
    mcmc_array = np.array(list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0], 0.5 * (v[2] - v[0])),
                                   zip(*np.percentile(samplerchain, [16, 50, 84], axis=0)))))

    ncosmo = 5
    theoryerror = np.zeros(ncosmo)
    meancosmo = mcmc_array[:ncosmo, :]
    for i, (t, m) in enumerate(zip(cosmosim, meancosmo)):
        if t > m[0] + m[1]:
            theoryerror[i] = t - (m[0] + m[1])
        elif t < m[0] - m[2]:
            theoryerror[i] = t - (m[0] - m[2])
    results = np.c_[cosmosim[:ncosmo], bestfit[:ncosmo], meancosmo, theoryerror]
    print(results)

    np.savetxt(os.path.join(CHAINPATH, "mcmcresults%sbox_%skmin_%skmax_%s_%s.txt") %
               (runtype, boxnumber, kmin, kmax, stringname), results, fmt='%.3f')

    print("Mean acceptance fraction: ", np.mean(sampler[0].acceptance_fraction))
    np.savetxt(os.path.join(CHAINPATH, 'AcceptanceFr%sbox_%skmin_%skmax_%s_%s.dat' %
                            (runtype, boxnumber, kmin, kmax, stringname)), sampler[0].acceptance_fraction, fmt='%.2f')
