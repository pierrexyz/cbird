import sys
from scipy import stats
import time
import pandas as pd
from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal
from multiprocessing import cpu_count
import scipy.optimize as op
import numpy as np
import emcee
import scipy.stats
import os.path as opa
import os
import fiber_col as fc

THIS_PATH  =  opa.dirname(__file__)
INPATH = opa.abspath(opa.join(THIS_PATH,'input/'))
OUTPATH = opa.abspath(opa.join(THIS_PATH,'output/'))
if not opa.isdir(OUTPATH): os.makedirs(OUTPATH)
CHAINPATH = opa.abspath(opa.join(OUTPATH,'chains/'))
if not opa.isdir(CHAINPATH): os.makedirs(CHAINPATH)

covPlanck2018 = np.array([[ 2.02650949e-04, -2.84195795e-05,  2.14797765e-05],
        [-2.84195795e-05,  5.43801869e-05, -3.93947626e-05],
        [ 2.14797765e-05, -3.93947626e-05,  2.90327720e-05]])

centralPlanck2018 = np.array([3.045, 0.3166, 0.6727])

C = 299792.458 # speed of light [km/s]
OG = 2.47282e-5 # omega_gamma = Omega_gamma h^2: normalized physical photon density today (T_cmb = 2.7255 (CLASS))
NUR = 3.046 #Nur: Number of ultra-relativistic species (CLASS):
ORAD = (1.+ NUR*7./8.*(4./11.)**(4./3.))*OG # omega_radiation
ZD = 1059.94 # Baryon-photon decoupling redshift (PLANCK 2018 TT,TE,EE+lowP+lensing (Table 4))
SIGMA_ZD = 0.30 # rd(zd+sigma)-dr(zd-sigma) < 0.2 sigma_rd: we take zd to be a delta function
RD = 147.09 # Sound horizon at decoupling [Mpc] (PLANCK 2018 TT,TE,EE+lowP+lensing (Table 4)):
SIGMA_RD = 0.26

def rs(Om, h , fb): # sound horizon at decoupling
    om = Om * h**2
    ob = om * fb/(fb+1.) # fb: ratio omega_b/omega_c
    R = 0.75 * ob/OG
    result = 2.*C/100./ np.sqrt(3.*R*om)* np.log( ( np.sqrt(1.+ZD+R) + np.sqrt((1.+ZD)*R*ORAD/om+R) ) / np.sqrt(1.+ZD) / (1.+np.sqrt(R*ORAD/om)) )
    return result

# rescale number density
nd = 4.5e-4
km = 0.7
knl = 0.7
shotnoiseprior = 100.

def gelman_rubin_convergence(withinchainvar, meanchain, n, Nchains, ndim):
    meanall  =  np.mean(meanchain, axis = 0)
    W  =  np.mean(withinchainvar, axis = 0)
    B  =  np.arange(ndim,dtype = np.float)
    for jj in range(0, ndim):
        B[jj]  =  0.
    for jj in range(0, Nchains):
        B  =  B + n*(meanall - meanchain[jj])**2/(Nchains-1.)
    estvar  =  (1. - 1./n)*W + B/n
    scalereduction  =  np.sqrt(estvar/W)

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
    Covbi = np.dot(Pi_data, np.dot(Cinv,Pi_data.T)) + 1./sigma**2 * np.identity(Pi_data.shape[0])
    return Covbi

def import_simspec_from_DataFrameCosmosim(simtype):
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
    
    return z_pk, cosmoref, Omega_mfid, fb

def import_data(simtype, boxnumber, kmin, kmax, kminbisp = 0, kmaxbisp = 0, ZONE = '', withhex = True):
    
    ##################################################
    ##### Loading covariance and Nkmu binning ##############
    if 'Challenge' in simtype:
        if 'Quarter' in simtype: # challenge quarter box
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype, boxnumber)))
            TableNkmu = np.loadtxt(opa.join(INPATH,'Binning/Nkmu%s%s.dat' % (simtype, boxnumber))).T
        elif 'Hybrid' in simtype: # challenge quarter covariance full power spectrum
            simtype = 'Challenge'
            TableNkmu = np.loadtxt(opa.join(INPATH,'Binning/Nkmu%s%s.dat' % (simtype, boxnumber))).T
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%('ChallengeQuarter', boxnumber)))
        elif 'PT' in simtype: # Japan PT challenge
            TableNkmu = None
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype, boxnumber)))
        else: # challenge full box
            #print('Using quarter covariance divided by 4.25 instead of full')
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%('ChallengeQuarter',boxnumber))) / 4.25
            TableNkmu = np.loadtxt(opa.join(INPATH,'Binning/Nkmu%s%s.dat' % (simtype, boxnumber))).T
    elif 'NseriesPatchy' in simtype:
        TableNkmu = None
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/CovNseriesPatchy.dat')) 
    else:
        TableNkmu = None
        if 'mean' in boxnumber:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%smean.dat'%(simtype,ZONE)))
        elif 'LightConeHector' in simtype or 'lowz' in simtype:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%sdatav5crossmultipole.dat'%(simtype,ZONE)))
        else:
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%sdata.dat'%(simtype,ZONE)))

    Bispdata = [] 
    masktriangle = []
    
    #################################################
    ##### Loading power spectrum data ####
    try:    kPS, PSdata, _ = np.loadtxt(opa.join(INPATH, 'DataSims/ps1D_%s%s_%s.dat' % (simtype, ZONE, boxnumber)), unpack=True)
    except: kPS, PSdata = np.loadtxt(opa.join(INPATH, 'DataSims/ps1D_%s%s_%s.dat' % (simtype, ZONE, boxnumber)), unpack=True)

    indexkred = np.argwhere((kPS < kmax) & (kPS > kmin))[:, 0]
    xdata = kPS[indexkred]
    ydata = PSdata[indexkred]

    kpred = xdata[:len(xdata)/3]

    if withhex:
        Covred = Full_Cov[indexkred.reshape((len(indexkred), 1)), indexkred]
    else:
        N = 2*len(ydata)/3
        Covred = Full_Cov[indexkred.reshape((len(indexkred), 1)), indexkred][:N, :N]
        xdata = xdata[:N]
        ydata = ydata[:N]

    masktrianglegrid = None

    if kmaxbisp > 0:
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s_Bispv5.dat'%(simtype,ZONE)))   
        #Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov_LightConeHectorNGC_Bispv5_nohex.dat'))
        Q1, Q2, Q3, Bispdata = np.loadtxt(opa.join(INPATH,'DataSims/Bispred_LightConeHector_%s_%s.dat'%(ZONE,boxnumber))).T
        KMAXBISP = 0.15 # kmax grid
        R1,R2,R3 = Q1[(Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)], Q2[(Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)], Q3[(Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)]
        masktriangle = (Q1>=kminbisp)&(Q1<=kmaxbisp)&(Q1<=Q2+Q3)&(Q1>=abs(Q2-Q3))&(Q2>=kminbisp)&(Q2<=kmaxbisp)&(Q3>=kminbisp)&(Q3<=kmaxbisp)
        masktrianglegrid = (R1>=kminbisp)&(R1<=kmaxbisp)&(R1<=R2+R3)&(R1>=abs(R2-R3))&(R2>=kminbisp)&(R2<=kmaxbisp)&(R3>=kminbisp)&(R3<=kmaxbisp)
        indextriangle = np.argwhere(masktriangle)[:,0]+kPS.shape[0]
        
        if withhex: indexkred = np.concatenate([indexkred, indextriangle])
        else: indexkred = np.concatenate([indexkred[:2*len(indexkred)/3], indextriangle])
        
        xdata = np.concatenate([xdata, np.arange(len(Bispdata[masktriangle])) ])
        ydata = np.concatenate([ydata, Bispdata[masktriangle]])
        Covred = Full_Cov[indexkred.reshape((len(indexkred), 1)), indexkred]

    Cinv = np.linalg.inv(Covred)
    chi2data = np.dot(ydata, np.dot(Cinv, ydata))
    Cinvdata = np.dot(ydata, Cinv)

    return kpred, chi2data, Cinvdata, Cinv, TableNkmu, masktrianglegrid, xdata, ydata, np.sqrt(np.diag(Covred))

def check_if_multipoles_k_array(setk):
    """Check if we have 3 identical sets of k in the same file"""
    return setk[int(len(setk) / 3)] == setk[0]


def get_grid(gridname, nbinsAs=80, nbinsOm=80, nbinsh=70, withBisp = False):

    thetatab = np.load(os.path.abspath(
        os.path.join(INPATH, 'GridsEFT/Tablecoord%s.npy' % gridname)))

    theta3D = thetatab.reshape((nbinsAs, nbinsOm, nbinsh, 3))

    lnAstab = theta3D[:, 0, 0, 0]
    Omtab = theta3D[0, :, 0, 1]
    htab = theta3D[0, 0, :, 2]

    lnAsmin = lnAstab.min()
    lnAsmax = lnAstab.max()
    Ommin = Omtab.min()
    Ommax = Omtab.max()
    hmin = htab.min()
    hmax = htab.max()

    TablePlin = np.load(os.path.abspath(
        os.path.join(INPATH, 'GridsEFT/TablePlin%s.npy' % gridname)))
    TablePloop = np.load(os.path.abspath(
        os.path.join(INPATH, 'GridsEFT/TablePloop%s.npy' % gridname)))

    Plininterp = scipy.interpolate.RegularGridInterpolator((lnAstab, Omtab, htab), TablePlin.reshape(
        (nbinsAs, nbinsOm, nbinsh, TablePlin.shape[-2], TablePlin.shape[-1])))
    Ploopinterp = scipy.interpolate.RegularGridInterpolator((lnAstab, Omtab, htab), TablePloop.reshape(
        (nbinsAs, nbinsOm, nbinsh, TablePloop.shape[-2], TablePloop.shape[-1])))

    if not withBisp:
        interpolations = [Plininterp, Ploopinterp]
    else:
        TableBisp = np.load(opa.abspath(opa.join(INPATH,'GridsEFT/TableBisp%s.npy'%gridname)))

        theta3D = thetatab.reshape((70, 70, 120, 3))

        lnAstab = theta3D[:, 0, 0, 0]
        Omtab = theta3D[0, :, 0, 1]
        htab = theta3D[0, 0, :, 2]

        lnAsmin2 = lnAstab.min()
        lnAsmax2 = lnAstab.max()
        Ommin2 = Omtab.min()
        Ommax2 = Omtab.max()
        hmin2 = htab.min()
        hmax2 = htab.max()

        if lnAsmin2 > lnAsmin: lnAsmin = lnAsmin2
        if lnAsmax2 < lnAsmax: lnAsmax = lnAsmax2
        if Ommin2 > Ommin: Ommin = Ommin2
        if Ommax2 < Ommax: Ommax = Ommax2
        if hmin2 > hmin: hmin = hmin2
        if hmax2 < hmax: hmax = hmax2

        Bispinterp = scipy.interpolate.RegularGridInterpolator((lnAstab,Omtab,htab),TableBisp.reshape((70, 70, 120,TableBisp.shape[-2],TableBisp.shape[-1])))
        interpolations = [Plininterp,Ploopinterp,Bispinterp]

    return lnAsmin, lnAsmax, Ommin, Ommax, hmin, hmax, interpolations

def computePS(bs, plin, ploop, kin, kmin, kmax, withhex = True, Puncorr = 0):
    plin0, plin2, plin4 = plin
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10 = bs
    ploop0, ploop2, ploop4 = ploop[:,:18,:]

    b2 = 1./np.sqrt(2.) * (c2 + c4)
    b4 = 1./np.sqrt(2.) * (c2 - c4)

    cvals = np.array([1., b1, b2, b3, b4, b1*b1, b1*b2, b1*b3, b1*b4, b2*b2, b2*b4, b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2])

    P0 = np.dot(cvals, ploop0) + plin0[0]+b1*plin0[1]+b1*b1*plin0[2] + b8 + b9/nd/km**2 * kin**2
    P2 = np.dot(cvals, ploop2) + plin2[0]+b1*plin2[1]+b1*b1*plin2[2]          + b10/nd/km**2 * kin**2
    P4 = np.dot(cvals, ploop4) + plin4[0]+b1*plin4[1]+b1*b1*plin4[2]                                   

    if Puncorr is not 0:
        P0 += Puncorr[0]
        P2 += Puncorr[1]

    kmask = np.where((kin > kmin) & (kin < kmax +0.0001))[0]

    if withhex: return np.array([P0[kmask], P2[kmask], P4[kmask]])
    else:       return np.array([P0[kmask], P2[kmask]])

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

def lnlike(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, Om_fid,
           interpolation_grid, sigma_prior = 100, marg_gaussian=False, model = 1, 
           masktriangle = None, withhex = True, Puncorr = 0):

    lnAs,Om,h,b1,c2,b3,c4,b5,b6,b7,b8,b9,b10,b11,e1 = match_para(theta, free_para, fix_para)

    # Import the power spectra interpolators on the grid
    if masktriangle is None: Plininterp, Ploopinterp = interpolation_grid
    else:   Plininterp,Ploopinterp,Bispinterp= interpolation_grid 

    Plin = Plininterp((lnAs, Om, h))
    Ploop = Ploopinterp((lnAs, Om, h))
    kfull = Plin[:, 0]
    if check_if_multipoles_k_array(kfull): kfull = kfull[:int(len(kfull) / 3)]

    Plin = np.swapaxes(Plin.reshape(3, len(kfull), Plin.shape[1]), axis1=1, axis2=2)[:, 1:, :]
    Ploop = np.swapaxes(Ploop.reshape(3, len(kfull), Ploop.shape[1]), axis1=1, axis2=2)[:, 1:, :]

    bs = np.array([b1,c2,b3,c4,b5,b6,b7,b8,b9,b10])
    modelX = computePS(bs, Plin, Ploop, kfull, kmin, kmax, withhex = withhex, Puncorr = Puncorr).reshape(-1)

    if masktriangle is None: bisp0 = None
    else:
        b2 = 1./np.sqrt(2.) * (c2 + c4)
        b4 = 1./np.sqrt(2.) * (c2 - c4)
        bval = np.array([1.,b1, b2, b4, 0., b1**2, b1*b2, b1*b4, b1**3, b1**2*b2, b1**2*b4, 0.]) # correction grid
        bisp0 = Bispinterp((lnAs,Om,h))[3:] # removing columns of triangle (k1, k2, k3)
        bisp = np.dot(bval, bisp0)[masktriangle]
        modelX = np.concatenate([ modelX, bisp ])

    if marg_gaussian:
        Pi = get_Pi_for_marg(Ploop, kfull, kmin, kmax, b1, model=model, bisp = bisp0, masktriangle = masktriangle, withhex = withhex)
        Covbi = get_Covbi_for_marg(Pi, Cinv, sigma=sigma_prior)
        Cinvbi = np.linalg.inv(Covbi)
        vectorbi = np.dot(modelX, np.dot(Cinv, Pi.T)) - np.dot(Cinvdata, Pi.T)
        chi2nomar = np.dot(modelX, np.dot(Cinv, modelX)) - 2.*np.dot(Cinvdata, modelX) + chi2data
        chi2mar = -np.dot(vectorbi, np.dot(Cinvbi, vectorbi)) + np.log(np.linalg.det(Covbi))
        chi2 = chi2mar + chi2nomar - 5*np.log(2.*np.pi)
    else:
        chi2 = np.dot(modelX - ydata, np.dot(Cinv, modelX - ydata))

    return -0.5 * chi2

def lnprior(theta, free_para, fix_para, bounds, withPlanck = False):
    value_array = match_para(theta, free_para, fix_para)
    lnAs,Om,h,b1,c2,b3,c4,b5,b6,b7,b8,b9,b10,b11,e1 = value_array

    withinprior = True
    for i in range(len(value_array)): withinprior = (withinprior) and (bounds[i][0] <= value_array[i] <= bounds[i][1])

    if withinprior:
        prior = - 0.5 * (  (b3/2.)**2 + (c4/2.)**2 + (b5/2.)**2 + (b6/4.)**2 + (b7/4.)**2 + (b8/400.)**2 + (b10/1.)**2 )
        if withPlanck: 
            return prior -0.5 * (RD - rs(Om, h, fb))**2 / SIGMA_RD**2
        elif withPlanck2018: 
            return prior + multivariate_normal.logpdf(value_array[:3], centralPlanck2018, covPlanck2018)
        else: return prior
    else: return -np.inf

def lnprob(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, Om_fid, 
    interpolation_grid, sigma_prior = 100, marg_gaussian=False, model = 1, masktriangle = None, withhex = True, Puncorr = 0):
    
    lp = lnprior(theta, free_para, fix_para, bounds, withPlanck = withPlanck)

    if not np.isfinite(lp): 
        dummy = -np.inf
    else:
        dummy = lp + lnlike(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, Om_fid,
                            interpolation_grid, sigma_prior=sigma_prior, 
                            marg_gaussian=marg_gaussian, model = model, masktriangle = masktriangle, withhex = withhex, Puncorr = Puncorr)
    return dummy


###########################################
#  Main program  ########################
###########################################

if __name__ == "__main__":

    ##################################################
    ##### Analysis choices: simulation, box, kmax, partially-marginalized likelihood or full likelihood, model, with bispectrum, ...

    withc4 = False
    withhex = False
    withkbin = False

    withPlanck = bool(int(sys.argv[8]))
    withPlanck2018 = False
    withfibercol = False

    simtype = str(sys.argv[1])
    boxnumber = sys.argv[2]
    kmin = 0.01
    kmax = float(sys.argv[3])
    marg_gaussian = int(sys.argv[4]) # 0 or 1 (True of False)
    model = int(sys.argv[5])
    ZONE = str(sys.argv[7])

    
    if "Challenge" in simtype: 
        ZONE = ''
        nd = 3.e-4
        if 'D' in boxnumber: 
            simname = 'ChallengeD'
            gridname = 'ChallengeD'
            RD = 147.253
        elif 'PT' in simtype: 
            simname = 'PTChallengeTrue'
            gridname = 'PTChallenge'
        else: 
            simname = 'ChallengeA'
            gridname = 'ChallengeA'
            RD = 147.68188

    if 'LightConeHector' in simtype:
        nd = 4.5e-4 # Gil-Marin Table 1
        gridname = 'LightConeHector' + ZONE

        simname = 'LightConeHector'
        if 'data' not in boxnumber: 
            RD = 147.652
            nd = 3.e-4

    if 'lowz' in simtype:
        simname = 'LightConeHector'
        gridname = 'lowz'
        nd = 4.e-4

    kminbisp = 0.04
    kmaxbisp = float(sys.argv[6])
    withBisp = False
    if kmaxbisp > 0: 
        withBisp = True
        withc4 = True

    runtype = simtype+ZONE
    if marg_gaussian: runtype += 'gaussMarg'
    if not withhex:   runtype += 'nohex'
    if withkbin:      runtype += 'withkbin'
    if not withc4:    runtype += 'noc4'

    runtype += 'kmin%s'%kmin

    if withPlanck: runtype += 'withPlanck'
    if withPlanck2018: runtype += 'withPlanck2018'

    if withBisp: runtype += 'withBisp%s'%kmaxbisp

    runtype += '-model%s'%model

    print ('Running %s'%runtype)

    ##################################################
    ##### Loading power spectrum data, covariance and Nkmu binning | loading simulation specification
    kpred, chi2data, Cinvdata, Cinv, TableNkmu, masktriangle, xdata, ydata, yerror = import_data(simtype, boxnumber, kmin, kmax, kminbisp, kmaxbisp, ZONE, withhex = withhex)

    simspec = import_simspec_from_DataFrameCosmosim(simname)
    z_pk, cosmosim, Om_fid, fb = simspec

    ##################################################
    ##### Parameters to minimize ######
    if marg_gaussian:
        free_para =  [True,True,True,
                            True,True,False,withc4,
                            False,False,False,
                            False,False,False,False,
                            False]
        a = 2. # emcee temperature scale factor
    else:
        a = 1.5
        if      model == 2:    free_model = [False, False, True, withBisp, withBisp] # paper model
        elif    model == 1:    free_model = [True, False, True, withBisp, withBisp]  # paper model + shot noise
        else:                  free_model = [False, False, False, withBisp, withBisp]   # paper model without stochastic terms

            free_para = [True, True, True, True, True, True, withc4, True, True, withhex] + free_model

    ndim  =  sum(free_para)
    fix_para  =  np.concatenate(([cosmosim[0], cosmosim[1], cosmosim[2]], [2.] + 11*[0])) # if free_para is false read the value in fix_para

    ##################################################
    ##### Loading Grid ######
    if 'LightConeHector' in simname or 'lowz' in simname:
        if 'small' in gridname:
            nbinsAs = 70
            nbinsOm = 48
            nbinsh = 60 
        else:
            nbinsAs = 70
            nbinsOm = 80
            nbinsh = 120
    elif 'Nseries' in simname:
        nbinsAs = 70
        nbinsOm = 48
        nbinsh = 60       
    else:
        nbinsAs = 90
        nbinsOm = 80
        nbinsh = 80

    lnAsmin, lnAsmax, Ommin, Ommax, hmin, hmax, interpolation_grid = get_grid(gridname, nbinsAs=nbinsAs, nbinsOm=nbinsOm, nbinsh=nbinsh, withBisp = withBisp)
    print("got grid %s!"%gridname)
    print ('bounds: ', lnAsmin, lnAsmax, Ommin, Ommax, hmin, hmax)
    #sys.exit()

    ##################################################
    ##### Uncorrelated part of fiber-collision correction ######
    if withfibercol is True:
        if withBisp: Plininterp, _, _ = interpolation_grid
        else: Plininterp, _ = interpolation_grid
        Plin = Plininterp((cosmosim[0], cosmosim[1], cosmosim[2]))
        kfull = Plin[:, 0]
        kfull = kfull[:int(len(kfull) / 3)]
        Puncorr = fc.dPuncorr(kfull)
    else:
        Puncorr = 0

    ##################################################
    ##### Uniform prior on the As, h, Om and b_i ######
    priorsup = 10.
    priorgauss = 4.
    bfmintab = np.concatenate(([lnAsmin, Ommin, hmin], [0., -4., -priorsup, -priorsup] + 8*[-priorsup])) # We require b_1 > 0
    bfmaxtab = np.concatenate(([lnAsmax, Ommax, hmax], [4., 4., priorsup, priorsup] + 8*[priorsup]))
    bounds = zip(bfmintab, bfmaxtab)

    onesigma = np.array([ 0.3, 0.04, 0.03, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1., 0.5]) # Guess for the \sigma, to help with the initial position of walkers
    tryini = np.array([3.12,  0.307,  0.669, 1.97,  1.56, -4.07, 5.02, 1.33, -6.39, -0.64, -0.21, 0, 0, 0, 0]) # initial guess

    ##################################################
    ##### find max likelihood ######
    try:
        tryini2 = [i for (i, v) in zip(tryini, free_para) if v]
        boundschi2 = [i for (i, v) in zip(bounds, free_para) if v]
        def chi2(theta): return -2. * lnlike(theta, chi2data, Cinvdata, Cinv, kmin, kmax, free_para, fix_para, bounds, Om_fid, interpolation_grid, 
           sigma_prior=priorgauss, marg_gaussian=marg_gaussian, model=model, masktriangle = masktriangle, withhex = withhex)

        result = op.minimize(chi2, tryini2, method='SLSQP', bounds=boundschi2, options={'maxiter': 1000})

        bestfit =  match_para(result["x"], free_para, fix_para)
        minchi2 = result["fun"]
        dof = len(ydata)#s-sum(free_para)
        pvalue = 1. - stats.chi2.cdf(minchi2, dof)

        print (bestfit)
        print( 'minchi2 /d.o.f. = %s / %s, p-value fit: %s' % ( minchi2, dof,  pvalue) )
        np.savetxt(os.path.join(OUTPATH, "minchi2%sbox_%skmax_%s.txt") %
                  (runtype, boxnumber, kmax), np.concatenate([bestfit, [minchi2, dof, pvalue]]))

        lnAs,Om,h,b1,c2,b3,c4,b5,b6,b7,b8,b9,b10,b11,e1 = bestfit
        if not withBisp: Plininterp, Ploopinterp = interpolation_grid
        else:   Plininterp,Ploopinterp,Bispinterp= interpolation_grid 
        Plin = Plininterp((lnAs, Om, h))
        Ploop = Ploopinterp((lnAs, Om, h))
        kfull = Plin[:, 0]
        if check_if_multipoles_k_array(kfull): kfull = kfull[:int(len(kfull) / 3)]
        Plin = np.swapaxes(Plin.reshape(3, len(kfull), Plin.shape[1]), axis1=1, axis2=2)[:, 1:, :]
        Ploop = np.swapaxes(Ploop.reshape(3, len(kfull), Ploop.shape[1]), axis1=1, axis2=2)[:, 1:, :]
        bs = np.array([b1,c2,b3,c4,b5,b6,b7,b8,b9,b10])
        bestfitPS = computePS(bs, Plin, Ploop, kfull, kmin, kmax, withhex = withhex).reshape(-1)
        if withBisp:
            b2 = 1./np.sqrt(2.) * (c2 + c4)
            b4 = 1./np.sqrt(2.) * (c2 - c4)
            bval = np.array([1.,b1, b2, b4, 0., b1**2, b1*b2, b1*b4, b1**3, b1**2*b2, b1**2*b4, 0.])
            bisp0 = Bispinterp((lnAs,Om,h))[3:] # removing columns of triangle (k1, k2, k3)
            bisp = np.dot(bval, bisp0)[masktriangle] #+ 100. * e1  / nd**2 * np.ones(sum(masktrianglegrid)) # applying masktriangle (triangle < kmaxbisp)
            bestfitPS = np.concatenate((bestfitPS, bisp), axis = 0)

            print (xdata.shape, ydata.shape, yerror.shape, bestfitPS.shape, np.array([ xdata, ydata, yerror, bestfitPS ]).shape)

        np.savetxt( opa.join(CHAINPATH, "minchi2_powerspectrum%sbox_%skmax_%s.txt") % (runtype,boxnumber,kmax) , np.array([ xdata, ydata, yerror, bestfitPS ]) )

    except:
        print ('min chi^2 not found!')

    #sys.exit(0)
    bestfit = tryini

    ###################################
    # run MCMC #######################
    ###################################

    # Start MCMC
    t0 = time.time()
    Nchains  =  4
    nwalkers  =  16*ndim
    nthreads = cpu_count()/Nchains

    # Set up the sampler.
    pos = []
    sampler = []
    rstate = np.random.get_state()

    def lnprobloaded(theta):
        return lnprob(theta, chi2data, Cinvdata, Cinv, kmin, kmax,
                      free_para, fix_para, bounds, Om_fid,
                      interpolation_grid,
                      sigma_prior=priorgauss, marg_gaussian=marg_gaussian, 
                      model=model, masktriangle = masktriangle, withhex = withhex, Puncorr = Puncorr)

    t_try = time.time()
    for jj in range(0, Nchains):
        initialpos = []
        for ii in xrange(nwalkers):
            accepted  =  False
            t_try = time.time()
            while (not accepted):
                trialfiducial  =  np.random.normal(loc = bestfit[free_para], scale =  onesigma[free_para])
                accepted  =  np.isfinite(lnprior(trialfiducial, free_para, fix_para,bounds))
            if accepted: initialpos.append(trialfiducial)
        pos.append(initialpos)
        sampler.append(emcee.EnsembleSampler(nwalkers, ndim, lnprobloaded, a = a, threads = nthreads))
    print("Found initial pos in %s seconds!" % str(time.time() - t_try))

    ##################################################
    epsilon = 0.003
    
    withinchainvar  =  np.zeros( (Nchains, ndim) )
    meanchain  =  np.zeros( (Nchains, ndim) )
    scalereduction  =  np.arange(ndim, dtype = np.float)
    for jj in range(0, ndim): scalereduction[jj]  =  2.
    minlength  =  6000
    chainstep  =  minlength
    itercounter  =  0
    loopcriteria  =  1
    ichaincheck  =  50
    ithin  =  1
    t1 = time.time()

    while loopcriteria:
        itercounter  =  itercounter + chainstep
        print("chain length  = ", itercounter, " minlength  = ", minlength)
        samplesJG = []
        for jj in range(0, Nchains):
            c = 0
            for result in sampler[jj].sample(pos[jj], iterations = chainstep, rstate0 = rstate, storechain = True, thin = ithin):
                pos[jj]  =  result[0]
                chainchi2  =  -2.*result[1]
                rstate  =  result[2]
            # we do the convergence test on the second half of the current chain (itercounter/2)
            chainsamples  =  sampler[jj].chain[:, itercounter/2:, :].reshape((-1, ndim))
            withinchainvar[jj]  =  np.var(chainsamples, axis = 0)
            meanchain[jj]  =  np.mean(chainsamples, axis = 0)
            samplesJG.append(chainsamples)
        scalereduction  =  gelman_rubin_convergence(withinchainvar, meanchain, itercounter/2, Nchains, ndim)
        print("scalereduction  =  ", scalereduction)
        loopcriteria  =  0
        for jj in range(0, ndim):
            if np.absolute(1-scalereduction[jj]) > epsilon: loopcriteria  =  1

        chainstep  =  ichaincheck

    print("Done.")

    ##################################################
    ##### Saving results #######################
    samplerchain = sampler[0].chain.reshape((-1, ndim))
    lnlike = sampler[0].lnprobability.reshape(-1)
    np.save(os.path.join(CHAINPATH, "samplerchain%sbox_%skmax_%s") % (runtype, boxnumber, kmax), samplerchain)
    np.save(os.path.join(CHAINPATH, "lnlikechain%sbox_%skmax_%s") % (runtype, boxnumber, kmax), lnlike)

    bestfit =  match_para(samplerchain[lnlike.argmax(),:], free_para, fix_para)
    minchi2 = -2.*lnlike.max()
    dof = len(xdata)-sum(free_para)
    pvalue = 1. - stats.chi2.cdf(minchi2, dof)
    np.savetxt( os.path.join(CHAINPATH, "bestfit-mcmc%sbox_%skmax_%s.txt") % (runtype, boxnumber, kmax), np.concatenate([ bestfit, [minchi2, dof, pvalue] ]), fmt = '%.3f' )

    # Compute the quantiles ; columns: expectation values, +sigma, -sigma, meansigma
    mcmc_array = np.array(list(map(lambda v: ( v[1], v[2] - v[1], v[1] - v[0], 0.5*(v[2]-v[0]) ),
        zip(*np.percentile(samplerchain, [16, 50, 84], axis=0)))))

    ncosmo = 3
    theoryerror = np.zeros(ncosmo)
    meancosmo = mcmc_array[:ncosmo,:]
    for i, (t, m) in enumerate(zip(cosmosim, meancosmo)):
        if t > m[0]+m[1]: theoryerror[i] = t - (m[0]+m[1])
        elif t < m[0]-m[2]: theoryerror[i] = t - (m[0]-m[2])
    results = np.c_[cosmosim[:ncosmo], bestfit[:ncosmo], meancosmo, theoryerror ]
    print (results)

    np.savetxt( os.path.join(CHAINPATH, "mcmcresults%sbox_%skmax_%s.txt") % (runtype, boxnumber, kmax), results, fmt = '%.3f' )

    # Saving bestfit power spectra
    lnAs,Om,h,b1,c2,b3,c4,b5,b6,b7,b8,b9,b10,b11,e1 = bestfit
    if not withBisp: Plininterp, Ploopinterp = interpolation_grid
    else:   Plininterp,Ploopinterp,Bispinterp= interpolation_grid 
    Plin = Plininterp((lnAs, Om, h))
    Ploop = Ploopinterp((lnAs, Om, h))
    kfull = Plin[:, 0]
    if check_if_multipoles_k_array(kfull): kfull = kfull[:int(len(kfull) / 3)]
    Plin = np.swapaxes(Plin.reshape(3, len(kfull), Plin.shape[1]), axis1=1, axis2=2)[:, 1:, :]
    Ploop = np.swapaxes(Ploop.reshape(3, len(kfull), Ploop.shape[1]), axis1=1, axis2=2)[:, 1:, :]
    bs = np.array([b1,c2,b3,c4,b5,b6,b7,b8,b9,b10])
    bestfitPS = computePS(bs, Plin, Ploop, kfull, kmin, kmax, withhex = withhex)
    if withBisp:
        b2 = 1./np.sqrt(2.) * (c2 + c4)
        b4 = 1./np.sqrt(2.) * (c2 - c4)
        bval = np.array([1.,b1, b2, b4, 0., b1**2, b1*b2, b1*b4, b1**3, b1**2*b2, b1**2*b4, 0.])
        bisp0 = Bispinterp((lnAs,Om,h))[3:] # removing columns of triangle (k1, k2, k3)
        bisp = np.dot(bval, bisp0)[masktrianglegrid] + 100. * e1  / nd**2 * np.ones(sum(masktrianglegrid)) # applying masktriangle (triangle < kmaxbisp)
        bestfitPS = np.concatenate([ bestfitPS, bisp ])
    
    np.save( opa.join(CHAINPATH, "powerspectrum%sbox_%skmax_%s.npy") % (runtype,boxnumber,kmax) , [ xdata, ydata, yerror, bestfitPS ] )

    # Estimate the integrated autocorrelation time for the time series in each parameter.
    print("Autocorrelation time:", sampler[0].get_autocorr_time(c=1))

    # Print out the mean acceptance fraction. In general, acceptance_fraction
    # has an entry for each walker so, in this case, it is a 250-dimensional vector.
    print("Mean acceptance fraction: ", np.mean(sampler[0].acceptance_fraction))
    np.savetxt(os.path.join(CHAINPATH, 'AcceptanceFr%sbox_%skmax_%s.dat' % (runtype, boxnumber, kmax)), sampler[0].acceptance_fraction, fmt = '%.2f' )
