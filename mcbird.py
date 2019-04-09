#!/usr/bin/env python

from __future__ import print_function

import os
import mpi4py
mpi4py.rc.initialize = False
mpi4py.rc.finalize = False
import emcee
import numpy as np
import scipy.stats
import scipy.optimize as op
from scipy.interpolate import interp1d
import os.path as opa
import time
from scipy import stats
import sys
#print("imported")
sys.stdout.flush()

THIS_PATH  =  opa.dirname(__file__)

import taylorbird as tbird
import tayloregg as tegg

INPATH = opa.abspath('/exports/pierre/EFTofBOSS/input/')
OUTPATH = opa.abspath(opa.join(THIS_PATH,'mcnest/')) 
if not opa.isdir(OUTPATH): raise Exception(OUTPATH + ' not there!')

nd = 3e-4
km = 0.7
knl = 0.7

def embed_Pi(ploop, masktriangle):
    #Embed the Pi into a bigger shape that has appropriate prepadding and postpadding so the bisp term can be easily added on
    #Dimension of Pi is [nterms, 3, 100]
    #return with dimension [nterms+1, 4, 100+nkbisp]
    nkbisp = sum(masktriangle)
    nkp = ploop.shape[1]
    #Cast into bigger array of dim p
    big_array = np.zeros(shape=(ploop.shape[0]+1,100+nkbisp)) 
    return big_array

def get_Pi_for_marg(Ploop, kpred, b1, kfull, bisp=None):
    nk = len(kfull)
    Onel0 = np.array([np.ones(nk),np.zeros(nk),np.zeros(nk)])
    kl0 = np.array([kfull,np.zeros(nk),np.zeros(nk)])
    kl2 = np.array([np.zeros(nk),kfull,np.zeros(nk)])

    Pi = np.array([Ploop[:,3,:]+b1*Ploop[:,7,:],
             (Ploop[:,15,:]+b1*Ploop[:,12,:]) / knl**2,
             (Ploop[:,16,:]+b1*Ploop[:,13,:]) / km**2,
             (Ploop[:,17,:]+b1*Ploop[:,14,:]) / km**2,
             Onel0 / nd,
             (kl0**2) / nd / km**2,
             kl2**2 / nd / km**2])

    if withBisp:
        Pi = np.array([Ploop[:,3,:]+b1*Ploop[:,7,:],
                 (Ploop[:,15,:]+b1*Ploop[:,12,:]) / knl**2,
                 (Ploop[:,16,:]+b1*Ploop[:,13,:]) / km**2,
                 (Ploop[:,17,:]+b1*Ploop[:,14,:]) / km**2])
                 #(kl0**2 + kl2**2) / nd / km**2,
                 #kl2**2 / nd / km**2])
    return interp1d(kfull, Pi)(kpred)
    

def get_Covbi_for_marg(Pi_data,Cinv,sigma=200):
    Covbi = np.dot(Pi_data,np.dot(Cinv,Pi_data.T))+ 1./sigma**2*np.identity(Pi_data.shape[0])
    return Covbi

def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0]    
            
def match_para(theta, free_para, fix_para):    
    value_array  =  np.empty(len(free_para),dtype = np.float)
    counter  =  0
    for i in range(len(free_para)):
        if(free_para[i]  ==  True):
            value_array[i] = theta[counter]
            counter +=  1
        else: value_array[i] = fix_para[i]
    #print (value_array)
    return value_array

def lnlike(theta,  kpred, chi2data, Cinvwdata, Cinvww, free_para, fix_para, bounds, cosmoref, tayloregg, sigma_prior = 100, marg_gaussian=False):
    if marg_gaussian:     # If marginalization, no variation over b3,b5,b6,b7,b8,b9,b10,b11. 
        #print(free_para)
        free_para = np.array(free_para)            
        free_para[[5,7,8,9,10,11,12,13,14,15]] = np.array([False]*10)
        fix_para[[5,7,8,9,10,11,12,13,14,15]] = np.array([0]*10)
        if withBisp:       #b8 is not marginalized over when bisp is present
            free_para[[5,7,8,9,11,12,13,14,15]] = np.array([False]*9)
            fix_para[[5,7,8,9,11,12,13,14,15]] = np.array([0]*9)
    #print(match_para(theta, free_para, fix_para))
    lnAs,Om,h,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,e1,e2 = match_para(theta, free_para, fix_para)
    
    '''
    if withBisp:
        if type(masktriangle) ==  type(None) or type(Bispdata) == type(None) or type(Bispinterp) ==  type(None):
            raise Exception('You want to use the bispectrum but forgot to provide a mask for the triangle or the data or the interpolation of the Bisp. Can be found in input/')
        if Cinv.shape[0] != xdata.shape + sum(masktriangle):
            raise Exception('You want to use the bispectrum but forgot to use the full covariance for power spectrum + Bisp. Can be found in input/Covariance')
    
        TermsBisp = Bispinterp((lnAs,Om,h))
        c2 = 0.5 * (b2 + b4) 
        c4 = 0.5 * (b2 - b4)
        bval = np.array([1.,b1,c2,c4,b1*b11/nd * 0.00952,b1**2,b1*c2,b1*c4,b1**3,b1**2*c2,b1**2*c4,b8**2/nd**2 * 0.00952**2]) # correction grid
        Bisp = 1.*np.dot(bval,TermsBisp[3:])
    '''
    bs = np.array([b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,e1,e2])
    cosmotarget = np.array([lnAs, Om, h, cosmoref[3], cosmoref[4]])  
    k, Plin, Ploop = tbird.taylorbird(cosmoref, cosmotarget, tayloregg)
    modelX = tbird.computePS(bs, Plin, Ploop, k, kpred).reshape(-1)

    '''
    if withBisp:
        modelX = np.concatenate([modelX, Bisp[masktriangle2]])
    '''
    
    if marg_gaussian:

        if withBisp:
            pass
            '''
            Pi = get_Pi_for_marg(Ploop, kpredb1, TermsBisp[3:])
            #Building a bigger Pi to include bispectrum term
            nparams = Pi_tot.shape[0]
            nkpred = Pi_tot.shape[1]
            nkbisp = sum(masktriangle)
            #print(Pi_tot.shape, TermsBisp[3:].shape, TermsBisp.shape)
            #Removing triangle contributions from the bispectrum expressions
            TermsBisp = TermsBisp[3:]
            #Applying mask and only getting b11 contribution
            bisp = TermsBisp[4][masktriangle2]
            
            newPi = np.zeros(shape=(nparams+1, nkpred+nkbisp))
            newPi[:nparams, :nkpred] = Pi_tot
            newPi[-1, nkpred:] = b1*bisp
            #total Pi is now the correctly embedded one
            Pi_tot = 1.*newPi 
            '''
        else: 
            Pi = get_Pi_for_marg(Ploop, kpred, b1, k)
            Pi_tot = 1.*Pi.reshape((Pi.shape[0],-1))

        Covbi = get_Covbi_for_marg(Pi_tot, Cinvww, sigma=sigma_prior)
        Cinvbi = np.linalg.inv(Covbi)
        vectorbi = np.dot( modelX, np.dot(Cinvww,Pi_tot.T) ) - np.dot(Cinvwdata, Pi_tot.T)
        chi2nomar = np.dot( modelX,np.dot(Cinvww,modelX) ) - 2*np.dot(Cinvwdata, modelX) + chi2data
        chi2mar = -np.dot(vectorbi, np.dot(Cinvbi,vectorbi)) + np.log(np.linalg.det(Covbi))
        chi2 = chi2mar + chi2nomar
    else:
         chi2 = np.dot( modelX, np.dot(Cinvww, modelX) ) - 2*np.dot(Cinvwdata, modelX) + chi2data       

    return -0.5*chi2

def lnprior(theta, free_para, fix_para, bounds):    
    value_array  =  match_para(theta, free_para, fix_para)
    withinprior = True
    for i in range(len(value_array)):
        withinprior = (withinprior) and (bounds[i][0] <= value_array[i] <= bounds[i][1]) 
    if withinprior:
        return 0.
    else:
        return -np.inf

def lnprob(theta, kpred, chi2data, Cinvwdata, Cinvww, free_para, fix_para, bounds, cosmoref, tayloregg, sigma_prior = 100, marg_gaussian=False):
    lp  =  lnprior(theta, free_para, fix_para, bounds)
    if np.isfinite(lp) == False :
        dummy  =  -np.inf
    dummy  =  lp + lnlike(theta, kpred, chi2data, Cinvwdata, Cinvww, free_para, fix_para, bounds, cosmoref, tayloregg, sigma_prior=sigma_prior, marg_gaussian=marg_gaussian)
    return dummy


    ############################################################################################################
    ###  Main program  ########################
    ############################################################################################################

if __name__ ==  "__main__":

    #print("started")
    boxnumber = sys.argv[1]
    kmin = 0.01
    kmax = float(sys.argv[2])
    simtype = sys.argv[3]
    marg_gaussian = int(sys.argv[4]) # 0 or 1 (True of False)
    stoch = float(sys.argv[5])
    ZONE = 'NGC'
    if "Challenge" in simtype:
        ZONE = ''
    runtype = simtype+ZONE
    if marg_gaussian:
        runtype += 'gaussMarg'

    withBisp = False 
    kminbisp = 0.04
    kmaxbisp = float(sys.argv[6])
    if kmaxbisp > 0:
        withBisp = True
        print('kmax bispectrum is bigger than zero, withBisp is %s'%withBisp)

    ##################################################
    ##### Loading covariance ##############
    if 'Challenge' in simtype:
        if 'D' in simtype:
            simname = 'ChallengeD'
        else:
            simname = 'ChallengeA'
        if 'Quarter' not in simtype:
            simtype_false = 'ChallengeQuarter%s'%boxnumber
            print('Using quarter covariance divided by 4.25 instead of full')
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simtype_false,ZONE))) / 4.25
        else:
            simlettered = 'ChallengeQuarter%s'%boxnumber
            Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s.dat'%(simlettered,ZONE)))
    else:
        simname = simtype
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%sdata.dat'%(simtype,ZONE)))

    Bispdata = [] 
    masktriangle = []
    if withBisp:
        runtype += 'withBispkmax%s'%kmaxbisp
        Full_Cov = np.loadtxt(opa.join(INPATH,'Covariance/Cov%s%s_Bisp.dat'%(simtype,ZONE)))   
        Q1,Q2,Q3,Bispdata = np.loadtxt(opa.join(INPATH,'DataSims/Bispred_LightConeHector_%s_%s.dat'%(ZONE,boxnumber))).T
        KMAXBISP = 0.15
        R1,R2,R3 = Q1[(Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)], Q2[(Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)], Q3[(Q1<KMAXBISP)&(Q2<KMAXBISP)&(Q3<KMAXBISP)]
        masktriangle = (Q1>=kminbisp)&(Q1<=kmaxbisp)&(Q1<=Q2+Q3)&(Q1>=abs(Q2-Q3))&(Q2>=kminbisp)&(Q2<=kmaxbisp)&(Q3>=kminbisp)&(Q3<=kmaxbisp)
        masktriangle2 = (R1>=kminbisp)&(R1<=kmaxbisp)&(R1<=R2+R3)&(R1>=abs(R2-R3))&(R2>=kminbisp)&(R2<=kmaxbisp)&(R3>=kminbisp)&(R3<=kmaxbisp)

    ##################################################
    ##### Loading power spectrum data ####
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

    chi2data = np.dot(ydata, np.dot(Cinv, ydata))
    Cinvwdata = np.dot(ydata, Cinvw)  

    ##################################################
    ##### Loading Taylor-egg ##############
    z_pk, cosmoref, cosmodelta = tegg.import_cosmoref_from_DataFrameCosmosim(simname)
    #tayloregg = tbird.load_taylor('%s-shifted'%simname)
    tayloregg = tbird.load_taylor('%s'%simname)


    ##################################################
    ##### Parameters to minimize ######
    if marg_gaussian:
        #b8 is free if bispectrum, b11 is never varied now (marged if bisp)
        free_para =  [True,True,True,True,True,False,True,False,False,False,withBisp,False,False,False,False,False]
        a = 1.8
        print("With marginalization")
    else:
        a = 1.5
        print("No marginalization")
        #free_para =  [True,True,True,True,True,True,True,True,True,True,True,False,False, withBisp,True,True]
        
        if stoch == 1:
            free_stoch = [False, True, withBisp, False, False]
        elif stoch == 2:
            free_stoch = [False, True, withBisp, False, True]
        elif stoch == 3:
            free_stoch = [False, True, withBisp, True, True]
        elif stoch ==4:
            free_stoch = [True, True, withBisp, True, True]
        else:
            free_stoch = [False, False, withBisp, False, False]

        free_para =  11*[True] + free_stoch
        runtype += '-shiftedstoch%s'%stoch

    all_true  =  np.concatenate(([cosmoref[0], cosmoref[1], cosmoref[2]], [2.] + 12*[0]))
    all_name  =  np.concatenate(([r'$A_s$',r'$\Omega_m$',r'$h$'],[r'$b_%s$'%(i+1) for i in range(13)]))

    ndim  =  sum(free_para)
    fix_para  =  all_true # if free_para is false read the value in fix_para
    free_true = all_true[free_para]
    free_name = all_name[free_para]

    ##################################################
    ##### Uniform prior on the As, h, Om and b_i ######
    priorsup = 20
    bmin = -priorsup
    bmax = priorsup
    b4sup = 100
    # We require b_1>0 and b_4 is large due to the PCA that we do between b2 and b4
    bmintab = [0, bmin, bmin, -b4sup, bmin, bmin, bmin, bmin, bmin, bmin, bmin, bmin, bmin] 
    bmaxtab = [bmax, bmax, bmax, b4sup, bmax, bmax, bmax, bmax, bmax, bmax, bmax, bmax, bmax]
    bfmintab = np.concatenate([[2., 0.1, 0.3], bmintab])
    bfmaxtab = np.concatenate([[4., 0.5, 1.], bmaxtab])
    bounds = zip(bfmintab, bfmaxtab)

    ##################################################
    ##### Find maximum likelihood ######
    chi2  =  lambda theta: -2 * lnlike(theta, kpred, chi2data, Cinvwdata, Cinvww, 
        free_para, fix_para, bounds, cosmoref, tayloregg, sigma_prior = priorsup, marg_gaussian=True) ### Keep marg_gaussian = True to be sure to converge

    result  =  op.minimize(chi2, all_true, method = 'SLSQP',bounds = bounds, options = {'maxiter':1000})
    all_ml  =  result["x"] 
    free_ml = all_ml[free_para]
    minchi2  =  result["fun"]
    #print(result)

    if type(masktriangle) == type(None):
        dof = len(xdata) - ndim
    else:
        dof = len(xdata) + sum(masktriangle) - ndim
    #print( 'marg: minchi2 /d.o.f. = %s / %s' % (minchi2, dof) )
    #np.savetxt(opa.join(OUTPATH,"minchi2%sbox_%skmax_%s.txt")%(runtype,boxnumber,kmax),np.concatenate([free_ml,[minchi2,dof]]))

    
    ##################################################
    ##### Find real maximum likelihood ######
    params_i = np.array([3.12,  0.307,  0.669, 1.97,  1.56, -4.07, 5.02, 1.33, -6.39, -0.64, -0.21, 0, 0, 0, 0, 0])
    #para_nomarg =  [True,True,True,True,True,True,True,True,True,True,True,True,True, withBisp,True,True]
    para_nomarg =  free_para#[True,True,True,True,True,True,True,True,True,True,True,False,True, withBisp,False,False]

    chi2_real  =  lambda theta: -2 * lnlike(theta, kpred, chi2data, Cinvwdata, Cinvww, 
        para_nomarg, params_i, bounds, cosmoref, tayloregg, sigma_prior = priorsup, marg_gaussian=False) ### Keep marg_gaussian = True to be sure to converge

    result_real  =  op.minimize(chi2_real, params_i, method = 'SLSQP',bounds = bounds, options = {'maxiter':10000})
    #print(result_real)
    print (result_real["x"])
    minchi2_real = result_real["fun"]
    dof_real = len(xdata)-sum(para_nomarg)
    print( 'nonmarg: minchi2 /d.o.f. = %s / %s, p-value fit: %s' % ( minchi2_real, dof_real, 1. - stats.chi2.cdf(minchi2_real, dof_real) ) )
    
    #sys.exit()
    ### FORCED VALUES
    free_ml = params_i[free_para]


    ##################################################
    ##### run MCMC #######################
    ichaincheck = 50
    ithin = 1
    itercounter = 0
    nwalkers = 300 * ndim
    chainstep = 4000
    burnin = 1000
    nthreads = 24

    # Set up the sampler.
    pos = []
    sampler = []
    rstate = np.random.get_state()

    def lnprobloaded(theta):
        return lnprob(theta, kpred, chi2data, Cinvwdata, Cinvww, free_para, fix_para, bounds, cosmoref, tayloregg, sigma_prior=priorsup, marg_gaussian=marg_gaussian)

    initialpos = []
    onesigma = np.array([ 0.25, 0.04, 0.02, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 1, 0.1, 0.1]) 
    t_try = time.time()
    for ii in xrange(nwalkers):
        accepted = False
        while (not accepted):
            trialfiducial = np.random.normal(loc=free_ml, scale=onesigma[free_para])
            accepted = np.isfinite(lnprior(trialfiducial, free_para, fix_para, bounds))
        if accepted:
            initialpos.append(trialfiducial)
    np.save(os.path.join(OUTPATH, "inipos%sbox_%skmax_%s") % (runtype, boxnumber, kmax), np.array(initialpos))
    print("Found initial pos in %s seconds!" % str(time.time() - t_try))

    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprobloaded, a=a, threads=nthreads)

    # Start MCMC
    print("Running MCMC...")

    # Do burn-in
    t_burnin = time.time()
    pos, prob, state = sampler.run_mcmc(initialpos, burnin, storechain=False)
    sampler.reset()
    print("Done burn-in step in %s seconds!" % str(time.time() - t_burnin))

    for i, result in enumerate(sampler.sample(pos, iterations=chainstep, storechain=True, thin=ithin)):
        if (i + 1) % 1000 == 0:
            pass
            #np.save(os.path.join(OUTPATH, "ChainsMidway", "samplerchainmid%sbox_%skmax_%s") %
            #        (runtype, boxnumber, kmax), sampler.chain)
            #np.save(os.path.join(OUTPATH, "ChainsMidway", "lnlikechainmid%sbox_%skmax_%s") %
            #        (runtype, boxnumber, kmax), sampler.lnprobability)

    # Print out the mean acceptance fraction. In general, acceptance_fraction
    # has an entry for each walker so, in this case, it is a 250-dimensional vector.
    print("Mean acceptance fraction: ", np.mean(sampler.acceptance_fraction))

    # Estimate the integrated autocorrelation time for the time series in each parameter.
    # print("Autocorrelation time:", sampler.get_autocorr_time())

    np.save(os.path.join(OUTPATH, "samplerchain%sbox_%skmax_%s") % (runtype, boxnumber, kmax), sampler.chain)
    np.save(os.path.join(OUTPATH, "lnlikechain%sbox_%skmax_%s") % (runtype, boxnumber, kmax), sampler.lnprobability)

    ##################################################
    ##### Compute the quantiles ####
    mcmc_array = list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0], v[4] - v[1], v[1] - v[3]),
                          zip(*np.percentile(sampler.chain.reshape((-1, ndim)),
                                             [15.86555, 50, 84.13445, 2.2775, 97.7225], axis=0))))

    np.savetxt(os.path.join(OUTPATH, "mcmcarray%sbox_%skmax_%s.txt") % (runtype, boxnumber, kmax), np.array(mcmc_array))

    print("Autocorrelation time:", sampler.get_autocorr_time(c=1))
    np.savetxt(os.path.join(OUTPATH, 'AcceptanceFr%sbox_%skmax_%s.dat' % (runtype, boxnumber, kmax)), sampler.acceptance_fraction)