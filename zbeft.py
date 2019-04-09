import numpy as np
import os.path as op
import zbetools, scipy.interpolate
import time

t0 = time.time()

BNMS = ['b1','b2','b3','b4', # bias parameters
        'b5','b6','b7', # counter-term parameters
        'b8','b9','b10']  # stochastic-term parameters

def coeffs(bnms,bvals,cnms,dmonly=False,kvals=None,kren=None,maxb=np.inf,verb=False):
    
    # maxb must be an integer!
    if type(maxb) is not int and maxb != np.inf: 
        raise TypeError('maxb must be an integer corresponding to the'+\
                        ' index of the last b in bnms that should be used!')

    # Multiply into the k-dependent coeffs in the DM-only case
    if dmonly and kvals is not None and kren is not None:
        bvals[BNMS.index('b5')] = bvals[BNMS.index('b5')]*np.log(kvals/kren)
        bvals[BNMS.index('b6')] = bvals[BNMS.index('b6')]*np.log(kvals/kren)
        bvals[BNMS.index('b7')] = bvals[BNMS.index('b7')]*np.log(kvals/kren)
    elif dmonly:
        raise Exception('Need to provide kvals and kren if dmonly!')

    # Set all coefficients to 1 to begin with
    cvals = [1.]*len(cnms)

    # Loop over coefficients and calculate the values by multiplication of bs
    for i,c in enumerate(cnms):

        # Multiply all the factors of the coefficient together
        bs = c.split('*')
        for b in bs:

            # Stop omit the term if it contains bias beyond the max desired
            if bnms.index(b)>maxb: cvals[i]=0.; break
            # This is not enough for lin, the inputs must also be for the linear, 
            # but should work for DM even with EFT inputs

            try:
                cvals[i] *= bvals[bnms.index(b)]
            except Exception as e:
                print bvals
                print bnms
                print b
                raise e
            if verb: print 'coeff*'+b+'=',bvals[bnms.index(b)]

    return cvals

def sum_zbEFT(pkterms, calphas, pklin=None, verb=False):
            
    # Pre-add the linear terms if they are given
    if pklin is None: 
        pkvals = pkterms[:,0]*0.
    else:
        pkvals = pklin.copy()
    
    # Sum over the EFT terms: SPT + counter-terms + stochastic terms
    for i,c in enumerate(calphas):

        if i==pkterms.shape[1]: break

        # Multiply each k-dependent term with the k-independent coefficient and add to the total
        pkvals += c*pkterms[:,i]

    if verb: print 'P_gg^l(k) = ', pkvals[-1], pkvals.shape

    return pkvals

def plks_zbEFT_fromterms(kvals, data, cnames, bvals, bnames, 
               maxb=np.inf, dmonly=False, kren=None, lintoo=False, lincols=3):
    # Multipy k-dependent terms with the coefficients and add them together
    
    # Calculate the values of the coefficients
    cas = coeffs(['1']+bnames,[1.]+list(bvals),cnames,dmonly,kvals,kren,maxb,verb=False)

    # Save summed power spectra into the final array
    pks = sum_zbEFT(data,cas,verb=False)

    if not lintoo:
        return kvals, pks
    else:
        return kvals, pks, sum_zbEFT(data[:,:lincols],cas[:lincols],verb=False)

def plks_zbEFT_fromfile(bvals, basepath='./'+zbetools.FBSNM, maxb=np.inf,
                        dmonly=False, kren=None, lintoo=False, 
                        kvals=None, mlps=[0,2,4],
                        eft=zbetools.FEFT, lin=zbetools.FLIN, bnms=BNMS, lincols=3):
    
    if dmonly and kren is None: raise Exception('If dmonly, need kren input!')

    # Get power spectra terms from the file
    mpls, kvals_f, data_f, names = zbetools.read_zbEFT(basepath, mlps=mlps, eft=eft, lin=lin)

    # Interpolate if neccessary (first double check the input )
    if kvals is not None and (len(kvals)!=len(kvals_f) or not np.any(kvals,kvals_f)):

        # Store the length of the available and the desired k arrays
        nl_f = zbetools.count_mlps(kvals_f)
        nk_f = len(kvals_f)/nl_f 
        nl = zbetools.count_mlps(kvals)
        nk = len(kvals)/nl 

        # Make sure the requested number of multipoles is self-consistent
        if (nl!=len(mlps) and nl!=1) or nl_f<nl:
            raise Exception('The requested number of multipoles is inconsistent! -> ' + str([nl,len(mlps),nl_f]))
        elif nl==1:
            kvals = np.concatenate([kvals]*len(mlps))
            nl = len(mlps)

        # Interpolate each multipole separately, because otherwise they will be mixed up/broken
        data = np.zeros((nk*nl,len(names)))
        for l in range(len(mlps)):

            # Save some stuff for easy reading
            imin = nk_f*l
            imax = nk_f*(l+1)
            k = kvals[nk*l:nk*(l+1)]

            # Now create interpolator and immediately query values and save to output data array
            try:
                data[nk*l:nk*(l+1),:] = scipy.interpolate.interp1d(kvals_f[imin:imax],data_f[imin:imax,:],axis=0)(k)
            except ValueError as e:
                message='\n You are asking for a range of ['+str(np.min(k))+', '+str(np.max(k))
                message+='], but we only have a range of ['+str(np.min(kvals_f[imin:imax]))+', '+str(np.max(kvals_f[imin:imax]))
                message+='] to interpolate.'
                raise type(e)(message)

    else:
        kvals = kvals_f
        data = data_f
        del kvals_f, data_f

    # Get the summed multipoles
    return plks_zbEFT_fromterms(kvals, data, names, bvals, bnms, maxb, dmonly, kren, lintoo, lincols)

def run_zbEFT(config):
    """ Runs Pierre Zhang's EFT code in z-space with biased tracers (written in collab
        with Leonardo Senatore and Ashley Perko.

        Inputs
        ------
        config : str or configobj.ConfigObj or dict
                either a string path to the full config file
                or a configobj.ConfigObj or a dictionary containing
                    the parameter names (predetermined - see zbetools.py) and values

        Outputs
        -------
        pathtoresult : string
            local path to the resulting files
    """

    # Get the individual configurations from the one big common config file/object/dict
    if type(config) == str:
        cclass,czbEFT,czbEFTw = zbetools.get_config(bigconfig_file=config)
    else:
        cclass,czbEFT,czbEFTw = zbetools.get_config(bigconfig=config)

    # Set it all up
    zbetools.setup_outputs(cclass,czbEFT,czbEFTw)
    zbetools.make_configfiles([cclass,czbEFT,czbEFTw], 
                            ['CLASS_configf','zbEFT_configf','zbEFTw_configf'], indfmeta=2)
    
    # Run CLASS (someday maybe use classy)
    command = [czbEFTw['CLASS_exe'],czbEFTw['CLASS_configf'], czbEFTw['CLASS_pre']]
    zbetools.runcommand(command, czbEFTw['logfile'], 
                        czbEFT['PathToLinearPowerSpectrum'], cwd=op.dirname(czbEFTw['CLASS_exe']))
    
    # Run zbEFT
    command = [czbEFTw['zbEFT_exe'],czbEFTw['zbEFT_configf']]
    zbetools.runcommand(command, czbEFTw['logfile'], 
                        czbEFT['PathToOutput'], cwd=czbEFTw['outpath'])

    return czbEFT['PathToOutput']
    
def run_zbEFTClassonly(config):
    """ Runs Pierre Zhang's EFT code in z-space with biased tracers (written in collab
        with Leonardo Senatore and Ashley Perko.

        Inputs
        ------
        config : str or configobj.ConfigObj or dict
                either a string path to the full config file
                or a configobj.ConfigObj or a dictionary containing
                    the parameter names (predetermined - see zbetools.py) and values

        Outputs
        -------
        pathtoresult : string
            local path to the resulting files
    """

    # Get the individual configurations from the one big common config file/object/dict
    if type(config) == str:
        cclass,czbEFT,czbEFTw = zbetools.get_config(bigconfig_file=config)
    else:
        cclass,czbEFT,czbEFTw = zbetools.get_config(bigconfig=config)

    # Set it all up
    zbetools.setup_outputs(cclass,czbEFT,czbEFTw)
    zbetools.make_configfiles([cclass,czbEFT,czbEFTw], 
                            ['CLASS_configf','zbEFT_configf','zbEFTw_configf'], indfmeta=2)
    
    # Run CLASS (someday maybe use classy)
    #command = [czbEFTw['CLASS_exe'],czbEFTw['CLASS_configf'], czbEFTw['CLASS_pre']]
    command = [czbEFTw['CLASS_exe'],czbEFTw['CLASS_configf']] ### Precision accuracy downgraded for faster computation -> OK for computing sigma8
    zbetools.runcommand(command, czbEFTw['logfile'], 
                        czbEFT['PathToLinearPowerSpectrum'], cwd=op.dirname(czbEFTw['CLASS_exe']))
    

    return czbEFT['PathToOutput']

def plks_zbEFT(pars, kvals, path=None):
    """ Calculates the power spectrum multipoles given a dictionary or ConfigObj of 
        input parameters and vector of multipoles.

        Inputs
        ------
        pars : dict or ConfigObj
            the configuration - see the global variables in zbetools.py for defaults
        kvals : 1d numpy array
            vector of k modes
        path : str
            a string to the folder containing the results of running the EFT code
            setting this to None (default) runs the EFT code from scratch
    """

    # Run EFT code if needed
    if path is None:
        path = run_zbEFT(pars)

    # Extract the EFT parameters into a simple array
    bvals = [pars['b1'], pars['b2'], pars['b3'], pars['b4'], 
             pars['b5'], pars['b6'], pars['b7'], 
             pars['b8'], pars['b9'], pars['b10']]

    # Get the multipoles from the EFT terms and parameters
    bp = op.abspath(op.join(path,zbetools.FBSNM))
    kvals, plks = plks_zbEFT_fromfile(bvals, basepath=bp, kvals=kvals, 
                                            maxb=10, dmonly=pars['DM'], 
                                            kren=pars['kren'], lintoo=False)
    return plks

if __name__=='__main__':

    # Need sys to read in command line argument
    import sys

    # Get the path to where the EFT code results are stored
    path = op.abspath(sys.argv[1])

    # Get the multipoles
    kvals = np.linspace(0.01,0.401,40)

    # If the input is just an ini file, run RedshiftBiasEFT_C to get the files
    path = run_zbEFT(path)

    print 'Wrapper took ' + str(round(time.time()-t0,2)) + ' seconds to run EFT code.\nResults of the zbEFT run saved to: '+path
