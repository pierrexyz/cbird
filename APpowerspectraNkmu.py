import numpy as np
import scipy
import scipy.interpolate as sp

def check_if_multipoles_k_array(setk):
    return setk[len(setk)/3] == setk[0] 

def changetoAPbinning(Pk,setkin,setkout,qperp,qpar,TableNkmu,l68 = None):
    _,kmean,mucent,nkmu = TableNkmu #import the data from the sims. mucent are central values, kmean the mean.

    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]

    # Add l=6,8 contribution
    if type(l68) != type(None):
        Pkloc = np.concatenate([Pk,l68])
    else:
        Pkloc = Pk

    Pkint = sp.interp1d(setkin,Pkloc,axis=-1,kind='cubic',bounds_error = False,fill_value = 'extrapolate')

    # Define the grid with the right kmax and kmin and reshape into (k,mu)
    kmin = setkout.min()
    kmax = setkout.max()
    kmeanx = kmean[(kmean>= kmin)&(kmean <= kmax)]
    mucentx = mucent[(kmean>= kmin)&(kmean <= kmax)]
   
    Nbink = len(kmeanx)/100
    Nbinmu = 100
    
    kgrid = kmeanx.reshape((Nbink,Nbinmu))   
    mugrid = mucentx.reshape((Nbink,Nbinmu))
    
    # Reshape N(k,mu) on the grid with right kmin and kmax
    nkmux = nkmu[(kmean>= kmin)&(kmean <= kmax)]
    nkgrid = nkmux.reshape((Nbink,Nbinmu)) 
    
    # Interpolate the mu part of N(k,mu)
    nkgridint = sp.interp1d(mugrid[0,:],nkgrid,axis = 1,kind = 'nearest',bounds_error = False,fill_value = 'extrapolate')
    
    # New array of mu with more points (better precision for the integration)
    muacc = np.linspace(0.,1.,1000)
    
    mugrid,kgrid = np.meshgrid(muacc,np.unique(kmeanx))  
    
    # AP factors
    F = float(qpar/qperp)
    k = kgrid/qperp*(1+mugrid**2*(F**-2-1))**0.5
    mup = mugrid/F*(1+mugrid**2*(F**-2-1))**-0.5
    
    # Goes from the multipoles back to P(k,mu) and apply AP
    if type(l68) == type(None):
        arrayLegendremup = nkgridint(muacc)*np.array([scipy.special.legendre(0)(mup),
                                scipy.special.legendre(2)(mup),
                                scipy.special.legendre(4)(mup)])
    else:
        #print ('l68 Legendre')
        arrayLegendremup = nkgridint(muacc)*np.array([scipy.special.legendre(0)(mup),
                                    scipy.special.legendre(2)(mup),
                                    scipy.special.legendre(4)(mup),
                                    scipy.special.legendre(6)(mup),
                                    scipy.special.legendre(8)(mup)])

    arrayLegendremugrid = np.array([2*(2*0+1.)/(2*qperp**2*qpar)*scipy.special.legendre(0)(mugrid),
                              2*(2*2.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(2)(mugrid),
                              2*(2*4.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(4)(mugrid)])

    Pkmu = np.einsum('lpkm,lkm->pkm', Pkint(k), arrayLegendremup)

    # Normalization for N(k,mu)dmu
    nk = np.trapz(nkgridint(muacc),x = muacc,axis = 1)
    
    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    Integrandmu = np.einsum('pkm,lkm->lpkm', Pkmu, arrayLegendremugrid)
    Pk_AP = np.trapz(Integrandmu, x = mugrid, axis = -1) / nk

    # interpolate on the wanted k-array for output
    Pk_AP_out =(sp.interp1d(np.unique(kmeanx),Pk_AP,axis=-1,bounds_error = False,fill_value = 'extrapolate'))(setkout) 
    
    return Pk_AP_out

def changetoAPnobinning(Pk, setkin, setkout, qperp, qpar, nbinsmu = 500):

    muacc = np.linspace(0.,1.,nbinsmu)

    # Check the k-arrays are in the right format (not concatenated for multipoles)
    if check_if_multipoles_k_array(setkin):
        setkin = setkin[:len(setkin)/3]
    if check_if_multipoles_k_array(setkout):
        setkout = setkout[:len(setkout)/3]

    # Interpolate the multipoles
    Pkint = scipy.interpolate.interp1d(setkin, Pk, axis=-1, kind='cubic', bounds_error = False, fill_value = 'extrapolate')

    # Define the grid with the right kmax and kmin and reshape into (k,mu)
    kgrid, mugrid = np.meshgrid(setkout, muacc, indexing='ij')
        
    # AP factors
    F = float(qpar/qperp)
    k = kgrid/qperp*(1+mugrid**2*(F**-2-1))**0.5
    mup = mugrid/F*(1+mugrid**2*(F**-2-1))**-0.5

   # Goes from the multipoles back to P(k,mu) and apply AP
    arrayLegendremup = np.array([scipy.special.legendre(0)(mup),
                                scipy.special.legendre(2)(mup),
                                scipy.special.legendre(4)(mup)]) 
        
    arrayLegendremugrid = np.array([2*(2*0+1.)/(2*qperp**2*qpar)*scipy.special.legendre(0)(mugrid),
                                  2*(2*2.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(2)(mugrid),
                                  2*(2*4.+1.)/(2*qperp**2*qpar)*scipy.special.legendre(4)(mugrid)])

    #print(k.shape, Pkint(k).shape, arrayLegendremugrid.shape) 
    # Pkint(k).shape: (multipoles, power spectra, ks, mus): (lpkm)
    Pkmu = np.einsum('lpkm,lkm->pkm', Pkint(k), arrayLegendremup)
    
    # Back to multipoles (factor of 2 because we integrate an even function from 0 to 1 instead of -1 to 1)
    #print(Pkmu.shape, arrayLegendremugrid.shape) 
    # Pkmu.shape: (power spectra, ks, mus): (pkm)
    Integrandmu = np.einsum('pkm,lkm->lpkm', Pkmu, arrayLegendremugrid)
    
    Pk_AP = np.trapz(Integrandmu, x = mugrid, axis = -1)
    #print (Pk_AP.shape)
    #Pk_AP.shape: (multipoles, power spectra, ks)

    return Pk_AP