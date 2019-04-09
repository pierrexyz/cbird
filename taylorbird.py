#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import os.path as opa
from scipy.interpolate import interp1d

import tayloregg

#THIS_PATH  =  '/home/pchonje/Documents/EFTofLSS/cbird/'
THIS_PATH = opa.dirname(opa.abspath(__file__)) 
OUTPATH = opa.abspath(opa.join(THIS_PATH, 'taylornest/'))

def load_taylor(which):
    k = np.load(opa.join(OUTPATH, "%s/k.npy"%which))
    Plin_ref = np.load(opa.join(OUTPATH, "%s/Plin_ref.npy"%which))
    Ploop_ref = np.load(opa.join(OUTPATH, "%s/Ploop_ref.npy"%which))
    fdlin = np.load(opa.join(OUTPATH, "%s/fdlin.npy"%which))
    sdlin = np.load(opa.join(OUTPATH, "%s/sdlin.npy"%which))
    fdloop = np.load(opa.join(OUTPATH, "%s/fdloop.npy"%which))
    sdloop =np.load(opa.join(OUTPATH, "%s/sdloop.npy"%which))
    return k, Plin_ref, Ploop_ref, fdlin, fdloop, sdlin, sdloop

def taylorbird(cosmoref, cosmotarget, tayloregg):
    k, Plin_ref, Ploop_ref, fdlin, fdloop, sdlin, sdloop = tayloregg
    delta = cosmotarget - cosmoref
    Plin = Plin_ref + np.einsum('clpk,c->lpk', fdlin, delta) + 0.5 * np.einsum('c,d,cdlpk->lpk', delta, delta, sdlin)
    Ploop = Ploop_ref + np.einsum('clpk,c->lpk', fdloop, delta) + 0.5 *  np.einsum('c,d,cdlpk->lpk', delta, delta, sdloop)
    #Plin = ( Plin_ref.T + fdlin.T.dot(delta) + 0.5 * sdlin.T.dot(delta).dot(delta) ).T
    #Ploop = ( Ploop_ref.T + fdloop.T.dot(delta) + 0.5 * sdloop.T.dot(delta).dot(delta) ).T
    return k, Plin, Ploop

nd = 3e-4
km = 0.7
knl = 0.7

def computePS(bs, datalin, dataloop, setkin, setkout):
    
    datalin0,datalin2,datalin4 = datalin
    data0,data2,data4 = dataloop[:,:18,:]
    b1,c2,b3,c4,b5,b6,b7,b8,b9,b10,e1,e2 = bs

    # Add a PCA between b2 and b4 to disantangle the degeneracy:
    b2 = 0.5 * (c2 + c4)
    b4 = 0.5 * (c2 - c4)

    # the columns of the Ploop data files.
    cvals = np.array([1., b1, b2, b3, b4, b1*b1, b1*b2, b1*b3, b1*b4, b2*b2, b2*b4, b4*b4, b1*b5/knl**2, b1*b6/km**2, b1*b7/km**2, b5/knl**2, b6/km**2 ,b7/km**2])

    P0 = interp1d( setkin, np.dot(cvals,data0) + datalin0[0]+b1*datalin0[1]+b1*b1*datalin0[2] + b1*e1*setkin**2/km**2*data0[13] )(setkout) + b8/nd + b9/nd/km**2 * setkout**2 
    P2 = interp1d( setkin, np.dot(cvals,data2) + datalin2[0]+b1*datalin2[1]+b1*b1*datalin2[2] + e2*setkin**2/km**2*data2[16] )(setkout) + b10/nd/km**2 * setkout**2 
    P4 = interp1d( setkin, np.dot(cvals,data4) + datalin4[0]+b1*datalin4[1]+b1*b1*datalin4[2] )(setkout)
    
    return np.array([P0,P2,P4])

def computeTaylorPS(bs, cosmoref, cosmotarget, tayloregg, kout):
    k, Plin, Ploop = taylorbird(cosmoref, cosmotarget, tayloregg)
    return computePS(bs, Plin, Ploop, k, kout)

if __name__ ==  "__main__":

    simtype = 'ChallengeA'
    z_pk, cosmoref, cosmodelta = tayloregg.import_cosmoref_from_DataFrameCosmosim(simtype)
    Om_fid = cosmoref[1]
    TableNkmu = np.loadtxt( "/home/pchonje/Documents/EFTofLSS/cbird/cbird-dev/input/Binning/NkmuChallengeA.dat" ).T

    bs = np.array([2.0573575829,   1.6107766447,   -1.4636660539,  -7.2081714231,  -1.4056024451,  -13.8955428899, -0.5069869722,  0, 0 , 0, 0, 0])

    #Plin_test = np.load(opa.join(OUTPATH, "Plin_test.npy"))
    #Ploop_test = np.load(opa.join(OUTPATH, "Ploop_test.npy"))
    #cosmotest = np.copy(cosmoref + cosmodelta)

    cosmotest = np.copy(cosmoref)

    cosmotest[0] += cosmodelta[0]
    cosmotest[2] += 6*cosmodelta[2]
    cosmotest[1] += 6*cosmodelta[1]

    egg = load_taylor('ChallengeA')
    k, Plin_taylor, Ploop_taylor = taylorbird(cosmoref, cosmotest, egg)
    PS_taylor = computePS(bs, Plin_taylor, Ploop_taylor, k, k)

    Plin_test, Ploop_test = tayloregg.CompPterms(z_pk, cosmotest, k, Om_fid, TableNkmu)
    PS_test = computePS(bs, Plin_test, Ploop_test, k, k)

    import matplotlib.pyplot as plt

    plt.figure(1)
    plt.semilogy(k, PS_taylor[0], 'g')
    plt.semilogy(k, PS_test[0], 'b')

    plt.figure(2)
    plt.plot(k, PS_taylor[0]/PS_test[0])

    plt.show()
