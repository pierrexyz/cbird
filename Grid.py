import numpy as np
import os
import plinear
import pnonlinear
import copy
from scipy.optimize import fsolve
import configobj as cfg

configfile = "config_CMASS.ini"
gridname = "z0p55-A_s-h-omega_cdm-omega_b-n_s-Sum_mnu"
freepar = ["A_s", "h", "omega_cdm", "omega_b", "n_s", "Sum_mnu"]
dx = np.array([0.1, 0.02, 0.04, 0.04, 0.02, 0.25])  # Relative error. DO NOT CHANGE FOR A FIXED GRID!
order = 3  # For the moment I keep this same for everything.

parref = cfg.ConfigObj(configfile)
redshift = float(parref["z_pk"])
valueref = np.array([float(parref[k]) for k in freepar])
delta = dx * valueref
squarecrd = [np.arange(-order, order + 1) for l in freepar]  # list -i_n, ..., +i_n where to evaluate each freepar
truecrd = [valueref[l] + delta[l] * np.arange(-order, order + 1) for l in range(len(freepar))]  # list -i_n, ..., +i_n where to evaluate each freepar
squaregrid = np.array(np.meshgrid(*squarecrd, indexing='ij'))
flattenedgrid = squaregrid.reshape([len(freepar), -1]).T
truegrid = np.array(np.meshgrid(*truecrd, indexing='ij'))

psdatadir = os.path.join("input", "DataSims") 


def get_masses(sum_masses, hierarchy='NH'):
    # a function returning the three masses given the Delta m^2, the total mass, and the hierarchy (e.g. 'IN' or 'IH')
    # Values are in the latest PDG
    # any string containing letter 'n' will be considered as refering to normal hierarchy
    if 'n' in hierarchy.lower():
        # Normal hierarchy massive neutrinos. Calculates the individual
        # neutrino masses from M_tot_NH and deletes M_tot_NH
        delta_m_squared_21 = 7.37e-5
        delta_m_squared_31 = 2.56e-3
        def m1_func(m1, M_tot):
            return M_tot**2 - (m1 + np.sqrt(m1**2 + delta_m_squared_21) + np.sqrt(m1**2 + delta_m_squared_31))**2
        m1, opt_output, success, output_message = fsolve(
            m1_func, sum_masses/3., (sum_masses), full_output=True, xtol=1e-04, maxfev=500)
        m1 = m1[0]
        m2 = (delta_m_squared_21 + m1**2.)**0.5
        m3 = (delta_m_squared_31 + m1**2.)**0.5
        return m1, m2, m3
    else:
        return None


# First step: compute PS. Tested, good
def CompPterms(pardict, kmin=0.01, kmax=0.31):
    """Given a parameter dictionary, a kmin (can be None) and a kmax (can be None),
    returns Plin, Ploop concatenated for multipoles, shape (lenk * 3, columns).
    The zeroth column are the k"""
    parlinear = copy.deepcopy(pardict)
    parnonlinear = copy.deepcopy(pardict)
    # print("As", parlinear["A_s"])
    if "ln10^{10}A_s" in parlinear.keys():
        print("we have the log here")
    m1, m2, m3 = get_masses(float(parlinear['Sum_mnu']))
    parlinear['m_ncdm'] = str(m1) + ', ' + str(m2) + ', ' + str(m3)
    this_plin = plinear.LinearPower(parlinear, np.logspace(np.log10(1.e-7), np.log10(2.01), 150))
    this_plin.compute()
    parnonlinear["PathToLinearPowerSpectrum"] = os.path.join(parnonlinear["PathToOutput"], "class_pk.dat")
    try:
        if not os.path.isdir(parnonlinear["PathToOutput"]):
            os.makedirs(parnonlinear["PathToOutput"])
    except IOError:
        print("Cannot create directory: %s" % parnonlinear["PathToOutput"])
    pnl = pnonlinear.NonLinearPower(parnonlinear, kmin, kmax)
    pnl.compute()
    Plin = pnl.get_Plin_resum()
    Ploop = pnl.get_Ploop_resum()
    return Plin, Ploop
