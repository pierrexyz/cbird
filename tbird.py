import numpy as np
import Grid

def load_pder(filename):
    allder = np.load(filename, allow_pickle=True)
    return allder


# It works well
def get_PSTaylor(dtheta, derivatives):
    # Shape of dtheta: number of free parameters
    # Shape of derivatives: tuple up to third derivative where each element has shape (num free par, multipoles, lenk, columns)
    t1 = np.einsum('p,pmkb->mkb', dtheta, derivatives[1])
    t2diag = np.einsum('p,pmkb->mkb', dtheta**2, derivatives[2])
    t2nondiag = np.sum([dtheta[d[0]] * dtheta[d[1]] * d[2] for d in derivatives[3]], axis=0)
    #t3diag = np.einsum('p,pmkb->mkb', dtheta**3, derivatives[4])
    t3diag = np.einsum('p,pmkb->mkb', np.array([dtheta[0]**3, 0, 0, 0, 0, 0 * dtheta[5]**3]), derivatives[4])
    #t3semidiag = np.sum([dtheta[d[0]]**2 * dtheta[d[1]] * d[2] for d in derivatives[5]], axis=0)
    #t3nondiag = np.sum([dtheta[d[0]] * dtheta[d[1]] * dtheta[d[2]] * d[3] for d in derivatives[6]], axis=0)
    #t4diag = np.einsum('p,pmkb->mkb', dtheta**4, derivatives[7])
    #t4semidiag1 = np.sum([dtheta[d[0]]**3 * dtheta[d[1]] * d[2] for d in derivatives[8]], axis=0)
    #t4semidiag2 = np.sum([dtheta[d[0]]**2 * dtheta[d[1]]**2 * d[2] for d in derivatives[9]], axis=0)
    #t4semidiag3 = np.sum([dtheta[d[0]]**2 * dtheta[d[1]] * dtheta[d[2]] * d[3] for d in derivatives[10]], axis=0)
    # t4nondiag = np.sum([dtheta[d[0]] * dtheta[d[1]] * dtheta[d[2]] * dtheta[d[3]] * d[4] for d in derivatives[11]], axis=0)
    # t5diag = np.einsum('p,pmkb->mkb', dtheta**5, derivatives[12])
    allPS = (derivatives[0] + t1 + 0.5 * t2diag + t2nondiag +
             t3diag / 6.) # + t3semidiag / 2. + t3nondiag)
             # t4diag / 24.  + t4semidiag1 / 6. + t4semidiag2 / 4. + t4semidiag3 / 2. + t4nondiag)  # + t5diag / 120.)
    return allPS


# Nothing much to test, helper function
def get_bias(d_fit):
    vlin = np.array([0., 1., d_fit['b1'], d_fit['b1']**2])
    vloop = np.array([0., 1., d_fit['b1'], d_fit['b2'], d_fit['b3'], d_fit['b4'],
                      d_fit['b1']**2, d_fit['b1'] * d_fit['b2'], d_fit['b1'] * d_fit['b3'], d_fit['b1'] * d_fit['b4'],
                      d_fit['b2']**2, d_fit['b2'] * d_fit['b4'], d_fit['b4']**2,
                      d_fit['b1'] * d_fit['b5'] / knl**2, d_fit['b1'] * d_fit['b6'] / km**2,
                      d_fit['b1'] * d_fit['b7'] / km**2, d_fit['b5'] / knl**2, d_fit['b6'] / km**2, d_fit['b7'] / km**2])
    return vlin, vloop


# Tested
def get_PSbias(plin, ploop, dfit):
    """Given a dictionary of biases, gets the k's and the full PS.
    The shape of the PS is (nmultipoles, len(k))"""
    vlin, vloop = get_bias(dfit)
    kin = plin[0, :, 0]
    PS = np.einsum('c,mkc->mk', vlin, plin) + np.einsum('c,mkc->mk', vloop, ploop)
    PS[0] = PS[0] + dfit['b8'] / nd + dfit['b9'] / nd / km**2 * kin**2
    PS[1] = PS[1] + dfit['b10'] / nd / km**2 * kin**2
    return kin, PS