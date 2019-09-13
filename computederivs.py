import numpy as np
import Grid
import os
import time
import copy
from itertools import combinations
import findiff
import sys

shape = Grid.truegrid.shape
parameters = copy.deepcopy(Grid.parref)
try:
    basedir = './'
    griddir = os.path.join(basedir, "GridsEFT")
except:
    basedir = './'
    griddir = 'grids'
    
gridname = sys.argv[1]

knl = 0.7
km = knl
nd = 4.e-4
# Shape of the crd is now lenpar, gridsize, ... gridsize
# Shape of the PS is now gridsize, ... gridsize, nmult, nk, columns (including the k)
# Since I have the dx, I don't need the coordinates really


# Tested, good
def get_grids(mydir, name=gridname, nmult=3, nout=2):
    # Coordinates have shape (len(freepar), 2 * order_1 + 1, ..., 2 * order_n + 1)
    # order_i is the number of points away from the origin for parameter i
    # The len(freepar) sub-arrays are the outputs of a meshgrid, which I feed to findiff
    # Power spectra needs to be reshaped.
    # crd = np.load(os.path.join(mydir, "Tablecoord_%s.npy" % name))  # Don't need this for uniform grid
    plin = np.load(os.path.join(mydir, "TablePlin_%s.npy" % name))
    # plin = plin.reshape((*shapecrd[1:], nmult, plin.shape[-2] // nmult, plin.shape[-1]))  # This won't work with Python 2 :(
    plin = plin.reshape((shape[1], shape[2], shape[3], shape[4], shape[5], shape[6], 3, plin.shape[-2] // 3, plin.shape[-1]))
    ploop = np.load(os.path.join(mydir, "TablePloop_%s.npy" % name))
    # ploop = ploop.reshape((*shapecrd[1:], nmult, ploop.shape[-2] // nmult, ploop.shape[-1]))  # This won't work with Python 2 :(
    ploop = ploop.reshape((shape[1], shape[2], shape[3], shape[4], shape[5], shape[6], 3, ploop.shape[-2] // 3, ploop.shape[-1]))
    # The output is not concatenated for multipoles since we remove the hexadecapole
    return plin[:, :, :, :, :, :, :nout, :, :], ploop[:, :, :, :, :, :, :nout, :, :]


# Tested, it works well
def get_pder(pi, dx, filename):
    """ Calculates the derivative aroud the Grid.valueref points. Do this only once.
    gridshape is 2 * order + 1, times the number of free parameters
    pi is of shape gridshape, n multipoles, k length, P columns (zeroth being k's)"""
    # Findiff syntax is Findiff((axis, delta of uniform grid along the axis, order of derivative, accuracy))
    t0 = time.time()
    lenpar = len(Grid.valueref)
    idx = Grid.order
    p0 = pi[idx, idx, idx, idx, idx, idx, :, :, :]
    t1 = time.time()
    print("Done p0 in %s sec" % str(t1 - t0))

    dpdx = np.array([findiff.FinDiff((i, dx[i], 1), acc=4)(pi)[idx, idx, idx, idx, idx, idx, :, :, :] for i in range(lenpar)])
    t0 = time.time()
    print("Done dpdx in %s sec" % str(t1 - t0))

    # Second derivatives
    d2pdx2 = np.array([findiff.FinDiff((i, dx[i], 2), acc=4)(pi)[idx, idx, idx, idx, idx, idx, :, :, :] for i in range(lenpar)])
    t1 = time.time()
    print("Done d2pdx2 in %s sec" % str(t1 - t0))

    d2pdxdy = np.array([[i, j, findiff.FinDiff((i, dx[i], 1), (j, dx[j], 1), acc=4)(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
                        for (i, j) in combinations(range(lenpar), 2)]) # if you get an error here, remove np.array([]) (except what's inside!)
    t0 = time.time()
    print("Done d2pdxdy in %s sec" % str(t1 - t0))

    # Third derivatives: we only need it for A_s, so I do this by hand
    d3pdx3 = np.array([findiff.FinDiff((i, dx[i], 3))(pi)[idx, idx, idx, idx, idx, idx, :, :, :] for i in range(lenpar)])
    t1 = time.time()
    print("Done d3pdx3 in %s sec" % str(t1 - t0))
    # d3pdAs3 = (0.5 * (pi[idx + 2, idx, idx, idx, idx, idx, :, :, :] - pi[idx - 2, idx, idx, idx, idx, idx, :, :, :]) + 
    #            pi[idx - 1, idx, idx, idx, idx, idx, :, :, :] - pi[idx + 1, idx, idx, idx, idx, idx, :, :, :]) / dx[0]**3
    # allzeros = np.zeros_like(d3pdAs3)
    # d3pdx3 = np.array([d3pdAs3, allzeros, allzeros, allzeros, allzeros, allzeros])

    # d3pdx2dy = np.array([[i, j, findiff.FinDiff((i, dx[i], 2), (j, dx[j], 1))(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
    #                      for (i, j) in combinations(range(lenpar), 2)])
    # t0 = time.time()
    # print("Done d3pdx2dy in %s sec" % str(t1 - t0))

    # d3pdxdydz = np.array([[i, j, k, findiff.FinDiff((i, dx[i], 1), (j, dx[j], 1), (k, dx[k], 1))(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
    #                      for (i, j, k) in combinations(range(lenpar), 3)])
    # t1 = time.time()
    # print("Done d3pdxdydz in %s sec" % str(t1 - t0))

    # # Fourth derivatives
    # d4pdx4 = np.array([findiff.FinDiff((i, dx[i], 4))(pi)[idx, idx, idx, idx, idx, idx, :, :, :] for i in range(lenpar)])
    # t0 = time.time()
    # print("Done d4pdx4 in %s sec" % str(t1 - t0))

    # d4pdx3dy = np.array([[i, j, findiff.FinDiff((i, dx[i], 3), (j, dx[j], 1))(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
    #                      for (i, j) in combinations(range(lenpar), 2)])
    # t1 = time.time()
    # print("Done d4pdx3dy in %s sec" % str(t1 - t0))

    # d4pdx2dy2 = np.array([[i, j, findiff.FinDiff((i, dx[i], 2), (j, dx[j], 2))(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
    #                      for (i, j) in combinations(range(lenpar), 2)])
    # t0 = time.time()
    # print("Done d4pdx2dy2 in %s sec" % str(t1 - t0))

    # d4pdx2dydz = np.array([[i, j, k, findiff.FinDiff((i, dx[i], 2), (j, dx[j], 1), (k, dx[k], 1))(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
    #                       for (i, j, k) in combinations(range(lenpar), 3)])
    # t1 = time.time()
    # print("Done d4pdx2dydz in %s sec" % str(t1 - t0))

    # d4pdxdydzdw = np.array([[i, j, k, l, findiff.FinDiff((i, dx[i], 1), (j, dx[j], 1),
    #                         (k, dx[k], 1), (l, dx[l], 1))(pi)[idx, idx, idx, idx, idx, idx, :, :, :]]
    #                         for (i, j, k, l) in combinations(range(lenpar), 4)])
    # t1 = time.time()
    # print("Done d4pdxdydzdw in %s sec" % str(t1 - t0))
    
    # Fifth derivatives, only diagonal
    # d5pdx5 = np.array([findiff.FinDiff((i, dx[i], 5))(pi)[idx, idx, idx, idx, idx, idx, :, :, :] for i in range(lenpar)])
    # t0 = time.time()
    # print("Done d5pdx5 in %s sec" % str(t1 - t0))

    allder = (p0, dpdx, d2pdx2, d2pdxdy, d3pdx3) #, d3pdx2dy, d3pdxdydz,
              # d4pdx4, d4pdx3dy, d4pdx2dy2, d4pdx2dydz, d4pdxdydzdw)
    np.save(filename, allder)
    return allder


if __name__ == "__main__":
    print("Let's start!")
    t0 = time.time()
    plingrid, ploopgrid = get_grids(griddir)
    print("Got grids in %s seconds" % str(time.time() - t0))
    dx = Grid.delta
    print("Calculate derivatives of linear PS")
    allderlin = get_pder(plingrid, dx, os.path.join(griddir, "DerPlin_%s.npy" % gridname))

    print("Calculate derivatives of loop PS")
    allderloop = get_pder(ploopgrid, dx, os.path.join(griddir, "DerPloop_%s.npy" % gridname))
