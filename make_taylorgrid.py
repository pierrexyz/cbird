from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

import numpy as np
import os
import sys
import Grid
import copy
import time

time.sleep(4)

basedir = './'
OUTPATH = os.path.join(basedir, "output")
outpk = os.path.join(basedir, "Pk")


ncores = size
nrun = int(sys.argv[1])
runs = int(sys.argv[2])
lenrun = int(len(Grid.flattenedgrid) / runs)
thetarun = Grid.flattenedgrid[nrun * lenrun:(nrun+1) * lenrun]
Ntot = len(thetarun)
sizered = Ntot/ncores

arrayred = thetarun[rank * sizered:(rank+1) * sizered]

freepar = Grid.freepar
print("lenrun, sizered", lenrun, sizered)
allPlin = []
allPloop = []
for i, theta in enumerate(arrayred):
    parameters = copy.deepcopy(Grid.parref)
    truetheta = Grid.valueref + theta * Grid.delta
    idx = nrun * lenrun + rank * sizered + i
    print("nrun, rank, i", nrun, rank, i)
    parameters["PathToOutput"] = os.path.join(OUTPATH, 'output' + str(nrun) + str(rank) + str(i))
    for k, var in enumerate(freepar):
        parameters[var] = truetheta[k]
    # parameters['h'] = 1/parameters['invh']
    Plin, Ploop = Grid.CompPterms(parameters)
    idxcol = np.full([Plin.shape[0], 1], idx)
    allPlin.append(np.hstack([Plin, idxcol]))
    allPloop.append(np.hstack([Ploop, idxcol]))
    if (i == 0) or ((i+1) % 100 == 0):
        print("theta check: ", Grid.flattenedgrid[idx], theta, truetheta)
        # np.save(os.path.join(outpk, "temp", "Plin_run%s_rank%si%s.npy" % (str(nrun), str(rank), str(i))), np.array(allPlin))
        # np.save(os.path.join(outpk, "temp", "Ploop_run%s_rank%si%s.npy" % (str(nrun), str(rank), str(i))), np.array(allPloop))
    np.save(os.path.join(outpk, "Plin_run%s_rank%s.npy" % (str(nrun), str(rank))), np.array(allPlin))
    np.save(os.path.join(outpk, "Ploop_run%s_rank%s.npy" % (str(nrun), str(rank))), np.array(allPloop))
