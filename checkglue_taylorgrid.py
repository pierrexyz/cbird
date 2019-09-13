import os
import sys
import numpy as np
import Grid

basedir = './'
pathpk = os.path.join(basedir, "Pk")
pathgrid = os.path.join(basedir, "GridsEFT")
gridname = 'taylorgrid-' + ('-').join(Grid.freepar)

nruns = int(sys.argv[1])
nranks = int(sys.argv[2])
lenbatch = int(sys.argv[3])
ntot = nruns * nranks

linfailed = []
loopfailed = []
for i in range(nruns):
    print(i)
    for j in range(nranks):
        checklin = os.path.isfile(os.path.join(
            pathpk, "Plin_run%d_rank%d.npy" % (i, j)))
        checkloop = os.path.isfile(os.path.join(
            pathpk, "Ploop_run%d_rank%d.npy" % (i, j)))
        if not checklin:
            print("Failed linear run %d, rank %d" % (i, j))
            linfailed.append((i, j))
        else:
            Plin = np.load(os.path.join(pathpk, "Plin_run%d_rank%d.npy" % (i, j)))
            if (lenbatch != len(Plin)):
                print("Failed length linear run %d, rank %d" % (i, j))
                linfailed.append((i, j))
        if not checkloop:
            print("Failed loop run %d, rank %d" % (i, j))
            loopfailed.append((i, j))
        else:
            Ploop = np.load(os.path.join(pathpk, "Ploop_run%d_rank%d.npy" % (i, j)))
            if (lenbatch != len(Ploop)):
                print("Failed length loop run %d, rank %d" % (i, j))
                loopfailed.append((i, j))

print("Linear failed: %d over %d, %f %%" %
      (len(linfailed), ntot, 100 * float(len(linfailed)) / ntot))
print("Loop failed: %d over %d, %f %%" %
      (len(loopfailed), ntot, 100 * float(len(loopfailed)) / ntot))

if (len(linfailed) + len(loopfailed)) > 0:
    np.save("allfailed.npy", np.array(linfailed))
    raise Exception("Some processes have failed!")

gridlin = []
gridloop = []
for i in range(nruns):
    print("Run ", i)
    for j in range(nranks):
        Plin = np.load(os.path.join(
            pathpk, "Plin_run%d_rank%d.npy" % (i, j)))
        Ploop = np.load(os.path.join(
            pathpk, "Ploop_run%d_rank%d.npy" % (i, j)))
        gridlin.append(Plin[:, :, :-1])
        gridloop.append(Ploop[:, :, :-1])
        # checklin = ((i * nranks + j) * lenbatch == int(Plin[0, 0, -1))
        # checkloop = ((i * nranks + j) * lenbatch == int(Ploop[0, 0, -1]))
        checklin = (lenbatch == len(Plin))
        checkloop = (lenbatch == len(Ploop))
        if not checklin:
            print("Problem in linear PS: ", i, j, (i * nranks + j) * lenbatch, Plin[0, 0, -1])
        if not checkloop:
            print("Problem in loop PS: ", i, j, (i * nranks + j) * lenbatch, Ploop[0, 0, -1])

np.save(os.path.join(pathgrid, "Tablecoord_%s.npy" % gridname), Grid.truegrid)
g1 = np.concatenate(gridlin)
np.save(os.path.join(pathgrid, "TablePlin_%s.npy" % gridname), g1)
g2 = np.concatenate(gridloop)
np.save(os.path.join(pathgrid, "TablePloop_%s.npy" % gridname), g2)
