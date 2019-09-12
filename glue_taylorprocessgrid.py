import os
import numpy as np
import Grid

modpath = os.environ["GROUP_HOME"]
pathpk = os.path.join(modpath, "modgrid")
pathgrid = os.path.join(modpath, "GridsEFT")
print(pathpk)
nruns = 343
lenbatch = 343
gridname = 'z0p32-' + ('-').join(Grid.freepar) # + '-APfcwinSGC' 

gridlin = []
gridloop = []
for i in range(nruns):
    print("Run ", i)
    Plin = np.load(os.path.join(pathpk, "Plin_run%s.npy" % str(i)))
    Ploop = np.load(os.path.join(pathpk, "Ploop_run%s.npy" % str(i)))
    print(Plin.shape, Ploop.shape)
    gridlin.append(Plin)
    gridloop.append(Ploop)
    checklin = (lenbatch == len(Plin))
    checkloop = (lenbatch == len(Ploop))
    if not checklin:
        print("Problema lin!", i, i * lenbatch, Plin[0, 0, -1])
    if not checkloop:
        print("Problema loop!", i, i * lenbatch, Ploop[0, 0, -1])
p1 = np.concatenate(gridlin)
np.save(os.path.join(pathgrid, "TablePlin_%s.npy" % gridname), p1)
p2 = np.concatenate(gridloop)
np.save(os.path.join(pathgrid, "TablePloop_%s.npy" % gridname), p2)

print("I'm done!")
