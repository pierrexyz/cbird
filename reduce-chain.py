import numpy as np
import os.path as opa

def reduce_chain(names, simtype, marg, stoch, kmaxbisp, foldername, kmax=0.25):
    c = 0
    chains = [0]*len(names)
    if kmaxbisp > 0:
        simtype += 'withBispkmax%s'%kmaxbisp
    if marg:
        simtype+='gaussMarg'

    simtype+= 'stoch%s.0'%stoch
    for box in names:
        chainpure = np.load( '%s/samplerchain%sbox_%skmax_%s.npy' % (folder, simtype, box, kmax) )
        stacked = chainpure[:,9*chainpure.shape[1]/10::1000,:].reshape((-1, 11+stoch)) 
        np.save( '%s/samples%sbox_%skmax_%s.npy' % (folder, simtype, box, kmax), stacked )
        lnlike = np.load( '%s/lnlikechain%sbox_%skmax_%s.npy' % (folder, simtype, box, kmax) )
        id1,idx=np.unravel_index(lnlike.argmax(), lnlike.shape)
        print -2*lnlike.max()
        bestfit = chainpure[id1, idx, :]
        np.savetxt( '%s/bestfit%sbox_%skmax_%s.txt' % (folder, simtype, box, kmax), bestfit )
        np.savetxt( '%s/maxlnlike%sbox_%skmax_%s.txt' % (folder, simtype, box, kmax), [-2*lnlike.max()] )
    return


###############
folder = '/exports/pierre/cbird/mcnest/'
marg = False
kmaxbisp = 0
###############
runtype = 'Challenge-shifted'

box = 'A'

stochoption = [0, 1, 2, 3, 4]

kmaxtab = [0.2, 0.25]
#kmaxtab = [0.15, 0.2, 0.25, 0.3]
for kmax in kmaxtab:
    for stoch in stochoption:
        reduce_chain(box, runtype, marg, stoch, kmaxbisp, folder, kmax=kmax)