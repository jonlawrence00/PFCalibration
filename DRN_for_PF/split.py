#!/bin/python3 -u

from Validation import Validation
import pickle
import numpy as np
import math
import awkward as ak


SPLIT = 0.8

folder = "2018_UL_Particle_MC_30M"
folder = "2018_Gamma_Jet/combPhotons"
folder = "2018_UL_Photon_Refined"
yayPhotons = True

def write(valid, train, pref):
    folderPath = '/home/rusack/shared/pickles/%s'%folder
    f = '%s/%s_valididx.pickle'%(folderPath,pref)
    print(f)
    with open(f, 'wb') as f:
        pickle.dump(valid.nonzero()[0], f)

    f = '%s/%s_trainidx.pickle'%(folderPath, pref)
    with open(f, 'wb') as f:
        pickle.dump(train.nonzero()[0], f)

    print("%s:"%pref)
    print("\t%d total points"%(np.sum(valid) + np.sum(train)))
    print("\t%d validation points"%np.sum(valid))
    print("\t%d training points"%np.sum(train))
    print("\tsplit of %f"%(np.sum(train)/(np.sum(train) + np.sum(valid)) ))
    print()

mod = Validation(".", folder, photons=yayPhotons)
mod._loadVariable('subdet')
mod._loadVariable('iZ')
#mod._loadVariable('R9')

mod.data['iZ'] = ak.to_numpy(mod.data['iZ'])

EB = mod.data['subdet'] == 1
EE = mod.data['subdet'] == 0
EEL = mod.data['iZ'] == -1
EER = mod.data['iZ'] == 1

#lowR9 = mod.data['R9'] < 0.96
#highR9 = ~lowR9

#EB_lowR9 = np.logical_and(EB, lowR9)
#EB_highR9 = np.logical_and(EB, highR9)

#EE_lowR9 = np.logical_and(EE, lowR9)
#EE_highR9 = np.logical_and(EE, highR9)

trainidxEB = np.random.choice(EB.nonzero()[0],
        (int)( np.rint(np.sum(EB)*SPLIT)),
        replace=False)
trainidxEE = np.random.choice(EE.nonzero()[0],
        (int)(np.rint(np.sum(EE)*SPLIT)),
        replace=False)
valid = np.ones(len(EE), dtype=bool)

valid[trainidxEB] = False
valid[trainidxEE] = False
train = ~valid

EBEEL = np.logical_or(EB, EEL)
EBEER = np.logical_or(EB, EER)

write(  valid, 
        train,
        'both')

#write(  np.logical_and(lowR9,valid), 
#        np.logical_and(lowR9,train),
#        'both_lowR9')

#write(  np.logical_and(highR9,valid), 
#        np.logical_and(highR9,train),
#        'both_highR9')

write(  np.logical_and(EB,valid), 
        np.logical_and(EB,train),
        'EB')

write(  np.logical_and(EE,valid), 
        np.logical_and(EE,train),
        'EE')

write(  np.logical_and(EBEEL,valid), 
        np.logical_and(EBEEL,train),
        'EBEEL')

write(  np.logical_and(EBEER,valid), 
        np.logical_and(EBEER,train),
        'EBEER')

#write(  np.logical_and(EE_lowR9,valid), 
#        np.logical_and(EE_lowR9,train),
#        'EE_lowR9')

#write(  np.logical_and(EB_lowR9,valid), 
#        np.logical_and(EB_lowR9,train),
#        'EB_lowR9')

#write(  np.logical_and(EE_highR9,valid), 
#        np.logical_and(EE_highR9,train),
#        'EE_highR9')

#write(  np.logical_and(EB_highR9,valid), 
       # np.logical_and(EB_highR9,train),
       # 'EB_highR9')
