#!/bin/python3 -u

from Validation import Validation
import pickle
import numpy as np
import math

folder = '2018_UL_ZEE_Data'

def write(idx, pref):
    fname = '%s/all_%s_valididx.pickle'%(folder, pref)
    print("dumping %s"%pref)
    print("\t%d poibnts"%len(idx))
    with open(fname, 'wb') as f:
        pickle.dump(idx, f)

mod = Validation(None, folder)
mod._loadVariable("R9")
mod._loadVariable('subdet')
mod._loadVariable('etareco')
mod._loadVariable('energy_ecal_mustache')
mod.data['pt'] = mod.data['energy_ecal_mustache']/np.cosh(mod.data['etareco'])

EB = mod.data['subdet'] == 0
EE = mod.data['subdet'] == 1

ptcut = mod.data['pt'] > 20
EB = np.logical_and(EB, ptcut)
EE = np.logical_and(EE, ptcut)

EBetacut = mod.data['etareco'] < 1.44
EB = np.logical_and(EBetacut, EB)

EEetacut = np.logical_and(mod.data['etareco'] > 1.57, mod.data['etareco'] < 2.5)
EE = np.logical_and(EEetacut, EE)

lowR9 = mod.data['R9'] < 0.96
highR9 = ~lowR9

EBEB = []
EBEE = []
EEEE = []

EBEB_low = []
EBEB_high = []

EBEE_low = []
EBEE_high = []

EEEE_low = []
EEEE_high = []

for i in range(len(EE)//2):
    e1 = 2*i
    e2 = e1+1
    if EB[e1] and EB[e2]: #EBEB
        EBEB.append(e1)
        EBEB.append(e2)

        if lowR9[e1] and lowR9[e2]:
            EBEB_low.append(e1)
            EBEB_low.append(e2)
        elif highR9[e1] and highR9[e2]:
            EBEB_high.append(e1)
            EBEB_high.append(e2)

    elif ((EB[e1] and EE[e2]) or (EB[e2] and EE[e1])):
        EBEE.append(e1)
        EBEE.append(e2)

        if lowR9[e1] and lowR9[e2]:
            EBEE_low.append(e1)
            EBEE_low.append(e2)
        elif highR9[e1] and highR9[e2]:
            EBEE_high.append(e1)
            EBEE_high.append(e2)

    elif EE[e1] and EE[e2]:
        EEEE.append(e1)
        EEEE.append(e2)

        if lowR9[e1] and lowR9[e2]:
            EEEE_low.append(e1)
            EEEE_low.append(e2)
        elif highR9[e1] and highR9[e2]:
            EEEE_high.append(e1)
            EEEE_high.append(e2)

write(EBEB,         'EBEB')
write(EBEB_low,     'EBEB_low')
write(EBEB_high,    'EBEB_high')

write(EBEE,         'EBEE')
write(EBEE_low,     'EBEE_low')
write(EBEE_high,    'EBEE_high')

write(EEEE,         'EEEE')
write(EEEE_low,     'EEEE_low')
write(EEEE_high,    'EEEE_high')
