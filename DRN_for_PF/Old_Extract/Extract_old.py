import ROOT    
from tqdm import tqdm
import numpy as np
from time import time
import pickle

######### Fixed values for scaling the input features v3 ##################
Eta_Min = -3.0
Eta_Max = 3.0

Phi_Min = -np.pi
Phi_Max = np.pi

Z_Min = -330.0
Z_Max = 330.0

X_Min = -150.0
X_Max = 150.0

Y_Min = -150.0
Y_Max = 150.0

Energy_Min = 0.0
Energy_Max = 300.

Scaled_ECAL_Min = 0.0
Scaled_ECAL_Max = 250.0

Scaled_ES_Min = 0.0
Scaled_ES_Max = 0.07

iEta_Min = -85
iEta_Max = 85

iX_Min = 1
iX_Max = 100

iY_Min = 1
iY_Max = 100

iPhi_Min = 1
iPhi_Max = 360

iZ_Min = -1
iZ_Max = 1

########## Conversion from nTuple events to data rows #####################

def Scale_feature(feature_value, feature_value_min, feature_value_max):
        Numerator = feature_value - feature_value_min
        Denominator = feature_value_max - feature_value_min
        if Denominator==0 :
                print("Caution: Max - Min of this input feature is 0!!! Changing to 1.")
                Denominator = 1
        return  Numerator/Denominator

def ElectronToFeatures(event, i, ES='scaled', coords='proj', fracs='yes'):
    '''
    event: the event
    i: the index of the electron (0 or 1)
    ES: one of 'yes', 'no', 'scaled'
        'yes' : use ES, same scale as ECAL
        'no' : don't use ES
        'scaled' : use ES, with its own scale
    coords: one of 'proj', 'cart', 'local'
        'proj' : eta, phi, z
        'cart' : x, y, z
        'local': ieta, iphi sign(Z) (barrel) or iX iY sign(Z) (endcaps)
            not compatible with use of ES.
            ie ES must be == 'no'
    fracs: one of 'yes', 'no', 'mult'
        'yes': inclue fracs
        'no' : don't include fracs
        'mult' : multiply energies by fracs
    '''

    if ES == 'no':
        useES = False
        scaled = False
    elif ES == 'yes':
        useES = True
        scaled = False
    elif ES == 'scaled':
        useES = True
        scaled = True

    if i==0:
        if coords!='local':
            if coords=='proj':
                eta = np.asarray(event.Hit_Eta_Ele1)
                phi = np.asarray(event.Hit_Phi_Ele1)
            else:
                x = np.asarray(event.Hit_X_Ele1)
                y = np.asarray(event.Hit_Y_Ele1)

        else:
            iEta = np.asarray(event.iEtaEle1)
            iPhi = np.asarray(event.iPhiEle1)

        z = np.asarray(event.Hit_Z_Ele1)

        energy = np.asarray(event.RecHitEnEle1)

        if fracs=='yes':
            fraction = np.asarray(event.RecHitFracEle1)
        elif fracs=='mult':
            energy *= np.asarray(event.RecHitFracEle1)

        if scaled:
            ES = np.zeros(energy.shape)

        if useES:
            if coords == 'proj':
                eta_ES = np.asarray(event.Hit_ES_Eta_Ele1)
                phi_ES = np.asarray(event.Hit_ES_Phi_Ele1)
            else:
                x_ES = np.asarray(event.Hit_ES_X_Ele1)
                y_ES = np.asarray(event.Hit_ES_Y_Ele1)

            z_ES = np.asarray(event.Hit_ES_Z_Ele1)
            energy_ES = np.asarray(event.ES_RecHitEnEle1)
            if fracs=='yes':
                fraction_ES = np.ones(energy_ES.shape)
            if scaled:
                ES_ES = np.ones(energy_ES.shape)
    else:
        if coords!='local':
            if coords == 'proj':
                eta = np.asarray(event.Hit_Eta_Ele2)
                phi = np.asarray(event.Hit_Phi_Ele2)
            else:
                x = np.asarray(event.Hit_X_Ele2)
                y = np.asarray(event.Hit_Y_Ele2)
        else:
            iEta = np.asarray(event.iEtaEle2)
            iPhi = np.asarray(event.iPhiEle2)

        z = np.asarray(event.Hit_Z_Ele2)

        energy = np.asarray(event.RecHitEnEle2)
        if fracs == 'yes':
            fraction = np.asarray(event.RecHitFracEle2)
        elif fracs == 'mult':
            energy *= np.asarray(event.RecHitFracEle2)
        if scaled:
            ES = np.zeros(energy.shape)

        if useES:
            if coords == 'proj':
                eta_ES = np.asarray(event.Hit_ES_Eta_Ele2)
                phi_ES = np.asarray(event.Hit_ES_Phi_Ele2)
            else:
                x_ES = np.asarray(event.Hit_ES_X_Ele2)
                y_ES = np.asarray(event.Hit_ES_Y_Ele2)

            z_ES = np.asarray(event.Hit_ES_Z_Ele2)
            energy_ES = np.asarray(event.ES_RecHitEnEle2)
            if fracs=='yes':
                fraction_ES = np.ones(energy_ES.shape)
            if scaled:
                ES_ES = np.ones(energy_ES.shape)

    if scaled:
        energy_scaled = Scale_feature(energy, Scaled_ECAL_Min, Scaled_ECAL_Max)
        if useES:
            energy_ES_scaled = Scale_feature(energy_ES, Scaled_ES_Min, Scaled_ES_Max)
            energy_scaled = np.concatenate((energy_scaled, energy_ES_scaled))
    else:
        if useES:
            energy = np.concatenate((energy, energy_ES))
        energy_scaled = Scale_feature(energy, Scaled_ECAL_Min, Scaled_ECAL_Max)

    if useES:
        if coords == 'proj':
            eta = np.concatenate((eta, eta_ES))
            phi = np.concatenate((phi, phi_ES))
        else:
            x = np.concatenate((x, x_ES))
            y = np.concatenate((y, y_ES))

        z = np.concatenate((z, z_ES))
        if fracs=='yes':
            fraction = np.concatenate((fraction, fraction_ES))
        if scaled:
            ES = np.concatenate((ES, ES_ES))

    if coords == 'proj':
        eta_scaled = Scale_feature(eta, Eta_Min, Eta_Max)
        phi_scaled = Scale_feature(phi, Phi_Min, Phi_Max)
        z_scaled = Scale_feature(z, Z_Min, Z_Max)
    elif coords == 'cart':
        x_scaled = Scale_feature(x, X_Min, X_Max)
        y_scaled = Scale_feature(y, Y_Min, Y_Max)
        z_scaled = Scale_feature(z, Z_Min, Z_Max)
    else:
        if np.abs(z[0]) < 300:
            iEta_scaled = Scale_feature(iEta, iEta_Min, iEta_Max)
            iPhi_scaled = Scale_feature(iPhi, iPhi_Min, iPhi_Max)
        else:
            iEta_scaled = Scale_feature(iEta, iX_Min, iX_Max)
            iPhi_scaled = Scale_feature(iPhi, iY_Min, iY_Max)

        iZ = np.sign(z)

    featlist = []
    if coords == 'proj':
        featlist += [eta_scaled, phi_scaled, z_scaled]
    elif coords == 'cart':
        featlist += [x_scaled, y_scaled, z_scaled]
    else:
        featlist += [iEta_scaled, iPhi_scaled, iZ]

    if fracs == 'yes':
        featlist += [energy_scaled, fraction]
    else:
        featlist += [energy_scaled]

    if scaled:
        featlist += [ES]

    train_feat = np.transpose(featlist)

    x = train_feat.astype(np.float32)

    return x

def extractRawE(event, i):
    if i == 1 and len(event.Hit_Eta_Ele1)==0:
        return event.Ele_SCRawE[0]
    else:
        return event.Ele_SCRawE[i]

def extractR9(event, i):
    if i == 1 and len(event.iEtaEle1)==0:
        return event.Ele_R9[0]
    else:
        return event.Ele_R9[i]

def extractBDT_track(event, i):
    if i == 1 and len(event.iEtaEle1)==0:
        return event.energy[0]/event.Ele_Gen_E[i]
    else:
        return event.energy[i]/event.Ele_Gen_E[i]

def extractMust(event, i):
    if i == 1 and len(event.iEtaEle1)==0:
        return event.energy_ecal_mustache[0]/event.Ele_Gen_E[i]
    else:
        return event.energy_ecal_mustache[i]/event.Ele_Gen_E[i]

def extractBDT(event, i): 
    if i == 1 and len(event.iEtaEle1)==0:
        return event.energy_ecal[0]/event.Ele_Gen_E[i]
    else:
        return event.energy_ecal[i]/event.Ele_Gen_E[i]

def extractHoE(event, i):
    if i == 1 and len(event.Hit_Eta_Ele1) == 0:
        return event.Ele_HadOverEm[0]
    else:
        return event.Ele_HadOverEm[i]

def extractNHit(event, i):
    if i==0:
        return len(event.iEtaEle1)
    else:
        return len(event.iEtaEle2)

def extractNHit_ES(event, i):
    if i==0:
        return len(event.Hit_ES_Eta_Ele1)
    else:
        return len(event.Hit_ES_Eta_Ele2)

def extractenergy(event, i):
    if i == 1 and len(event.iEtaEle1)==0:
        return event.energy[0]
    else:
        return event.energy[i]

def extractPt2(event, i):
    if i == 1 and len(event.pt)==1:
        return event.pt[0]
    else:
        return event.pt[i]

def extractenergy_ecal(event, i):
    if i == 1 and len(event.iEtaEle1)==0:
        return event.energy_ecal[0]
    else:
        return event.energy_ecal[i]
    
def extractenergy_ecal_mustache(event, i):
    if i == 1 and len(event.iEtaEle1)==0:
        return event.energy_ecal_mustache[0]
    else:
        return event.energy_ecal_mustache[i]

def extract_subdet(event, i):
    if i==0:
        if np.all(np.abs(event.Hit_Z_Ele1) < 300):
            return 0
        elif np.all(np.abs(event.Hit_Z_Ele1) > 300):
            return 1
        else:
            return 2
                
    else:
        if np.all(np.abs(event.Hit_Z_Ele2) < 300):
            return 0
        elif np.all(np.abs(event.Hit_Z_Ele2) > 300):
            return 1
        else:
            return 2

toTargets = {
    'trueE' : lambda event, i: event.Ele_Gen_E[i],
    'log' : lambda event, i: np.log(event.Ele_Gen_E[i]),
    'ratio' : lambda event, i: event.Ele_Gen_E[i]/extractRawE(event, i),
    'logratio' : lambda event, i: np.log(event.Ele_Gen_E[i]/extractRawE(event, i)),
    'ratiolog' : lambda event, i: np.log(event.Ele_Gen_E[i])/np.log(extractRawE(event, i)),
    'ratioflip' : lambda event, i: extractRawE(event, i)/event.Ele_Gen_E[i],
    'logratioflip' : lambda event, i: np.log(extractRawE(event, i)/event.Ele_Gen_E[i]),
    'ratiologflip' : lambda event, i: np.log(extractRawE(event, i))/np.log(event.Ele_Gen_E[i]),
    'logratiolog' : lambda event, i: np.log(toTargets['ratiolog'](event, i)),
    'logratiologflip' : lambda event, i: np.log(toTargets['ratiologflip'](event, i))
}

toRows = {
    'rho' : lambda event, i: event.rho,
    'phi' : lambda event, i: event.Ele_Gen_Phi[i],
    'eta' : lambda event, i: event.Ele_Gen_Eta[i],
    'pt' : lambda event, i: event.Ele_Gen_Pt[i],
    'pt2' : lambda event, i: extractPt2(event, i),
    'trueE' : lambda event, i: event.Ele_Gen_E[i],
    'rawE' : lambda event, i: extractRawE(event, i),
    'R9' : lambda event, i: extractR9(event, i),
    'BDT': lambda event, i: extractBDT(event, i),
    'BDT_track': lambda event, i: extractBDT_track(event, i),
    'nHit' : lambda event, i: extractNHit(event, i),
    'nHit_ES' : lambda event, i: extractNHit_ES(event, i),
    'HoE' : lambda event, i: extractHoE(event, i),
    'must' : lambda event, i: extractMust(event, i),

    'energy' : lambda event, i: extractenergy(event,i),
    'energy_ecal' : lambda event, i: extractenergy_ecal(event,i),
    'energy_ecal_mustache' : lambda event, i: extractenergy_ecal_mustache(event,i),

    'etareco' : lambda event, i: event.eta[i],
    'phireco' : lambda event, i: event.phi[i],

    'flag' : lambda event, i: len(event.energy),

    'data_f1': lambda event, i : extract_data_v1_f1(event, i),
    'data_f2': lambda event, i : extract_data_v1_f2(event, i),

    'iEta' : lambda event, i : np.asarray(event.iEtaEle1,dtype=int) if i==0 else np.asarray(event.iEtaEle2,dtype=int),
    'iPhi' : lambda event, i : np.asarray(event.iPhiEle1,dtype=int) if i==0 else np.asarray(event.iPhiEle2,dtype=int),
    'iZ' : lambda event, i : np.sign(np.asarray(event.Hit_Z_Ele1) if i==0 else np.asarray(event.Hit_Z_Ele2)),

    'subdet' : lambda event, i : extract_subdet(event, i),
    'dR' : None
}


###########################################################################

########### Cuts ##########################################################

def cut_Zee(event, i):
    return event.nElectrons==2

def cut_gun(event, i):
    return event.Ele_Gen_E[i] > 5\
            and event.Ele_Gen_E[i] < 300 \
            and (event.Hit_Eta_Ele1.size()>0 if i==0 else event.Hit_Eta_Ele2.size()>0) 

cuts = {
    "Zee": cut_Zee,
    "gun": cut_gun
}

###########################################################################

################ Extraction ###############################################
def getTree(f = "nTuples/2018_UL_Particle_MC_v3/nTupleMC_Merged.root",T=None):
    f = ROOT.TFile(f)
    if T is None:
        T = f.Get('nTuplelize/T')
    else: 
        T = f.Get(T)
    return f, T

def extract(name, varname, N, cut, toRow, folder, f, T=None):
    f, T = getTree(f,T)

    print("Extracting %s..."%varname)

    i=0
    dataset=[]
    for event in tqdm(T):
        i=i+1
        if N>0 and i==N:
            break
        
        if cut(event, 0):
            dataset.append(toRow(event, 0))

        if cut(event, 1):
            dataset.append(toRow(event, 1))

    print("\tDone.")
    print("\tThere are",len(dataset),"entries")
    print("\tAn entry looks like",dataset[0])
    print()

    print("\tDumping..")
    t0 = time()
    filename = '%s/%s_%s.pickle'%(folder, varname, name)
    with open(filename, 'wb') as f:
        pickle.dump(dataset, f)
    t1 = time()
    print("\tDone. Took %.3f seconds"%(t1-t0))
    print()
    print()

def match(event, i):
    minDR = 999
    idx = -1

    phiTs = event.Ele_Gen_Phi
    etaTs = event.Ele_Gen_Eta

    if len(phiTs) < 2:
        return idx, minDR

    phi = toRows['phireco'](event, i)
    eta = toRows['etareco'](event, i)

    for j in range(len(phiTs)):
        dphi = phi - phiTs[j]
        if np.abs(dphi) <= np.pi:
            dphi = np.abs(dphi)
        else:
            dphi = 2*np.pi - np.abs(dphi)
        deta = np.abs(eta - etaTs[j])

        dR = np.sqrt(dphi*dphi + deta*deta)

        if dR < minDR:
            minDR = dR
            idx = j

    return idx, minDR

def extractGM(name, varname, N, cut, toRow, folder, f, T=None):
    f, T = getTree(f,T)

    print("Extracting genmatched %s..."%varname)

    i=0
    dataset=[]
    for event in tqdm(T):
        i=i+1
        if N>0 and i==N:
            break
        
        if cut(event, 0):
            idx, dR = match(event, 0)
            if varname != 'dR' and idx>=0:
                dataset.append(toRow(event, idx))
            elif idx>=0:
                dataset.append(dR)
            else:
                dataset.append(999)

        if cut(event, 1):
            idx, dR = match(event, 1)
            if varname != 'dR' and idx>=0:
                dataset.append(toRow(event, idx))
            elif idx>=0:
                dataset.append(dR)
            else:
                dataset.append(999)

    print("\tDone.")
    print("\tThere are",len(dataset),"entries")
    print("\tAn entry looks like",dataset[0])
    print()

    print("\tDumping..")
    t0 = time()
    filename = '%s/%s_%s.pickle'%(folder, varname, name)
    with open(filename, 'wb') as f:
        pickle.dump(dataset, f)
    t1 = time()
    print("\tDone. Took %.3f seconds"%(t1-t0))
    print()
    print()

def extractFeatures(name, ES, coords, fracs, N, cut, folder, 
        f = 'nTuples/2018_UL_Particle_MC_v3/nTupleMC_Merged.root'):

    cut = cuts[cut]
    varname = 'features_%sES_%s_%sfrac'%(ES,coords,fracs)
    toRow = lambda event, i: ElectronToFeatures(event, i, ES, coords, fracs)

    extract(name, varname, N, cut, toRow, folder, f)

def extractTargets(name, target, N, cut, folder,
        f = 'nTuples/2018_UL_Particle_MC_v3/nTupleMC_Merged.root'):

    cut = cuts[cut]
    varname = 'targets_%s'%target
    toRow = toTargets[target]

    extract(name, varname, N, cut, toRow, folder, f)

def extractVariables(name, varname, N, cut, folder, 
        f = 'nTuples/2018_UL_Particle_MC_v3/nTupleMC_Merged.root', T = None):

    cut = cuts[cut]
    toRow = toRows[varname]

    extract(name, varname, N, cut, toRow, folder, f, T)

def extractGenMatched(name, varname, N, cut, folder, 
        f = 'nTuples/2018_UL_Particle_MC_v3/nTupleMC_Merged.root', T = None):

    cut = cuts[cut]
    toRow = toRows[varname]

    extractGM(name, varname, N, cut, toRow, folder, f, T)

def deltaR(eta1, eta2, phi1, phi2):
    deta = eta1-eta2

    dphi = phi1-phi2 if np.abs(phi1-phi2) < np.pi else 2*np.pi - np.abs(phi1-phi2)

    return np.sqrt(deta*deta + dphi*dphi)

###########################################################################

import sys

def test(fname, N):
    sys.stdout = open(fname, 'w')

    #f = "nTuples/2018_UL_Particle_MC_v4/nTupleMC_Merged.root"
    f = ROOT.TFile("nTuples/2018_UL_Particle_MC_v4/nTupleMC_Merged.root")
    T = f.Get('nTuplelize/T')
    i = 0
    print("trueE","energy","energy_ecal","energy_ecal_mustache","eta","flag")
    for event in T:
        i+=1
        if i>=N:
            sys.stdout.close()
            return

        if len(event.energy) != len(event.energy_ecal) or len(event.energy) != len(event.energy_ecal_mustache) or len(event.energy_ecal) != len(event.energy_ecal_mustache):
            print("\n\n\n!!!!! BAD BAD BAD \n\n\n")
        if len(event.energy)==2:
            flag = 0
            for ele in range(2):
                if event.Ele_Gen_Eta[ele] <1.4442:
                    print(event.Ele_Gen_E[ele], 
                            event.energy[ele], 
                            event.energy_ecal[ele], 
                            event.energy_ecal_mustache[ele],
                            event.Ele_Gen_Eta[ele],
                            flag)
        elif len(event.energy)==1:
            flag = 1
            if len(event.Hit_X_Ele1)==0:
                ele = 1
            else:
                ele = 0

            print(event.Ele_Gen_E[ele], 
                    event.energy[0], 
                    event.energy_ecal[0], 
                    event.energy_ecal_mustache[0],
                    event.Ele_Gen_Eta[ele],
                    flag)

    sys.stdout.close()
