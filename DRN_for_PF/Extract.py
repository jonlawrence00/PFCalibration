import uproot 
#from numba import jit
import numpy as np
import awkward as ak
from time import time
import pickle
import tqdm
from torch_geometric.data import Data
import torch

MISSING = -999999

########################################
# HGCAL Values                         #
########################################
##### Full
#HGCAL_X_Min =-261.7752
#HGCAL_X_Max = 261.77518
#HGCAL_Y_Min = -261.77524
#HGCAL_Y_Max = 261.77518
#HGCAL_Z_Min = -522.798
#HGCAL_Z_Max = 522.798
#HGCAL_Min = 0.1 #0.10000002
#HGCAL_Max = 425.68573#425.68573 #431.37378

##### trans_pos
HGCAL_X_Min =-261.77518
HGCAL_X_Max = 261.77518
HGCAL_Y_Min = -261.77518
HGCAL_Y_Max = 261.77518
HGCAL_Z_Min = 134.62964
HGCAL_Z_Max = 522.798
HGCAL_Min = 0.10000003
HGCAL_Max = 359.28555

##### Full pos
#HGCAL_X_Min = -261.77518
#HGCAL_X_Max = 261.77518
#HGCAL_Y_Min = -261.77518
#HGCAL_Y_Max = 261.77518
#HGCAL_Z_Min = -32.03129
#HGCAL_Z_Max = 522.798
#HGCAL_Min = 0.10000003
#HGCAL_Max = 388.53964

##### ec_out pos
#HGCAL_X_Min = -94.21816
#HGCAL_X_Max = 94.21816
#HGCAL_Y_Min = -94.218155
#HGCAL_Y_Max = 94.218155
#HGCAL_Z_Min = 317.5769
#HGCAL_Z_Max = 522.798
#HGCAL_Min = 0.10000199
#HGCAL_Max = 271.34967

##### ec_in pos 
#HGCAL_X_Min = -232.83936
#HGCAL_X_Max = 232.83936
#HGCAL_Y_Min = -232.83936
#HGCAL_Y_Max = 232.83936
#HGCAL_Z_Min = 265.24896
#HGCAL_Z_Max = 522.798
#HGCAL_Min = 0.10000035
#HGCAL_Max = 254.90266


#HGCAL_Max = 350.79248
# HGCAL_Max_EE =3200
# HGCAL_Max_FH = 2200
# HGCAL_MAX_AH =800
Eta_Max = 2.5
Eta_Min = -2.5
#Eta_Min = -1.549997
#Eta_Max = 1.549997


########################################
# ECAL Values                          #
########################################

X_Min = -150
X_Max = 150

Y_Min = -150
Y_Max = 150

Z_Min = -330
Z_Max = 330

Phi_Min = -np.pi
Phi_Max = np.pi

iEta_Min = -85
iEta_Max = 85

iPhi_Min = 1
iPhi_Max = 360

iX_Min = 1
iX_Max = 100

iY_Min = 1
iY_Max = 100

ECAL_Min = 0
ECAL_Max = 250

def rescale(feature, minval, maxval):
    top = feature-minval
    bot = maxval-minval
    return top/bot

def dphi(phi1, phi2):
    dphi = np.abs(phi1-phi2)
    gt = dphi > np.pi
    dphi[gt] = 2*np.pi - dphi[gt]
    return dphi

def dR(eta1, eta2, phi1, phi2):
    dp = dphi(phi1, phi2)
    de = np.abs(eta1-eta2)

    return np.sqrt(dp*dp + de*de)

def cartfeat_HGCAL(x,y,z, En):
    #frac =  ((z<54)*0.0105) + (np.logical_and(z>54, z<154)*0.0789) + ((z>154)*0.0316)
    #((z<54)*0.035) + ((z>54)*0.095) #(np.logical_and(z>54, z<154)*0.0789) + ((z>154)*0.0316)
    #    En = En*frac
    # if(z<54):
    #     HGCAL_Max=3200
    # elif(z>54 and z<154):
    #     HGCAL_Max=2200
    # elif(z>154):
    #     HGCAL_Max=850
    #HGCAL_Max = (z<54)*3300 +  (np.logical_and(z>54, z<154)*2500) + ((z>154)*900)
    E = rescale(En, HGCAL_Min, HGCAL_Max)
    x = rescale(x, HGCAL_X_Min, HGCAL_X_Max)
    y = rescale(y, HGCAL_Y_Min, HGCAL_Y_Max)
    z = rescale(z, HGCAL_Z_Min, HGCAL_Z_Max)
    #eta = rescale(eta, Eta_Min, Eta_Max)
    return ak.concatenate((x[:,:,None], y[:,:,None], z[:,:,None], E[:,:,None]), -1)
    #return ak.concatenate((z[:,:,None], E[:,:,None]), -1)

def cartfeat(x, y, z, En ,frac, det=None):
    E = rescale(En*frac, ECAL_Min, ECAL_Max)
    x = rescale(x, X_Min, X_Max)
    y = rescale(y, Y_Min, Y_Max)
    z = rescale(z, Z_Min, Z_Max)

    if det is None:
        return ak.concatenate((x[:,:,None], y[:,:,None], z[:,:,None], E[:,:,None]), -1)
    else:
        return ak.concatenate((x[:,:,None], y[:,:,None], z[:,:,None], E[:,:,None], det[:,:,None]), -1)

def projfeat(eta, phi, z, En ,frac, det=None):
    E = rescale(En*frac, ECAL_Min, ECAL_Max)
    eta = rescale(eta, Eta_Min, Eta_Max)
    phi = rescale(phi, Phi_Min, Phi_Max)
    z = rescale(z, Z_Min, Z_Max)

    if det is None:
        return ak.concatenate((eta[:,:,None], phi[:,:,None], z[:,:,None], E[:,:,None]), -1)
    else:
        return ak.concatenate((eta[:,:,None], phi[:,:,None], z[:,:,None], E[:,:,None], det[:,:,None]), -1)

def localfeat(i1, i2, z, En ,frac, det=None):
    '''
    In the barrel:
        i1 = iEta
        i2 = iPhi
    In the endcaps:
        i1 = iX
        i2 = iY
    '''

    if det is not None:
        print("Error: local coordinates not defined for ES")
        return

    E = rescale(En*frac, ECAL_Min, ECAL_Max)

    Zfirst = ak.firsts(z)
    barrel = np.abs(Zfirst) < 300 #this is 1 if we are in the barrel, 0 in the endcaps
    
    xmax = barrel * iEta_Max + ~barrel * iX_Max
    xmin = barrel * iEta_Min + ~barrel * iX_Min

    ymax = barrel * iPhi_Max + ~barrel * iY_Max
    ymin = barrel * iPhi_Min + ~barrel * iY_Min

    x = rescale(i1, xmin, xmax)
    y = rescale(i2, ymin, ymax)

    whichEE = 2*(Zfirst > 300) - 1 #+1 to the right of 0, -1 to the left of 0

    iZ = whichEE * ~barrel #0 in the barrel, -1 in left EE, +1 in right EE

    iZ, _ = ak.broadcast_arrays(iZ, x)

    return ak.concatenate((x[:,:,None], y[:,:,None], iZ[:,:,None], E[:,:,None]), -1)

def torchify(feat, graph_x = None):
    data = [Data(x = torch.from_numpy(ak.to_numpy(ele).astype(np.float32))) for ele in feat]
    if graph_x is not None:
        for d, gx in zip(data, graph_x):
            d.graph_x = gx
    return data

def npify(feat):
    t0 = time()
    data = [ak.to_numpy(ele) for ele in feat]
    print("took %f"%(time()-t0))
    return data

varlists = {
    'BDTvars': ['Pho_R9', #'Pho_S4', S4 is not populated
                 'Pho_SigIEIE', 'Pho_SigIPhiIPhi',
                 'Pho_SCEtaW', 'Pho_SCPhiW',
                 #'Pho_CovIEtaIEta', 'Pho_CovIEtaIPhi','Pho_ESSigRR', not populated
                 'Pho_SCRawE', 
                 'Pho_SC_ESEnByRawE', 'Pho_HadOverEm',
                 'eta', 'phi',
                 'Pho_Gen_Eta', 'Pho_Gen_Phi',
                 'iEtaPho1', 'iEtaPho2', 'Hit_Z_Pho1', 'Hit_Z_Pho2', "Pho_Gen_E"],
    'gun_pho': ['nPhotons', 
                'Pho_Gen_E', 'Pho_Gen_Eta', 'Pho_Gen_Phi', 
                'Pho_SCRawE', 'pt', 'eta', 'phi',
                'Pho_R9', 'Pho_HadOverEm', 'rho',
                'iEtaPho1', 'iEtaPho2',
                'iPhiPho1', 'iPhiPho2',
                'Hit_ES_Eta_Pho1', 'Hit_ES_Eta_Pho2',
                'Hit_ES_Phi_Pho1', 'Hit_ES_Phi_Pho2',
                'Hit_ES_X_Pho1', 'Hit_ES_X_Pho2',
                'Hit_ES_Y_Pho1', 'Hit_ES_Y_Pho2',
                'Hit_ES_Z_Pho1', 'Hit_ES_Z_Pho2',
                'ES_RecHitEnPho1', 'ES_RecHitEnPho2',
                'Hit_Eta_Pho1', 'Hit_Eta_Pho2',
                'Hit_Phi_Pho1', 'Hit_Phi_Pho2',
                'Hit_X_Pho1', 'Hit_X_Pho2',
                'Hit_Y_Pho1', 'Hit_Y_Pho2',
                'Hit_Z_Pho1', 'Hit_Z_Pho2',
                'RecHitEnPho1', 'RecHitEnPho2',
                'RecHitFracPho1', 'RecHitFracPho2',
                #'passLooseId', 'passMediumId', 'passTightId',
                'energy'],
    'Hgg': ['nPhotons', 
                'Pho_SCRawE', 'eta', 'phi',
                'Pho_R9', 'Pho_HadOverEm', 'rho',
                'iEtaPho1', 'iEtaPho2',
                'iPhiPho1', 'iPhiPho2',
                'Hit_ES_Eta_Pho1', 'Hit_ES_Eta_Pho2',
                'Hit_ES_Phi_Pho1', 'Hit_ES_Phi_Pho2',
                'Hit_ES_X_Pho1', 'Hit_ES_X_Pho2',
                'Hit_ES_Y_Pho1', 'Hit_ES_Y_Pho2',
                'Hit_ES_Z_Pho1', 'Hit_ES_Z_Pho2',
                'ES_RecHitEnPho1', 'ES_RecHitEnPho2',
                'Hit_Eta_Pho1', 'Hit_Eta_Pho2',
                'Hit_Phi_Pho1', 'Hit_Phi_Pho2',
                'Hit_X_Pho1', 'Hit_X_Pho2',
                'Hit_Y_Pho1', 'Hit_Y_Pho2',
                'Hit_Z_Pho1', 'Hit_Z_Pho2',
                'RecHitEnPho1', 'RecHitEnPho2',
                'RecHitFracPho1', 'RecHitFracPho2',
                #'passLooseId', 'passMediumId', 'passTightId',
                'energy'],

    'gun_30M': ['nElectrons', 
                'Ele_Gen_E', 'Ele_Gen_Eta', 'Ele_Gen_Phi', 
                'Ele_SCRawE', 'eta', 'phi',
                'Ele_R9', 'Ele_HadOverEm', 'rho',
                'iEtaEle1', 'iEtaEle2',
                'iPhiEle1', 'iPhiEle2',
                'Hit_ES_Eta_Ele1', 'Hit_ES_Eta_Ele2',
                'Hit_ES_Phi_Ele1', 'Hit_ES_Phi_Ele2',
                'Hit_ES_X_Ele1', 'Hit_ES_X_Ele2',
                'Hit_ES_Y_Ele1', 'Hit_ES_Y_Ele2',
                'Hit_ES_Z_Ele1', 'Hit_ES_Z_Ele2',
                'ES_RecHitEnEle1', 'ES_RecHitEnEle2',
                'Hit_Eta_Ele1', 'Hit_Eta_Ele2',
                'Hit_Phi_Ele1', 'Hit_Phi_Ele2',
                'Hit_X_Ele1', 'Hit_X_Ele2',
                'Hit_Y_Ele1', 'Hit_Y_Ele2',
                'Hit_Z_Ele1', 'Hit_Z_Ele2',
                'RecHitEnEle1', 'RecHitEnEle2',
                'RecHitFracEle1', 'RecHitFracEle2',
                'passLooseId', 'passMediumId', 'passTightId',
                'energy_ecal_mustache'],
    'gun_v3': ['nElectrons', 
                'Ele_Gen_E', 'Ele_Gen_Eta', 'Ele_Gen_Phi', 
                'Ele_SCRawE', 'eta', 'phi',
                'Ele_R9', 'Ele_HadOverEm', 'rho',
                'iEtaEle1', 'iEtaEle2',
                'iPhiEle1', 'iPhiEle2',
                'Hit_ES_Eta_Ele1', 'Hit_ES_Eta_Ele2',
                'Hit_ES_Phi_Ele1', 'Hit_ES_Phi_Ele2',
                'Hit_ES_X_Ele1', 'Hit_ES_X_Ele2',
                'Hit_ES_Y_Ele1', 'Hit_ES_Y_Ele2',
                'Hit_ES_Z_Ele1', 'Hit_ES_Z_Ele2',
                'ES_RecHitEnEle1', 'ES_RecHitEnEle2',
                'Hit_Eta_Ele1', 'Hit_Eta_Ele2',
                'Hit_Phi_Ele1', 'Hit_Phi_Ele2',
                'Hit_X_Ele1', 'Hit_X_Ele2',
                'Hit_Y_Ele1', 'Hit_Y_Ele2',
                'Hit_Z_Ele1', 'Hit_Z_Ele2',
                'RecHitEnEle1', 'RecHitEnEle2',
                'RecHitFracEle1', 'RecHitFracEle2'],
    'Zee_data': ['nElectrons', 
                'Ele_SCRawE', 'eta', 'phi',
                'Ele_R9', 'Ele_HadOverEm', 'rho',
                'iEtaEle1', 'iEtaEle2',
                'iPhiEle1', 'iPhiEle2',
                'Hit_ES_Eta_Ele1', 'Hit_ES_Eta_Ele2',
                'Hit_ES_Phi_Ele1', 'Hit_ES_Phi_Ele2',
                'Hit_ES_X_Ele1', 'Hit_ES_X_Ele2',
                'Hit_ES_Y_Ele1', 'Hit_ES_Y_Ele2',
                'Hit_ES_Z_Ele1', 'Hit_ES_Z_Ele2',
                'ES_RecHitEnEle1', 'ES_RecHitEnEle2',
                'Hit_Eta_Ele1', 'Hit_Eta_Ele2',
                'Hit_Phi_Ele1', 'Hit_Phi_Ele2',
                'Hit_X_Ele1', 'Hit_X_Ele2',
                'Hit_Y_Ele1', 'Hit_Y_Ele2',
                'Hit_Z_Ele1', 'Hit_Z_Ele2',
                'RecHitEnEle1', 'RecHitEnEle2',
                'RecHitFracEle1', 'RecHitFracEle2',
                'energy_ecal_mustache'],
    'Zee_MC' : ['nElectrons', 
                'Ele_Gen_E', 'Ele_Gen_Eta', 'Ele_Gen_Phi', 
                'Ele_SCRawE', 'eta', 'phi',
                'Ele_R9', 'Ele_HadOverEm', 'rho',
                'iEtaEle1', 'iEtaEle2',
                'iPhiEle1', 'iPhiEle2',
                'Hit_ES_Eta_Ele1', 'Hit_ES_Eta_Ele2',
                'Hit_ES_Phi_Ele1', 'Hit_ES_Phi_Ele2',
                'Hit_ES_X_Ele1', 'Hit_ES_X_Ele2',
                'Hit_ES_Y_Ele1', 'Hit_ES_Y_Ele2',
                'Hit_ES_Z_Ele1', 'Hit_ES_Z_Ele2',
                'ES_RecHitEnEle1', 'ES_RecHitEnEle2',
                'Hit_Eta_Ele1', 'Hit_Eta_Ele2',
                'Hit_Phi_Ele1', 'Hit_Phi_Ele2',
                'Hit_X_Ele1', 'Hit_X_Ele2',
                'Hit_Y_Ele1', 'Hit_Y_Ele2',
                'Hit_Z_Ele1', 'Hit_Z_Ele2',
                'RecHitEnEle1', 'RecHitEnEle2',
                'RecHitFracEle1', 'RecHitFracEle2',
                'energy_ecal_mustache'],

} 

gun_readcut = 'nElectrons>0'
gun_pho_readcut = 'nPhotons>0'
Zee_readcut = 'nElectrons==2'

readcuts = {
    'gun_30M': gun_readcut,
    'gun_v3': gun_readcut,

    'Zee_data': Zee_readcut,
    'Zee_MC': gun_readcut,

    'gun_pho' : gun_pho_readcut,

    'BDTvars' : gun_pho_readcut,

    'Hgg' : 'nPhotons==2',
}

def gun_savecut(result):
    return np.logical_and(result['Ele_Gen_E'] < 300, result['Ele_Gen_E'] > 5)

def gun_pho_savecut(result):
    return np.logical_and(result['Pho_Gen_E'] < 300, result['Pho_Gen_E'] > 5)

def Zee_savecut(result):
    return np.ones(result['phi'].shape, dtype=bool)


savecuts = {
    'gun_30M': gun_savecut,
    'gun_v3' : gun_savecut,

    'Zee_data': Zee_savecut,
    'Zee_MC': Zee_savecut,
    
    'gun_pho' : gun_pho_savecut,
    'BDTvars' : gun_pho_savecut,

    'Hgg' : Zee_savecut,
}

hasgen = {
    'gun_30M' : True,
    'gun_v3' : True,

    'Zee_data' : False,
    'Zee_MC' : True,

    'gun_pho' : True,
    'BDTvars' : True,

    'Hgg': False,
}

isEle = {
    'gun_30M' : True,
    'gun_v3' : True,

    'Zee_data': True,
    'Zee_MC': True,

    'gun_pho' : False,
    'BDTvars' : False,

    'Hgg' : False,
}

class Extract:
    def __init__(self, outfolder, path, treeName='nTuplelize/T'):
        if path is not None:
            #path = '~/shared/nTuples/%s'%path
            self.tree = uproot.open("%s:%s"%(path, treeName))

        self.outfolder = outfolder

    def get_subdet(self):
        print("Getting subdet")

        t0 = time()
        with open("%s/Hit_Z.pickle"%self.outfolder, 'rb') as f:
            Z = pickle.load(f)
        print("\tLoaded Hit_Z in %0.2f seconds"%(time()-t0))
        t0 = time()
    
        subdet = np.abs(ak.to_numpy(ak.firsts(Z))) < 300

        print("dumping...")
        with open("%s/subdet.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(subdet,f)
        print('done')

    def build_localfeat(self, ES=False, scaled=False):
        if ES:
            print("Error: no local coords for ES")
            return 

        print("Building localfeat")
        t0 = time()
        with open("%s/iEta.pickle"%self.outfolder, 'rb') as f:
            iEta = pickle.load(f)
        print("\tLoaded iEta in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/iPhi.pickle"%self.outfolder, 'rb') as f:
            iPhi = pickle.load(f)
        print("\tLoaded iPhi in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/Hit_Z.pickle"%self.outfolder, 'rb') as f:
            Z = pickle.load(f)
        print("\tLoaded Hit_Z in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/RecHitEn.pickle"%self.outfolder, 'rb') as f:
            En = pickle.load(f)
        print("\tLoaded En in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/RecHitFrac.pickle"%self.outfolder, 'rb') as f:
            frac = pickle.load(f)
        print("\tLoaded Frac in %0.2f seconds"%(time()-t0))
        t0 = time()

        if ES:
            with open("%s/iEta.pickle"%self.outfolder, 'rb') as f:
                iEta = pickle.load(f)
            print("\tLoaded iEta in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/iPhi.pickle"%self.outfolder, 'rb') as f:
                iPhi = pickle.load(f)
            print("\tLoaded iPhi in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/Hit_Z.pickle"%self.outfolder, 'rb') as f:
                Z = pickle.load(f)
            print("\tLoaded Hit_Z in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/RecHitEn.pickle"%self.outfolder, 'rb') as f:
                En = pickle.load(f)
            print("\tLoaded En in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/RecHitFrac.pickle"%self.outfolder, 'rb') as f:
                frac = pickle.load(f)
            print("\tLoaded Frac in %0.2f seconds"%(time()-t0))
            t0 = time()


    
        lf = localfeat(iEta, iPhi, Z, En, frac)
        print("\tMake localfeat in %0.2f seconds"%(time()-t0))
        t0 = time()

        lf = torchify(lf)
        print("\tTorchified in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/localfeat.pickle"%self.outfolder, 'wb') as f:
            torch.save(lf, f, pickle_protocol=4)
        print("\tDumped in %0.2f seconds"%(time()-t0))

    def build_projfeat(self, ES=False, scaled=False):
        print("Building projfeat")
        t0 = time()
        with open("%s/Hit_Eta.pickle"%self.outfolder, 'rb') as f:
            Eta = pickle.load(f)
        print("\tLoaded Eta in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/Hit_Phi.pickle"%self.outfolder, 'rb') as f:
            Phi = pickle.load(f)
        print("\tLoaded Phi in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/Hit_Z.pickle"%self.outfolder, 'rb') as f:
            Z = pickle.load(f)
        print("\tLoaded Hit_Z in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/RecHitEn.pickle"%self.outfolder, 'rb') as f:
            En = pickle.load(f)
        print("\tLoaded En in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/RecHitFrac.pickle"%self.outfolder, 'rb') as f:
            frac = pickle.load(f)
        print("\tLoaded Frac in %0.2f seconds"%(time()-t0))
        t0 = time()

        if ES:
            with open("%s/Hit_ES_Eta.pickle"%self.outfolder, 'rb') as f:
                ES_Eta = pickle.load(f)
            print("\tLoaded ES_Eta in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/Hit_ES_Phi.pickle"%self.outfolder, 'rb') as f:
                ES_Phi = pickle.load(f)
            print("\tLoaded ES_Phi in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/Hit_ES_Z.pickle"%self.outfolder, 'rb') as f:
                ES_Z = pickle.load(f)
            print("\tLoaded ES_Z in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/ES_RecHitEn.pickle"%self.outfolder, 'rb') as f:
                ES_En = pickle.load(f)
            print("\tLoaded ES_En in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/ES_RecHitFrac.pickle"%self.outfolder, 'rb') as f:
                ES_frac = pickle.load(f)
            print("\tLoaded ES_Frac in %0.2f seconds"%(time()-t0))
            t0 = time()

            if scaled:
                ES_En = ES_En*3500
                fname = 'projfeat_ES_scaled'
            else:
                fname = 'projfeat_ES'

            ES = ak.ones_like(ES_Eta)
            ECAL = ak.ones_like(Eta)

            Eta = ak.concatenate( (Eta, ES_Eta), axis=1)
            Phi = ak.concatenate( (Phi, ES_Phi), axis=1)
            Z = ak.concatenate( (Z, ES_Z), axis=1)
            En = ak.concatenate( (En, ES_En), axis=1)
            frac = ak.concatenate( (frac, ES_frac), axis=1)
            det = ak.concatenate( (ECAL, ES), axis=1)
        else:
            fname = 'projfeat'
            det = None

        pf = projfeat(Eta, Phi, Z, En, frac, det)
        print("\tMake projfeat in %0.2f seconds"%(time()-t0))
        t0 = time()

        lf = torchify(lf)
        print("\tTorchified in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/%s.pickle"%(self.outfolder, fname), 'wb') as f:
            torch.save(pf, f, pickle_protocol=4)
        print("\tDumped in %0.2f seconds"%(time()-t0))

    def build_cartfeat(self, ES=False, scaled=False, graph_features = None):
        print("Building cartfeat")
        t0 = time()
        with open("%s/Hit_X.pickle"%self.outfolder, 'rb') as f:
            X = pickle.load(f)
        print("\tLoaded X in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/Hit_Y.pickle"%self.outfolder, 'rb') as f:
            Y = pickle.load(f)
        print("\tLoaded Y in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/Hit_Z.pickle"%self.outfolder, 'rb') as f:
            Z = pickle.load(f)
        print("\tLoaded Hit_Z in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/RecHitEn.pickle"%self.outfolder, 'rb') as f:
            En = pickle.load(f)
        print("\tLoaded En in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/RecHitFrac.pickle"%self.outfolder, 'rb') as f:
            frac = pickle.load(f)
        print("\tLoaded Frac in %0.2f seconds"%(time()-t0))
        t0 = time()

        if ES:
            with open("%s/Hit_ES_X.pickle"%self.outfolder, 'rb') as f:
                ES_X = pickle.load(f)
            print("\tLoaded ES_X in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/Hit_ES_Y.pickle"%self.outfolder, 'rb') as f:
                ES_Y = pickle.load(f)
            print("\tLoaded ES_Y in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/Hit_ES_Z.pickle"%self.outfolder, 'rb') as f:
                ES_Z = pickle.load(f)
            print("\tLoaded ES_Z in %0.2f seconds"%(time()-t0))
            t0 = time()

            with open("%s/ES_RecHitEn.pickle"%self.outfolder, 'rb') as f:
                ES_En = pickle.load(f)
            print("\tLoaded ES_En in %0.2f seconds"%(time()-t0))
            t0 = time()

            ES_frac = ak.ones_like(ES_En)

            if scaled:
                ES_En = ES_En*3500
                fname = 'cartfeat_ES_scaled'
            else:
                fname = 'cartfeat_ES'

            ES = ak.ones_like(ES_En)
            ECAL = ak.ones_like(En)

            X = ak.concatenate( (X, ES_X), axis=1)
            Y = ak.concatenate( (Y, ES_Y), axis=1)
            Z = ak.concatenate( (Z, ES_Z), axis=1)
            En = ak.concatenate( (En, ES_En), axis=1)
            frac = ak.concatenate( (frac, ES_frac), axis=1)
            det = ak.concatenate( (ECAL, ES), axis=1)
        else:
            det = None
            fname = 'cartfeat'

        graph_x = None
        if graph_features is not None:
            graph_x = []
            for var in graph_features:
                t0 = time()
                with open("%s/%s.pickle"%(self.outfolder, var), 'rb') as f:
                    graph_x.append(pickle.load(f))
                print("\tLoaded %s in %0.2f seconds"%(var, time()-t0))
                fname += "_%s"%var
            graph_x = np.concatenate(graph_x, 1)


        cf = cartfeat(X, Y, Z, En, frac, det)
        print("\tMake cartfeat in %0.2f seconds"%(time()-t0))
        t0 = time()

        cf = torchify(cf, graph_x)
        print("\tTorchified in %0.2f seconds"%(time()-t0))
        t0 = time()

        with open("%s/%s.pickle"%(self.outfolder, fname), 'wb') as f:
            torch.save(cf, f, pickle_protocol=4)
        print("\tDumped in %0.2f seconds"%(time()-t0))

    def add_graph_features(self, coords, ES, scaled, graph_features):
        if type(graph_features)==str:
            graph_features = [graph_features]

        fname = "%s/%sfeat"%(self.outfolder, coords)
        if ES:
            fname += "_ES"
            if scaled:
                fname += "_scaled"

        print("Adding features",graph_features,"to",fname)

        graph_x = []
        suffix = ''
        for var in graph_features:
            t0 = time()
            with open("%s/%s.pickle"%(self.outfolder, var), 'rb') as f:
                graph_x += [pickle.load(f)]
            print("\tLoaded %s in %0.2f seconds"%(var, time()-t0))
            suffix += "_%s"%var
        if len(graph_x)==1:
            graph_x = graph_x[0]
        else:
            graph_x = np.stack(graph_x, 1)

        t0 = time()
        data = torch.load("%s.pickle"%fname)
        print("\tLoaded node features in %0.2f seconds"%(time()-t0))
        
        for d, gx in zip(data, graph_x):
            d.graph_x = torch.tensor(gx)

        fname+=suffix

        t0 = time()
        with open("%s.pickle"%fname, 'wb') as f:
            torch.save(data, f, pickle_protocol=4)
        print("\tDumped in %0.2f seconds"%(time()-t0))

        return data

    @staticmethod
    def gen_match(phigen, phireco, etagen, etareco, threshold= 0.05):
        idxs = ak.argcartesian( (etagen, etareco), axis = 1)

        #awkward insists that I index the cartesitna product pairs with '0' and '1' rather than ints
        genetas = etagen[idxs[:,:,'0']]
        recoetas = etareco[idxs[:,:,'1']]

        genphis = phigen[idxs[:,:,'0']]
        recophis = phireco[idxs[:,:,'1']]

        dphis = np.abs(genphis - recophis)
        gt = dphis > np.pi
        #you can't assign to awkward arrays in place, so this is an inefficient hack
        dphis = gt * (2*np.pi - dphis) + (1 - gt)*(dphis) 

        detas = np.abs(genetas - recoetas)

        dR2s = dphis*dphis + detas*detas

        matched = dR2s < threshold*threshold

        return idxs[matched] #these are (gen, reco) index pairs

    @staticmethod
    def gen_unmatch(phigen, phireco, etagen, etareco, threshold=0.4):
        idxs = ak.argcartesian( (etagen, etareco), axis = 1)

        #awkward insists that I index the cartesitna product pairs with '0' and '1' rather than ints
        genetas = etagen[idxs[:,:,'0']]
        recoetas = etareco[idxs[:,:,'1']]

        genphis = phigen[idxs[:,:,'0']]
        recophis = phireco[idxs[:,:,'1']]

        dphis = np.abs(genphis - recophis)
        gt = dphis > np.pi
        #you can't assign to awkward arrays in place, so this is an inefficient hack
        dphis = gt * (2*np.pi - dphis) + (1 - gt)*(dphis) 

        detas = np.abs(genetas - recoetas)

        dR2s = dphis*dphis + detas*detas

        matched = dR2s < threshold*threshold

        matched_idxs = idxs[matched]
        reco_idxs = ak.firsts(matched_idxs[:,:,'1'])
        
        has0 = ak.fill_none(reco_idxs == 0, False)
        has1 = ak.fill_none(reco_idxs == 1, False)

        unmatched = np.ones( (len(recoetas), 2), dtype=bool)
        unmatched[has0, 0] = False
        unmatched[has1, 1] = False

        return unmatched

    def readfakes(self):
        reco  = ['Pho_Gen_Eta', 'Pho_Gen_Phi', 
                 'pt','eta', 'phi',
                 'energy', "Pho_SCRawE",
                 'Pho_R9']
        hits = ['Hit_X_Pho1', 'Hit_X_Pho2',
                    'Hit_Y_Pho1', 'Hit_Y_Pho2',
                    'Hit_Z_Pho1', 'Hit_Z_Pho2',
                    'Hit_ES_X_Pho1', 'Hit_ES_X_Pho2',
                    'Hit_ES_Y_Pho1', 'Hit_ES_Y_Pho2',
                    'Hit_ES_Z_Pho1', 'Hit_ES_Z_Pho2',
                    'ES_RecHitEnPho1', 'ES_RecHitEnPho2',
                    'RecHitEnPho1', 'RecHitEnPho2',
                    'RecHitFracPho1', 'RecHitFracPho2']

        varnames = reco+hits

        hits_trimmed = ['Hit_X_', 'Hit_Y_', 'Hit_Z_', 
                       'Hit_ES_X_', 'Hit_ES_Y_', 'Hit_ES_Z_',
                       'RecHitEn', 'RecHitFrac',
                       'ES_RecHitEn']

        arrs = self.tree.arrays(varnames, 'nPhotons==2')

        unmatched = self.gen_unmatch(arrs['Pho_Gen_Phi'], arrs['phi'], arrs['Pho_Gen_Eta'], arrs['eta'])
        Pho0 = unmatched[:,0]
        Pho1 = unmatched[:,1]

        result = {}

        reco = ['pt','eta','phi','energy', 'Pho_SCRawE']
        for var in reco:
            arrs[var] = ak.to_regular(arrs[var])
            result[var] = ak.to_numpy(ak.concatenate( (arrs[var][Pho0,0], arrs[var][Pho1,1])))
            print(var,result[var].shape)

        #return result

        for var in hits_trimmed:
            Pho0Name = var+'Pho1'
            Pho1Name = var+'Pho2'
            bettername = var
            if var[-1]=='_':
                bettername = var[:-1]
            
            result[bettername] = ak.concatenate( (arrs[Pho0Name][Pho0], arrs[Pho1Name][Pho1]) )

        result['subdet'] = np.abs(ak.to_numpy(ak.firsts(result['Hit_Z']))) < 300

        print("Dumping...");
        for var in result.keys():
            t0 = time()
            varname = var
            with open("%s/%s.pickle"%(self.outfolder, varname), 'wb') as f:
                pickle.dump(result[var], f, protocol = 4)
            print("\tDumping %s took %f seconds"%(varname, time()-t0))

        return result

    def readHGCAL(self):
        print("Reading in HGCAL branches:")

        t0=time()
        print("Reading Hcal rechit_x...")
        HcalRechit_X = self.tree['HcalRechit_posx'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

                
        t0=time()
        print("Reading Hcal rechit_y...")
        HcalRechit_Y = self.tree['HcalRechit_posy'].array()
        print("\ttook %0.3f seconds"%(time()-t0))
        
        t0=time()
        print("Reading Hcal rechit_z...")
        HcalRechit_Z = self.tree['HcalRechit_posz'].array()#comb_rechit_z_trimAhcal'].array()#combined_rechit_z'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        
        t0=time()
        print("Reading Ecal rechit_x...")
        EcalRechit_X = self.tree['EcalRechit_posx'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Ecal rechit_y...")
        EcalRechit_Y = self.tree['EcalRechit_posy'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Ecal rechit_z...")
        EcalRechit_Z = self.tree['EcalRechit_posz'].array()#comb_rechit_z_trimAhcal'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Hcal rechit_energy...")
        HcalrecHitEn = self.tree['HcalRechit_E'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Hcal cluster eta ...")
        HcalPFclustereta = self.tree['HcalPFclustereta'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Ecal cluster eta ...")
        EcalPFclustereta = self.tree['EcalPFclustereta'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Ecal rechit_energy...")
        EcalrecHitEn = self.tree['EcalRechit_E'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading trueBeanEnergy...")
        trueE = self.tree['true'].array()#trueBeamEnergy'].array()
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Ecal Energy...")
        #ecal = self.tree['ecalEn'].array()#trueBeamEnergy'].array()                                                
        ecal = self.tree['ecal'].array()#trueBeamEnergy'].array()                                                
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading Hcal Energy...")
        #hcal = self.tree['hcalEn'].array()#trueBeamEnergy'].array()                                                
        hcal = self.tree['hcal'].array()#trueBeamEnergy'].array()                                                
        print("\ttook %0.3f seconds"%(time()-t0))


        t0=time()
        print("Reading trueBeam eta...")
        eta = self.tree['eta'].array()#trueBeamEnergy'].array()                                                 
        print("\ttook %0.3f seconds"%(time()-t0))

        t0=time()
        print("Reading BeanEnergy...")
        beamEn = self.tree['HcalRechit_totE'].array()#beamEnergy'].array()
        print("\ttook %0.3f seconds"%(time()-t0))
       ############################################## 
       ## define boolean indexing array to select particular elements out of true energy
        #b = np.logical_and(trueE_orig>=20,trueE_orig<=300)
        #b =trueE_orig>=20
        #b =trueE_orig<=300
        t0 = time() 
        # Hit_X = Hit_X_orig
        # Hit_Y = Hit_Y_orig
        # Hit_Z = Hit_Z_orig
        # recHitEn=recHitEn_orig
        # trueE=trueE_orig
        # print(trueE, "type")
        # print(trueE_orig,"type-2")
        # SsLocation=SsLocation_orig
        # beamEn=beamEn_orig
        print()
        print("Building feature matrices...")
        #frac =((Hit_Z<54)*0.0105) + (np.logical_and(Hit_Z>54, Hit_Z<154)*0.0812) + ((Hit_Z>154)*0.12574)                             
        #recHitEn = recHitEn*frac
       
        recHitEn=ak.concatenate((HcalrecHitEn, EcalrecHitEn), axis=1)
        Hit_X=ak.concatenate((HcalRechit_X, EcalRechit_X), axis=1)
        Hit_Y=ak.concatenate((HcalRechit_Y, EcalRechit_Y), axis=1)
        Hit_Z=ak.concatenate((HcalRechit_Z, EcalRechit_Z), axis=1)
        clusterEta=ak.concatenate((HcalPFclustereta, EcalPFclustereta), axis=1)

        t0 = time()
        print("Max of Hit X : ",ak.max(Hit_X)," | Min of Hit X : ",ak.min(Hit_X))
        print("Max of Hit Y : ",ak.max(Hit_Y)," | Min of Hit Y : ",ak.min(Hit_Y))
        print("Max of Hit Z : ",ak.max(Hit_Z)," | Min of Hit Z : ",ak.min(Hit_Z))
        print("Max of recHitEn : ",ak.max(recHitEn)," | Min of recHitEn : ",ak.min(recHitEn))
        print("Max of Eta : ",ak.max(clusterEta)," | Min of Eta : ",ak.min(clusterEta))

        print("True E : ",trueE)
        trueEH = trueE[ecal>0]
        trueH = trueE[ecal==0]
        print("trueEH len : ",len(trueEH))
        print("trueH len : ",len(trueH))
        print("Hit X len : ",len(Hit_X))
        print("trueE len : ",len(trueE))
        print("Raw E len : ",len(hcal))
        print("Rechit E len : ",len(recHitEn))

        cf = cartfeat_HGCAL(Hit_X, Hit_Y, Hit_Z, recHitEn)
        #cf = cartfeat_HGCAL( Hit_Z, recHitEn)
        print("\tbuilding matrices took %0.3f seconds"%(time()-t0))
        
        t0 = time()
        cf = torchify(cf)  
        print("\tcasting to torch objects took %0.3f seconds"%(time()-t0))

        print()

        print("Building targets...")
        t0 = time()
        rawE = ak.sum(recHitEn, axis=1)
       
        #ratio = trueE/rawE
        #trueE = trueE[trueE<150]
#        print(.where(rawE<=0),1, rawE)
#        rawE = rawE[eta<1]
#        trueE= trueE[eta<1]
#        print("========= After condition ============")
        print(" recHitEn : ",recHitEn)
        print(" clusterEta : ",clusterEta)
        print("Length of Hit X : ",len(Hit_X))
        print("Length of Hit Y : ",len(Hit_Y))
        print("Length of Hit Z : ",len(Hit_Z))
        print("Length of Eta : ",len(clusterEta))

        print("Raw E : ",rawE)
        rawEH = rawE[ecal>0]
        rawH =rawE[ecal==0]
        print("RawEH len : ",len(rawEH))
        print("RawH len : ",len(rawH))

        ratio = trueE/rawE
        ratioflip = rawE/trueE        
        print(ratioflip)
        logratioflip = np.log(ratioflip)
        print("\tTook %0.3f seconds"%(time()-t0))


        print()

        print("Dumping:")
        t0=time()
        with open("%s/HcalRechit_X.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(HcalRechit_X, f, protocol = 4)
        print("\tDumped HcalRechit_X in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/HcalRechit_Y.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(HcalRechit_Y, f, protocol = 4)
        print("\tDumped HcalRechit_Y in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/HcalRechit_Z.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(HcalRechit_Z, f, protocol = 4)
        print("\tDumped HcalRechit_Z in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/EcalRechit_X.pickle"%self.outfolder, 'wb') as f:
                pickle.dump(EcalRechit_X, f, protocol = 4)
        print("\tDumped EcalRechit_X in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/EcalRechit_Y.pickle"%self.outfolder, 'wb') as f:
                pickle.dump(EcalRechit_Y, f, protocol = 4)
                print("\tDumped EcalRechit_Y in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/EcalRechit_Z.pickle"%self.outfolder, 'wb') as f:
                pickle.dump(EcalRechit_Z, f, protocol = 4)
        print("\tDumped EcalRechit_Z in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/HcalrecHitEn.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(HcalrecHitEn, f, protocol = 4)
        print("\tDumped HcalrecHitEn in %0.3f seconds"%(time()-t0))
        
        t0=time()
        with open("%s/EcalrecHitEn.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(EcalrecHitEn, f, protocol = 4)
        print("\tDumped EcalrecHitEn in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/clusterEta.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(clusterEta, f, protocol = 4)
        print("\tDumped clusterEta in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/recHitEn.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(recHitEn, f, protocol = 4)
        print("\tDumped recHitEn in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/Hit_X.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(Hit_X, f, protocol = 4)
        print("\tDumped Hit_X in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/Hit_Y.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(Hit_Y, f, protocol = 4)
        print("\tDumped Hit_Y in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/Hit_Z.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(Hit_Z, f, protocol = 4)
        print("\tDumped Hit_Z in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/trueE.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(trueE, f, protocol = 4)
        print("\tDumped trueE in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/ecal.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(ecal, f, protocol = 4)
        print("\tDumped ecal in %0.3f seconds"%(time()-t0))

        with open("%s/hcal.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(hcal, f, protocol = 4)
        print("\tDumped hcal in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/rawEH.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(rawEH, f, protocol = 4)
        print("\tDumped rawEH in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/rawH.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(rawH, f, protocol = 4)
        print("\tDumped rawH in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/trueEH.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(trueEH, f, protocol = 4)
        print("\tDumped trueEH in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/trueH.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(trueH, f, protocol = 4)
        print("\tDumped trueH in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/eta.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(eta, f, protocol = 4)
        print("\tDumped eta in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/beamEn.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(beamEn, f, protocol = 4)
        print("\tDumped beamEn in %0.3f seconds"%(time()-t0))
        t0=time()
        with open("%s/rawE.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(rawE, f, protocol = 4)
        print("\tDumped rawE in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/trueE_target.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(trueE, f, protocol = 4)
        print("\tDumped trueE target in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/ratio_target.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(ratio, f, protocol = 4)
        print("\tDumped ratio target in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/ratioflip_target.pickle"%self.outfolder, 'wb') as f:
            pickle.dump(ratioflip, f, protocol = 4)
        print("\tDumped ratioflip target in %0.3f seconds"%(time()-t0))

#        t0=time()
#        with open("%s/logratioflip_target.pickle"%self.outfolder, 'wb') as f:
#            pickle.dump(logratioflip, f, protocol = 4)
#        print("\tDumped logratioflip target in %0.3f seconds"%(time()-t0))

        t0=time()
        with open("%s/cartfeat.pickle"%self.outfolder, 'wb') as f:
            torch.save(cf, f, pickle_protocol = 4)
        print("\tDumped features in %0.3f seconds"%(time()-t0))

        print()
        return    

    def read(self, kind, N=None):
        varnames = varlists[kind]
        readcut = readcuts[kind]

        t0 = time()
        print("Reading in %s..."%kind)
        arrs = self.tree.arrays(varnames, readcut, entry_stop = N)


        gen = []
        reco = []
        event = []
        hits = []

        result = {}

        for var in arrs.fields:
            if var[-4:-1] == 'Ele' or var[-4:-1] == 'Pho': #hit level information 
                                    #different from reco: branches named "ele1" "ele2"
                name = var[:-4]
                hits.append(name) 
                continue
            elif var[:7] == 'Ele_Gen' or var[:7] == 'Pho_Gen': #gen level information
                gen.append(var)
            elif var == 'rho' or var == 'nElectrons' or var == 'nPhotons': #event level information
                event.append(var)
            else: #reco level information
                reco.append(var) 

        print("\tio took %f seconds"%(time()-t0))

        if hasgen[kind]:
            t0=time()
            if isEle[kind]:
                matched_idxs = self.gen_match(arrs['Ele_Gen_Phi'], arrs['phi'],
                                              arrs['Ele_Gen_Eta'], arrs['eta'])
            else:
                matched_idxs = self.gen_match(arrs['Pho_Gen_Phi'], arrs['phi'],
                                              arrs['Pho_Gen_Eta'], arrs['eta'])
            gen_idxs = matched_idxs[:,:,'0']
            reco_idxs = matched_idxs[:,:,'1']
            print("\tgen matching took %f seconds"%(time()-t0))

            t0 = time()
            
            for var in gen:
                arrs[var] = arrs[var][gen_idxs]

            for var in reco:
                print(var, arrs[var].type)
                arrs[var] = arrs[var][reco_idxs]

            print("\tapplying gen matching took %f seconds"%(time()-t0))

        t0 = time()

        #it can happen that there is exactly 1 reco electron, but it is identified as Ele2
        #I have no idea why, and it's super annoying, but here we are
        if not isEle[kind]:
            noEle1 = ak.num(arrs['iEtaPho1']) == 0
        else:
            noEle1 = ak.num(arrs['iEtaEle1']) == 0

        if hasgen[kind]:
            Ele1 = np.logical_and(reco_idxs == 0, ~noEle1) #recoidx==0 and there is an ele1
        else:
            Ele1 = ak.local_index(arrs['phi']) == 0

        Ele2 = ~Ele1  

        eventEle1 = ak.any(Ele1, axis=1)
        eventEle2 = ak.any(Ele2, axis=1)

        for var in gen + reco: #per-particle information, flattened 
            result[var] = ak.to_numpy(
                    ak.concatenate( 
                        (ak.flatten(arrs[var][Ele1]), 
                         ak.flatten(arrs[var][Ele2]) )))

        for var in event: #per-event information, broadcasted and flattened
            result[var] = ak.to_numpy(
                    ak.concatenate( 
                        (arrs[var][eventEle1], 
                         arrs[var][eventEle2]) ))

        for var in hits: #hit information, flattened
                         #note that this stays in awkward array format, while everything else is np
            if isEle[kind]:
                nameEle1 = var + 'Ele1'
                nameEle2 = var + 'Ele2'
            else:
                nameEle1 = var + 'Pho1'
                nameEle2 = var + 'Pho2'
            if var[-1] == '_':
                name = var[:-1]
            else:
                name = var
            result[name] = ak.concatenate( (arrs[nameEle1][eventEle1], arrs[nameEle2][eventEle2]) )

        print("\tbroadcasting and flattening took %f seconds"%(time()-t0))


        t0 = time()

        eventEle1 = ak.to_numpy(eventEle1)
        eventEle2 = ak.to_numpy(eventEle2) 

        #event idx
        #usefuly mostly for troubleshooting
        result['eventidx']= np.concatenate( (eventEle1.nonzero()[0], eventEle2.nonzero()[0]) )

        #hit subdetector
        #1: barrel 0: endcaps
        result['subdet'] = np.abs(ak.to_numpy(ak.firsts(result['Hit_Z']))) < 300

        print("\tdetermening aux features took %f seconds"%(time()-t0))

        t0 = time()

        savecut = savecuts[kind](result)

        for var in result.keys():
            result[var] = result[var][savecut]

        print("\tapplying savecut took %f seconds"%(time()-t0))
        
        print("Dumping...");
        for var in result.keys():
            t0 = time()
            varname = var
            with open("%s/%s.pickle"%(self.outfolder, varname), 'wb') as f:
                pickle.dump(result[var], f, protocol = 4)
            print("\tDumping %s took %f seconds"%(varname, time()-t0))

        t0 = time()

        print("Building cartesian features..")
        cf = cartfeat(result['Hit_X'], result['Hit_Y'], result['Hit_Z'], result['RecHitEn'], result['RecHitFrac'])
        print("\tBuilding features took %f seconds"%(time()-t0))
        t0 = time()
        result['cartfeat'] = torchify(cf)
        print("\tTorchifying took %f seconds"%(time()-t0))
        t0 = time()
        with open("%s/cartfeat.pickle"%(self.outfolder), 'wb') as f:
            torch.save(result['cartfeat'], f, pickle_protocol = 4)
        print("\tDumping took %f seconds"%(time()-t0))

        print("Building projective features..")
        pf = projfeat(result['Hit_Eta'], result['Hit_Phi'], result['Hit_Z'], result['RecHitEn'], result['RecHitFrac'])
        print("\tBuilding features took %f seconds"%(time()-t0))
        t0 = time()
        result['projfeat'] = torchify(pf)
        print("\tTorchifying took %f seconds"%(time()-t0))
        t0 = time()
        with open("%s/projfeat.pickle"%(self.outfolder), 'wb') as f:
            torch.save(result['projfeat'], f, pickle_protocol = 4)
        print("\tDumping took %f seconds"%(time()-t0))

        print("Building local features..")
        lf = localfeat(result['iEta'], result['iPhi'], result['Hit_Z'], result['RecHitEn'], result['RecHitFrac'])
        print("\tBuilding features took %f seconds"%(time()-t0))
        t0 = time()
        result['localfeat'] = torchify(lf)
        print("\tTorchifying took %f seconds"%(time()-t0))
        t0 = time()
        with open("%s/localfeat.pickle"%(self.outfolder), 'wb') as f:
            torch.save(result['localfeat'], f, pickle_protocol = 4)
        print("\tDumping took %f seconds"%(time()-t0))

        print()
        return result 
