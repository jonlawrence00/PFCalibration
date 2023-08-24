#!/bin/python3 -u
import pickle

class tranche:
    def __init__(self, folder, subtranches = None):
        self.folder = folder

        if subtranches is None:
            self.read()

        else:
            self.combine(subtranches)
            self.write()

    def combine(self, subtranches):
        print("combining...")
        self.E = []
        self.etareco = []
        #self.featcart = []
        #self.featproj = []
        #self.featlocal = []
        self.phireco = []
        self.R9 = []
        self.rawE = []
        self.subdet = []

        for tr in subtranches:
            self.E += tr.E
            self.etareco += tr.etareco
            #self.featcart += tr.featcart 
            #self.featlocal += tr.featlocal
            #self.featproj += tr.featproj
            self.phireco += tr.phireco
            self.R9 += tr.R9
            self.rawE += tr.rawE
            self.subdet += tr.subdet

    def read(self):
        print("reading from %s..."%self.folder)
        with open("%s/energy_ecal_mustache_all.pickle"%self.folder, 'rb') as f:
            self.E = pickle.load(f)

        with open("%s/etareco_all.pickle"%self.folder, 'rb') as f:
            self.etareco = pickle.load(f)

        #with open("%s/features_noES_cart_multfrac_all.pickle"%self.folder, 'rb') as f:
        #    self.featcart = pickle.load(f)

        #with open("%s/features_noES_proj_multfrac_all.pickle"%self.folder, 'rb') as f:
        #    self.featproj = pickle.load(f)

        #with open("%s/features_noES_local_multfrac_all.pickle"%self.folder, 'rb') as f:
        #    self.featlocal = pickle.load(f)

        with open("%s/phireco_all.pickle"%self.folder, 'rb') as f:
            self.phireco = pickle.load(f)

        with open("%s/R9_all.pickle"%self.folder, 'rb') as f:
            self.R9 = pickle.load(f)

        with open("%s/rawE_all.pickle"%self.folder, 'rb') as f:
            self.rawE = pickle.load(f)

        with open("%s/subdet_all.pickle"%self.folder, 'rb') as f:
             self.subdet = pickle.load(f)

    def write(self):
        print("writing to %s..."%self.folder)
        with open("%s/energy_ecal_mustache_all.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.E, f)

        with open("%s/etareco_all.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.etareco, f)

        #with open("%s/features_noES_cart_multfrac_all.pickle"%self.folder, 'wb') as f:
        #    pickle.dump(self.featcart, f)

        #with open("%s/features_noES_proj_multfrac_all.pickle"%self.folder, 'wb') as f:
        #    pickle.dump(self.featproj, f)

        #with open("%s/features_noES_local_multfrac_all.pickle"%self.folder, 'wb') as f:
        #    pickle.dump(self.featlocal, f)

        with open("%s/phireco_all.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.phireco, f)

        with open("%s/R9_all.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.R9, f)

        with open("%s/rawE_all.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.rawE, f)

        with open("%s/subdet_all.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.subdet, f)

T0 = tranche('2018_UL_ZEE_Data/T0')
T1 = tranche('2018_UL_ZEE_Data/T1')
T2 = tranche('2018_UL_ZEE_Data/T2')
T3 = tranche('2018_UL_ZEE_Data/T3')
T4 = tranche('2018_UL_ZEE_Data/T4')
T5 = tranche('2018_UL_ZEE_Data/T5')

big = tranche('2018_UL_ZEE_Data', [T0,T1,T2,T3,T4,T5])
