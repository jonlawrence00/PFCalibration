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
        self.pred = []

        for tr in subtranches:
            self.pred += tr.pred

    def read(self):
        print("reading from %s..."%self.folder)

        with open("%s/pred.pickle"%self.folder, 'rb') as f:
             self.pred = pickle.load(f)

    def write(self):
        print("writing to %s..."%self.folder)

        with open("%s/pred.pickle"%self.folder, 'wb') as f:
            pickle.dump(self.pred, f)

T0 = tranche('ZeeDataT0')
T1 = tranche('ZeeDataT1')
T2 = tranche('ZeeDataT2')
T3 = tranche('ZeeDataT3')
T4 = tranche('ZeeDataT4')
T5 = tranche('ZeeDataT5')

big = tranche('ZeeData', [T0,T1,T2,T3,T4,T5])
