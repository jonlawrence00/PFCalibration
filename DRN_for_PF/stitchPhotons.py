#import the stuff
import pandas as pd #dataframes etc
import matplotlib.pyplot as plt #plotting
import pickle
import numpy as np
import os, sys
#list the pickles
fakeDir = "/home/rusack/shared/pickles/2018_Gamma_Jet/fakePhotons/"
realDir = "/home/rusack/shared/pickles/2018_Gamma_Jet/realPhotons/"
combDir = "/home/rusack/shared/pickles/2018_Gamma_Jet/combPhotons/"

fakePickles = os.listdir(fakeDir)
realPickles = os.listdir(realDir)
#every fake pickle has a real counterpart
for ipickle in fakePickles:
    print(ipickle+":")
    print("opening...",end='')
    fakePickle = open(fakeDir+ipickle, "rb")
    realPickle = open(realDir+ipickle, "rb")
    #extract arrays
    print("extracting...",end='')
    fakeArray = pickle.load(fakePickle)
    realArray = pickle.load(realPickle)
    #make combined array
    print("combining...",end='')
    combArray = np.concatenate((fakeArray,realArray))
    #pickle combined array
    print("dumping")
    pickle.dump(combArray, open(combDir+ipickle, "wb"))
    #now we count and make the label pickle
print("=================================")
print("counting...  ",end="")
print(len(fakeArray),"fakes", len(realArray),"reals")
fakeArray = np.zeros(len(fakeArray))
realArray = np.ones(len(realArray))

combArray = np.concatenate((fakeArray, realArray))
print("...dumping...")

pickle.dump(combArray, open(combDir+"labels_target.pickle", "wb"))
