import Extract

ex = Extract.Extract("/home/rusack/shared/pickles/realPhotons/", "/home/rusack/shared/nTuples/2018_Gamma_Jet/nTupleMC_Pt-20to40.root")
ex.read("gun_pho")

ex = Extract.Extract("/home/rusack/shared/pickles/fakePhotons/", "/home/rusack/shared/nTuples/2018_Gamma_Jet/nTupleMC_Pt-20to40.root")

ex.readfakes()
