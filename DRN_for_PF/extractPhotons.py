import Extract

ex = Extract.Extract("/home/rusack/shared/pickles/2018_Gamma_Jet/realPhotons/", "/home/rusack/evans908/shared/nTuples/2018_Gamma_Jet/nTupleMC_Pt20toInf.root")
ex.read("gun_pho")

ex = Extract.Extract("/home/rusack/shared/pickles/2018_Gamma_Jet/fakePhotons/", "/home/rusack/evans908/shared/nTuples/2018_Gamma_Jet/nTupleMC_Pt20toInf.root")

ex.readfakes()
