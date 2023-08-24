class tranche:
    def __init__(self, folder):
        self.folder = folder
        self.data = {}
        self.feat = {}

        self.loadAll()

    def loadVariable(self, var):
        with open("%s/%s_all.pickle"%(self.folder, var)) as f:
            self.data[var] = pickle.load(f)

    def loadFeat(self, coords):
        with open("%s/features_noES_%s_multfrac_all.pickle"%(self.folder, coords)) as f:
            self.feat[coords] = pickle.load(f)

    def loadAll(self):
        for var in ['energy_ecal_mustache', 'etareco', 'phireco', 'R9', 'rawE', 'subdet']:
            self.loadVariable(var)

        for coords in ['cart', 'proj', 'local']:
            self.loadFeat(coords)
