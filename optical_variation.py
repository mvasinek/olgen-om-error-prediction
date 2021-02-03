
#otherwise stretch factor
class OpticalVariation:
    def __init__(self, mols, sf):
        self._mols = mols
        self._sf = sf

    def processOne(self):
        d = []
        for mol in self._mols:
            pos = mol.pos

            new_pos = []
            for p in pos:
                stretched_p = int(p*self._sf)                
                new_pos.append(stretched_p)                

            mol.pos = new_pos
            mol.length = int(mol.length*self._sf)

            d += self.convertToDist(mol.pos)

        keys, f0 = self.distF2(d)
        return self.distP2(keys, f0, len(d))

    def convertToDist(self, pos):
        ds = []
        for i in range(len(pos)-1):
            d = int(pos[i+1] - pos[i])            
            ds.append(d)

        return ds

    def distF2(self, ds):
        fs = {}

        for d in ds:
            if d in fs:
                fs[d] += 1
            else:
                fs[d] = 1

        keys = []
        for key in fs:
            keys.append(key)

        keys.sort()

        return (keys, fs)
    
    def distP2(self, keys, f0, dlen):
        p0 = [0]*25001
        for key in keys:
            if key <= 25000:
                p0[key] = f0[key]/dlen
        return p0