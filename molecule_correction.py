import random

from genome import OMGenome
from bnx import BNX, BNXData

class MoleculeCorrection:
    def __init__(self, gpath, bnxpath):
        self._gpath = gpath
        self._bnxpath = bnxpath

    def correct(self, ds):
        emols = self._eMols()
        pnot = emols / len(ds)
        dout = []

        for d in ds:
            if random.random() > pnot:
                dout.append(d)
        
        return d

    def _eMols(self):
        g = OMGenome(self._gpath)

        L = g.intervalLength()
        e_mols = None
        
        if self._bnxpath.endswith(".bnx"):
            e_mols = BNX.ExpectedMolsPerGenome(L, BNX.MolsToDist(BNX.CollectMolLengths(self._bnxpath)))
        elif self._bnxpath.endswith(".bnxdata"):
            e_mols = BNXData.ExpectedMolsPerGenome(L, BNXData.MolsToDist(BNXData.CollectMolLengths(self._bnxpath)))

        return (e_mols, g.labelsCount())