"""
Vyuzit tridy OMGenome a insertion_site a missing site pro vygenerovani 
umelych dat.
"""

from genome import OMGenome
from insertion_site import InsertionSite
from missing_site import MissingSite
from sizing_error import SizingError

import random
import sys

from multiprocessing import Queue,Process, freeze_support

class Molecule:
    def __init__(self, mol_beg, mol_length):
        self._sites = []
        #molecule beginning with respect to reference genome
        self._mol_beg = mol_beg
        self._mol_length = mol_length

    def addSite(self, site_pos):
        #convert relatively, site pos with respect to reference genome
        self._sites.append(site_pos - self._mol_beg)

    def serialize(self):
        #molecule_id - mol_length - run_id - number_of_sites - *sites
        s = "0 " + str(self._mol_length) + " 0 " + str(len(self._sites))

        for site in self._sites:
            s += " " + str(site)

        return s

    @staticmethod
    def ToString(molecule):
        return molecule.serialize()

class Simulator:
    def __init__(self, gpath, min_mol_len=15000, max_mol_len=1000000):
        self._gpath = gpath        
        self._min_len = min_mol_len
        self._max_len = max_mol_len

    def generate(self, outpath, irate, mrate, coverage):
        #nacist genom
        omg = OMGenome(self._gpath)

        mols = []

        for i in range(coverage):
            print("%d out of %d" % (i,coverage))
            queue = Queue()
            #spustit missing
            
            mi = MissingSite(None, mrate)
            mi.pos = omg.pos
            mi.chrlen = omg.chrlen
            mi.setPD(25000*[mrate])
            #mi.setPD(mi.bezierPD(0.7757, 0.2415, 0.3287, 0.1576, 1225))            
            
            p = Process(target = mi.processOne, args=(queue,))
            p.start()        
            
            pos_mi = queue.get()        
            p.join()
            
            #spustit inserce
            ins = InsertionSite(None, omg.chrlen, pos_mi)               
            ins.pos = pos_mi
            ins.expectedInsertions(irate)

            p = Process(target = ins.processOne, args=(queue,))
            p.start()
            pos_ins = queue.get()
            p.join()
            #nahodne rozkouskovat do molekul
            #zapsat do souboru

            siz = SizingError(None, 100, omg.chrlen, pos_ins)

            siz.pos = pos_ins        
            p = Process(target = siz.processOne, args=(queue,))
            p.start()        
            pos_siz = queue.get()        
            p.join()
            
            mols += self._intoMolecules(pos_siz, omg.chrlen, coverage)
        
        s_mols = map(Molecule.ToString, mols)

        f = open(outpath, "w")
        f.write("\n".join(s_mols))
        f.close()

    def _intoMolecules(self, pos, chr_len, coverage):
        mols = []
        
        #process pos - chromosome by chromosome            
        chr_id = 0
        for c_pos in pos:
            #find the first position on chromosome
            d = c_pos[0]
            if d > self._max_len:
                d = self._max_len
                
            #pujde napric celym chromozomem a na zaver vyfiltruju molekuly aby mely alespon dve znacky
            #pruchod chromozomem
            current_pos = 0
            mol_len = random.randint(self._min_len, self._max_len)

            pos_id = 0

            while(current_pos + mol_len < chr_len[chr_id]):
                #vygeneruj molekulu projdi pozice a zkontroluj ktere jsou uvnitr molekuly
                mol = Molecule(current_pos, mol_len)
                mol_end = current_pos + mol_len
                while pos_id < len(c_pos) and c_pos[pos_id] < mol_end:
                    mol.addSite(c_pos[pos_id])
                    pos_id += 1

                mols.append(mol)
                #zaktualizuj pozice
                current_pos += mol_len
                #delka nove molekuly
                mol_len  = random.randint(self._min_len, self._max_len)

            chr_id += 1

        return mols

if __name__ == '__main__':
    freeze_support()
    s = Simulator(sys.argv[1])
    s.generate(sys.argv[2], 0.05,0.9,400)