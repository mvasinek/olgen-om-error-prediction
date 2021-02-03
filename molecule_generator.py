import copy
import random
import sys
import time

from bnx import BNX
from genome import OMGenome

from missing_site import EnzymeMissing

class Molecule:
    def __init__(self, pos, length, chr_origin=None, chr_origin_pos=None):
        self.pos = pos
        self.length = length
        self.chr_origin = chr_origin
        self.chr_origin_pos = chr_origin_pos

class MoleculeGenerator:
    def __init__(self, bnx_path, genome, coverage):
        self._mol_lens = BNX.CollectMolLengths(bnx_path)
        self._p_mols = self.__lengthsDistribution()
        #self.__storeLengthsDist()
        #BNX.CollectMolLengthsSites(bnx_path)
        self._genome = genome

        self.coverage = coverage

    def __lengthsDistribution(self):
        lim = 0
        for l in self._mol_lens:
            if l>lim:
                lim = int(l)

        #print("maximum molecule length: ", lim)
        f = [0]*(lim+1)
        for l in self._mol_lens:
            f[int(l)] += 1

        p = [0]*lim
        for i in range(lim):
            p[i] += f[i] / len(self._mol_lens)

        return p        

    def __storeLengthsDist(self):
        f = open("molecule_length.csv","w")

        for i in range(len(self._p_mols)):
            if i % 375 == 0:
                f.write("%d\t%.8f\n" % (i, self._p_mols[i]))
        f.close()

    def generate(self, erange, vmiss, sex):
        #pass through genome and create molecules 
        """
        1 - obtain distribution of mol lens - probabilities
        2 - generate length and take positions
        3 - scale positions from zero-based coordinate system
        """
        p_mols = self._p_mols
        c_lens = self._genome.chrlen

        molecules = []

        arr = []
        for i in range(len(p_mols)):
            arr.append(i)            

        #nachystame si kopii genome pos
        orig_pos = copy.deepcopy(self._genome.pos)

        for cov_id in range(self.coverage):            
            c_id = 0                        
            #pridame missing
            new_pos = None
            if vmiss != 0:
                mi = EnzymeMissing(None, erange, vmiss)
                mi.pos = orig_pos
                mi.chrlen = self._genome.chrlen
                new_pos = mi.processOne(None)
            else:
                new_pos = copy.deepcopy(orig_pos)

                    
            #zpracujeme do molekul            
            mol_lens = random.choices(arr, p_mols, k=27000)
            mol_len_id = 0
            for c_pos in new_pos:
                #female doesnt have chr24
                if sex == "female" and c_id == 23:
                    break

                if sex == "male":                    
                    if c_id == 22 and random.random() > 0.5:
                        c_id += 1
                        continue
                    if c_id == 23 and random.random() > 0.5:
                        c_id += 1
                        break

                cur = 0
                pos_id = 0
                c_len = c_lens[c_id]

                while cur < c_len:
                    #faze generovani molekul -> nejdrive delka                
                    mol_len = mol_lens[mol_len_id%27000] + random.randint(0, 375)
                    #pozice uvnitr chromozomu
                    start_pos = cur                

                    end_pos = start_pos + mol_len

                    #projit pozice v genomu
                    start_chr_pos = pos_id
                    while pos_id < len(c_pos) and c_pos[pos_id] < end_pos:
                        pos_id += 1

                    mol_positions = []
                    for j in range(start_chr_pos, pos_id):
                        mol_pos = c_pos[j] - start_pos
                        mol_positions.append(mol_pos)

                    #pridej molekulu                
                    molecules.append(Molecule(mol_positions, mol_len, c_id+1, start_pos))
                    
                    #novy start molekuly
                    cur += mol_len
                    mol_len_id += 1
                
                c_id += 1                        
        
        return molecules

    """
    distribution of the number of sites in mols
    """
    def sitesDist(self, mols, include_one=False):
        max_s = 0
        sites = []
        for mol in mols:
            if len(mol.pos) > max_s:
                max_s = len(mol.pos)

            sites.append(len(mol.pos))

        if max_s < 100:
            max_s = 100

        f = [0]*(max_s+1)

        for mol in mols:
            f[len(mol.pos)] += 1       
        
        if not include_one:
            sites_number = len(sites) - f[0] - f[1]
        else:
            sites_number = len(sites) - f[0]

        f[0] = 0
        if not include_one:
            f[1] = 0

        p = [0]*(max_s+1)

        if sites_number == 0:
            print(sites_number, len(mols), max_s)

        if sites_number != 0:
            for i in range(max_s+1):
                p[i] = f[i] / sites_number

        return p

    @staticmethod
    def SitesDistBNX(sites_a):
        max_s = 0        
        for sites_no in sites_a:
            if sites_no > max_s:
                max_s = sites_no

        f = [0]*(max_s+1)

        for site_no in sites_a:
            f[site_no] += 1

        p = [0]*(max_s+1)

        for i in range(max_s+1):
            p[i] = f[i] / len(sites_a)

        return p

    @staticmethod
    def StoreDist(mg, mols, bnx_sites_dist):
        f = open("artificial.csv","w")
        sites_dist = mg.sitesDist(mols)

        for i in range(len(sites_dist)):
            f.write("%d\t%.8f\t%.8f\n" % (i,sites_dist[i], bnx_sites_dist[i]))

        f.close()

    @staticmethod
    def StoreMolLengthDist(mols):
        lim = 0
        for mol in mols:
            if mol.length>lim:
                lim = int(mol.length)
        
        f = [0]*(lim+1)
        for mol in mols:
            f[int(mol.length)] += 1

        p = [0]*lim
        for i in range(lim):
            p[i] += f[i] / len(mols)

        fo = open("simul-mol-lens.csv","w")
        for i in range(lim):
            if i % 375 != 0:
                continue
            fo.write("%d\t%.6f\n" % (i,p[i]))
        fo.close()