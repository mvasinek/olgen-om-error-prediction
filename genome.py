import random
import sys
import time

from bnx import BNX, BNXData

"""
Class responsible for representation of reference genome fragment lengths distribution.
- class assumes human genome with 24 chromosomes.
"""
class OMGenome:
    def __init__(self, g_path):
        self.genome_path = g_path
        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        
        d = self.convertToDist(self.pos)        
        self.dlen = len(d)
        #self.max_d = self.maxDistance(self.d)
        #self.f0 = self.distF(self.d, self.max_d)
        keys, f0 = self.distF2(d)        
        #self.p0 = self.distP(self.f0, self.max_d, self.dlen)
        self.p0 = self.distP2(keys, f0, self.dlen)
        #self.portions = self.weightedLen2(keys, f0)        

    """
    Parse positions from reference cmap file.
    """
    def parsePositions(self):
        f = open(self.genome_path, "r")
        lines = f.readlines()
        f.close()

        pos = [[]]        
        chr_id = 0
        
        for line in lines:
            if line.startswith("#"):
                continue

            line_values = line.strip().split("\t")
            self.chrlen[int(line_values[0])-1] = int(float(line_values[1]))                

            #end of last chromosome marker.
            if int(line_values[3]) > int(line_values[2]):
                if line_values[0] == "24":
                    break
                pos.append([])
                chr_id += 1
                continue
            else:
                pos[chr_id].append(int(float(line_values[5])))

        return pos

    """
    Returns genome size
    """
    def genomeSize(self):
        gs = 0
        for csize in self.chrlen:
            gs += csize
        return gs

    """
    Returns number of labels in the genome.
    """
    def labelsCount(self):
        lc = 0
        for chrpos in self.pos:
            lc += len(chrpos)
        return lc

    def insertionRate(self, fpr):
        gsize = self.genomeSize()
        lc = self.labelsCount()

        return (lc*fpr)/(gsize-lc)

    """
    Converts positions of restriction sites into distribution of fragment lengths.
    """
    def convertToDist(self, positions):
        d = []
        for chr_pos in positions:
            if len(chr_pos) == 0:
                continue

            pos_l = chr_pos[0]

            for i in range(1,len(chr_pos)):
                d.append(chr_pos[i] - pos_l)                
                pos_l = chr_pos[i]

        return d

    """
    Finds the largest distance between restriction sites.
    """
    def maxDistance(self, ds):
        maximum = 0
        for d in ds:
            if d > maximum:
                maximum = d
        return maximum

    """
    Counts frequency of fragment lengths.
    """
    def distF(self, dist, maxd):
        f0 = [0]*(maxd+1)
        for d in dist:
            f0[d] += 1
        return f0

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


    """
    Counts probability distribution of fragment lengths.
    """
    def distP(self, f0, maxd, dlen):
        p0 = [0]*(maxd+1)
        for i in range(maxd+1):
            if dlen > 0:
                p0[i] = f0[i]/dlen
        return p0

    def distP2(self, keys, f0, dlen):
        p0 = [0]*25001
        for key in keys:
            if key <= 25000:
                p0[key] = f0[key]/dlen
        return p0

    """
    correction from distribution
    """
    #def correctionFromDistribution(self, p, emols, g_path, genome=None, g_pos=None):
    def correctionFromDistribution(self, p, emols, genome):
        #ratio of removed intervals with respect to total labels count
        l = genome.labelsCount()
        r = emols / l

        z = l / (l*(1-r))

        #portions of each fragment length on the total genome length
        #portions = self.weightedLen(genome)
        portions = genome.portions

        #out distribution
        p_out = [0]*25000

        for i in range(25000):
            try:
                p_out[i] = abs((1/z) * (p[i] - r*portions[i]))
            except:
                print("err correction: ", len(p_out), i, len(p), len(portions))

        """
        p_rem = 0
        for i in range(750,938):
            p_rem += p_out[i]
            p_out[i] = 0

        
        factor = 1/(1- p_rem)

        for i in range(938, 25000):
            p_out[i] *= factor
        """

        return p_out

    """
    returns distribution of distances given positions
    """
    def distFromPositions(self, pos, correction=False, b_path=None, g_path=None):
        d = self.convertToDist(pos)
        if correction:            
            d = self.correct(d, b_path, g_path)

        dlen = len(d)
        max_d = self.maxDistance(d)
        f0 = self.distF(d, max_d)
        return self.distP(f0, max_d, dlen)

    """
    Function returns portion of bp each interval covers in genome.
    """
    def weightedLen(self, genome):
        #unique non-zero intervals
        #f0 = genome.f0
        d = self.convertToDist(self.pos)
        max_d = self.maxDistance(d)
        f0 = self.distF(d, max_d)
        uniq_f = []

        for l in range(len(f0)):
            if f0[l] != 0:
                uniq_f.append(l)

        w_sum = 0

        portion = [0]*len(f0)
        for l in uniq_f:
            portion[l] = l*f0[l]
            w_sum += portion[l]
            #w_sum += l*f0[l]

        
        for l in uniq_f:
            #portion[l] = l*f0[l] / w_sum
            portion[l] /= w_sum

        return portion

    def weightedLen2(self, keys, f0):
        portion = [0]*25001
        w_sum = 0

        for l in keys:
            mul = l*f0[l]
            if l <= 25000:
                portion[l] = mul
            w_sum += mul           
        
        for l in keys:    
            if l <= 25000:        
                portion[l] /= w_sum

        return portion

    def pCorrections(self, emols, portions, dlen):
        eCorrections = [0]*len(portions)
        for l in range(len(portions)):
            eCorrections[l] = emols*portions[l]/dlen
        return eCorrections

    def correct(self, ds, b_path, g_path):        
        emols,genome = self._eMols(b_path, g_path)

        portions = self.weightedLen(genome)
        #p_corrs = self.pCorrections(emols, portions, len(ds))
        
        dout = []
        out = 0

        l = []
        p_reduced = []
        for i in range(len(portions)):
            if portions[i] > 0:
                l.append(i)
                p_reduced.append(portions[i])

        while out < emols:
            #generate random interval
            interval = random.choices(l, p_reduced,k=1)[0]
            try:
                ds.remove(interval)
                out += 1
            except:
                continue
        
        return ds

    def _eMols(self, bnxpath, g_path, g_pos=None):
        g = None
        L = None

        if g_pos == None:
            g = OMGenome(g_path)
            L = g.intervalLength()
        else:
            L = g.intervalLength(g_pos)        
        
        e_mols = None
        
        if bnxpath.endswith(".bnx"):
            t0 = time.time()
            mol_lens = BNX.CollectMolLengths(bnxpath)            
            t1 = time.time()
            print("Mol collect in: %.4f" % (t1-t0))
            mol_dist = BNX.MolsToDist(mol_lens)
            t2 = time.time()
            print("Mol dist in: %.4f" % (t2-t1))
            e_mols = BNX.ExpectedMolsPerGenome(L, mol_dist)            
            t3 = time.time()
            print("Expected mols in: %.4f" % (t3-t2))
        elif bnxpath.endswith(".data"):
            m_lens = BNXData.CollectMolLengths(bnxpath)            
            mols_dist = BNXData.MolsToDist(m_lens)            
            e_mols = BNXData.ExpectedMolsPerGenome(L, mols_dist)

        return (e_mols, g)

    """
    Length of intervals covered by neighbouring labels. Without telomers where no labels are present.
    """
    def intervalLength(self, pos=None):
        if pos == None:
            pos = self.pos
        d = 0

        for chr_pos in pos:
            d += chr_pos[-1] - chr_pos[0]

        return d