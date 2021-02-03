import math
import os
import random
import sys
import numpy

from multiprocessing import Queue, Process, freeze_support

from genome import OMGenome

class SizingError(OMGenome):
    def __init__(self, gpath, sigma, chrlen = None, pos=None):
        if gpath != None:
            super().__init__(gpath)

        self.genome_path = gpath
        
        self.chrlen  = chrlen
        self.pos = pos
        self.lim = 25000

        self.sigma = int(sigma)

        self.gsize = self.genomeSize()
        self.lc = self.labelsCount()

    """
    Function returns fragment length modified by normal distribution.
    """
    def __dNormalError(self, d0):
        return int(random.normalvariate(d0, self.sigma))

    def _dError(self, d0):
        return self.__dNormalError(d0)

    def _linError(self, d0):
        return d0 + random.randint(-1*self.sigma, self.sigma)

    def processOne(self, queue):
        p1 = []
        for chr_id in range(len(self.pos)):
            p1.append([])

            current_chr = self.pos[chr_id]

            last = 0
            pos_id = 0
            while pos_id < len(current_chr):
                position = current_chr[pos_id]
                if pos_id > 0:
                    diff = position - current_chr[pos_id - 1]
                else:
                    diff = position

                #new_diff = self._dError(diff)
                new_diff = self._linError(diff)
                new_pos = last + new_diff
                if pos_id < len(current_chr) - 1:
                    if new_pos > current_chr[pos_id + 1] - 6:
                        new_pos = current_chr[pos_id + 1] - 6
                    
                if new_pos < last + 6:
                    new_pos = last + 6

                if new_pos > self.chrlen[chr_id] - 6:
                    #print(chr_id, new_pos, self.chrlen[chr_id])
                    new_pos = self.chrlen[chr_id] - 6
                
                p1[chr_id].append(new_pos)
                last = new_pos
                pos_id += 1

            
            """
            for i in range(len(p1[chr_id])-1):
                if p1[chr_id][i] >= p1[chr_id][i+1]:
                    print("sizing error: ")
                    print(i)
                    print(len(p1[chr_id]))
                    print(p1[chr_id][i])
                    print(p1[chr_id][i+1])                    
                    #print("Error: ", p1[chr_id][i], p1[chr_id][i+1], chr_id, len(p1[chr_id]), current_chr[pos_id-1], current_chr[pos_id])
            """


        queue.put(p1)

"""
Li et. al - Towards a More Accurate Error Model for BioNano Optical Maps
"""
class SizingErrorLi(SizingError):
    #based on the Li et. al paper
    LAPLACE_COEFS = [(0.858181, 0.180196), (0.980760, 0.071176), (1.003354, 0.052800), (1.00482, 0.042428)]
    def _dLaplaceError(self, r):
        """
        sk = o / r 
        o - observed fragment length
        r - reference fragment length
        o = r*sk
        """
        sk = None
        if r < 2400:
            sk = numpy.random.laplace(LAPLACE_COEFFS[0][0],LAPLACE_COEFFS[0][1])
        elif r >= 2400 and r < 3600:
            sk = numpy.random.laplace(LAPLACE_COEFFS[1][0],LAPLACE_COEFFS[1][1])
        elif r >= 3600 and r < 4800:
            sk = numpy.random.laplace(LAPLACE_COEFFS[2][0],LAPLACE_COEFFS[2][1])
        elif r >= 4800:
            sk = numpy.random.laplace(LAPLACE_COEFFS[3][0],LAPLACE_COEFFS[3][1])

        return int(sk*r)

    def _dError(self, d0):
        return self._dLaplaceError(d0)

class SizingErrorLaplace(SizingError):
    def __init__(self, gpath, levels, mb_tupple, chrlen=None, pos=None):
        super().__init__(gpath, 0, chrlen, pos)

        #levels is an list of upper bounds for each level
        self._levels = levels
        self._mb = mb_tupple

    def _dLaplaceError(self, r):
        sk = None

        l_id = 0        
        for l in self._levels:
            if r < l:                
                break
            
            l_id += 1

        m = self._mb[l_id][0]
        b = self._mb[l_id][1]

        sk = numpy.random.laplace(m, b)

        return int(sk*r)

    def _dError(self, d0):
        return self._dLaplaceError(d0)