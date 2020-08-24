import math
import os
import random
import sys

from multiprocessing import Queue, Process, freeze_support

from genome import OMGenome

class InsertionSite(OMGenome):
    def __init__(self, gpath, chrlen = None, pos=None):
        if gpath != None:
            super().__init__(gpath)

        self.chrlen  = chrlen
        self.pos = pos
        self.lim = 25000

        self.gsize = self.genomeSize()
        self.lc = self.labelsCount()
        self.ir = 0
        self.pfp = 0
        self.efp = 0

    """
    Spocitame ocekavanou cetnost fp znacek.
    - ffp = pocetZnacek*ir
    - pravdepodobnost ze se na konkretni pozici vyskytuje znacka pfp = ffp/gsize-pocetZnacek
    """
    def insertionRate(self, ir):
        self.ir = ir
        self.pfp = (self.lc*self.ir)/(self.gsize-self.lc)

    """
    Vraci ocekavanou cetnost fp
    """
    def expectedInsertions(self, ir):
       self.efp = ir*self.lc         

    def processOne(self, queue):
        p1 = []
        paux = []
        for chr_id in range(24):
            paux.append([])
            #ocekavany pocet fp na chromozom
            echr_fp = int(round(self.efp*self.chrlen[chr_id]/self.gsize))

            chr_fp = []
            for i in range(echr_fp):
                chr_fp.append(random.randint(1,self.chrlen[chr_id]))

            chr_fp.sort()

            #pomala varianta spoji pole a setridi je
            p1_c = chr_fp + self.pos[chr_id]
            p1_c.sort()

            #projdeme p1_c a odstranime duplicity
            i = 0
            lc = len(p1_c) - 1
            p1_a = []
            while i < lc:
                if p1_c[i] != p1_c[i+1]:
                    p1_a.append(p1_c[i])
                i += 1
                
            #posledni prvek pridavame vzdy
            p1_a.append(p1_c[lc])

            p1.append(p1_a)

        queue.put(p1)
