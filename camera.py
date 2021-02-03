import math
import os
import random
import sys
import time

from scipy.stats import binom
from scipy.stats import norm

from multiprocessing import Process, freeze_support

from genome import OMGenome


"""
Class responsible for conversion of fragment lengths into pixel image and back into operating resolution.
Methods mimics function of capturing system.
"""
class Camera(OMGenome):
    def __init__(self, g_path, pos=None, chrlen=None, resolution=450, ebpp=500):
        #expected base pairs per pixel
        self.bpp = resolution
        #operating resolution - names differs from the one in the article
        self.ebpp = ebpp
        self.genome_path = g_path

        if chrlen == None:
            self.chrlen = 24*[0]
        else:
            self.chrlen = chrlen

        self.weights = Camera.triangle(self.bpp)
        self.w_pop = [i for i in range(2*self.bpp)]

        if pos == None:
            self.d = self.discreteTransform(self.parsePositions())
        else:
            self.d = self.discreteTransform(pos)
        
        self.dlen = len(self.d)
        self.max_d = self.maxDistance(self.d)
        self.f0 = self.distF(self.d, self.max_d)
        self.p0ti = self.distP(self.f0, self.max_d, self.dlen)        

    """
    Creates pixel image from positions in a sequence.
    """
    def pixelImage(self, pos, off):
        c_id = 0
        pi = []        
        bo = self.bpp - off
        #bo = off        
        for chromosome in pos:            
            c_len = self.chrlen[c_id]

            pi_chr = (math.ceil((c_len + off) / self.bpp) + 1)*[0]

            for p in chromosome:
                x = (p+bo) // self.bpp
                pi_chr[x] = 1

            pi.append(pi_chr)

            c_id += 1            

        return pi

    """
    returns list of distances between bright pixels
    """
    def imageToDist(self, pi):
        d = []
        for c_id in range(len(pi)):
            d_chrom = []
            chromosome = pi[c_id]            
            
            #find first non/zero entry
            start_pos = 0            
            for i in range(len(chromosome)):
                if chromosome[i] != 0:
                    start_pos = i
                    break

            #compute all distances
            last_pos = start_pos
            for i in range(start_pos+1, len(chromosome)):
                if chromosome[i] != 0:                    
                    d_chrom.append(int(i - last_pos))
                    last_pos = i

            d.append(d_chrom)

        return d

    def distWithStretch(self, ds):
        #for each distance compute triangular guess
        d_stretch = []        

        for c_id in range(len(ds)):
            d_chrom = ds[c_id]
            d_chr_s = []

            for d in d_chrom:
                d_chr_s.append(self.triangularGuess(d, self.bpp))

            d_stretch.append(d_chr_s)

        return d_stretch

    @staticmethod
    def triangle(res):
        f = [0]*(2*res)
        t = 0
        for i in range(1,res+1):
            for j in range(res+1, res*2+1):                
                f[j-i] += 1
                t += 1

        for i in range(2*res):
            f[i] /= t

        return f

    """
    gives a random length provided by triangular distribution
    """    
    def triangularGuess(self, d, bpp):        
        return random.choices(self.w_pop, self.weights)[0] + (d-1)*bpp
    
    def autonoiseDistance(self, ds):    
        a_ds = []
        for c_ds in ds:
            #look for clusters with distance below 450
            a_d = []

            a_d_id = -1
            last_off = 0
            i = 0

            while i < len(c_ds):
                t = 0
                while i < len(c_ds) and c_ds[i] < 450:
                    t += c_ds[i]
                    i += 1

                if t == 0:
                    #no autonoise
                    a_d.append(c_ds[i] + last_off)
                    a_d_id += 1
                    i += 1
                    last_off = 0
                else:
                    last_off = int(t/2)
                    #autonoise distance exists already                    
                    if a_d_id != -1:                        
                        #update last appended value with t/2                        
                        a_d[a_d_id] += last_off

            a_ds.append(a_d)

        return a_ds

    """
    Compute distances over several offsets.
    """
    def discreteTransform(self, positions):        
        #pixel image for selected offsets
        d = []
        mdiv = int(self.bpp / 1)
        for offset in range(self.bpp):
            if offset % mdiv == 0:                
                pi = self.pixelImage(positions, offset)                
                #conversion of pixel image into distance between bright pixels
                ds = self.imageToDist(pi)                
                #convert distances based on triangular distribution
                stretched_d = self.distWithStretch(ds)
                #autonoise it
                a_d = self.autonoiseDistance(stretched_d)
                #self.storeDist(sum(a_d, []))
                #linearize
                d += sum(a_d, [])
            
        return d

    def storeDist(self, ds):
        f = 25000*[0]


        o = open("pi_ds.txt","w")
        
        for d in ds:            
            o.write("%.2f\n" % d)
            if d < 25000:
                f[d] += 1
        o.close()

        o = open("pi_real_dist.csv","w")
        for i in range(0,25000):
            o.write("%d\t%.8f\n" % (i, f[i]))
        o.close()


    def storePixelDist(self, pi):
        f = 1000*[0]

        for chr in pi:
            c_len = len(chr)

            ones = []
            for i in range(c_len):
                if chr[i] == 1:
                    ones.append(i)                    
        
            for i in range(len(ones)-1):
                d = ones[i+1] - ones[i]
                if d < 1000:
                    f[d] += 1

        o = open("pixel_dist2.csv","w")
        for i in range(0,1000):
            o.write("%d\t%d\n" % (i, f[i]))
        o.close()

        