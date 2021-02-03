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
class PhasedResolution(OMGenome):
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

        if pos == None:
            self.d = self.discreteTransform(self.parsePositions())
        else:
            self.d = self.discreteTransform(pos)
        
        self.dlen = len(self.d)
        self.max_d = self.maxDistance(self.d)
        self.f0 = self.distF(self.d, self.max_d)
        self.p0 = self.distP(self.f0, self.max_d, self.dlen)
        
        self.p0ti = self.ebppTriangleInterpolation()        

    """
    Bitmap for phased image.
    """
    def phasedImage(self, off):
        phased_i = []
        
        for l in self.chrlen:
            phased_i.append((math.ceil((l + off) / self.bpp) + 1)*[0])
            
        return phased_i

    """
    Creates pixel image from positions in a sequence.
    """
    def pixelImage(self, pos, off):
        c_id = 0
        pi = []
        #proc???
        bo = self.bpp - off
        #bo = off
        pi_list = []
        for chromosome in pos:
            pi_list.append([])
            c_len = self.chrlen[c_id]

            pi_chr = (math.ceil((c_len + off) / self.bpp) + 1)*[0]

            i = 0
            for p in chromosome:
                x = (p+bo) // self.bpp
                if pi_chr[x] == 0:
                    pi_list[c_id].append(i)
                
                pi_chr[x] = 1               

                i += 1

            pi.append(pi_chr)

            c_id += 1            

        return (pi, pi_list)

    """
    Converts pixel map into distribution of bright pixel distaces.
    """
    def pixelToDist(self, pi, phased_i):
        d = []
        for c_id in range(len(pi)):
            chromosome = pi[c_id]
            phased = phased_i[c_id]
            
            #find first non/zero entry
            last_pos = 0
            last_ipos = 0
            for i in range(len(chromosome)):
                if chromosome[i] != 0:
                    last_ipos = i
                    #if pixel is not phase shifted then store i
                    if phased[i] == 0:
                        last_pos = i
                    #if pixel is phased shifted then store i + 0.5
                    else:
                        last_pos = i + 0.5
                    break

            #compute all distances
            for i in range(last_ipos+1, len(chromosome)):
                if chromosome[i] != 0:                    
                    daux = None
                    #not a phase shift d*2 => even values of d are without phase shit or with double phase shift, ensure it is integer.
                    if phased[i] == 0:
                        daux = 2*(i - last_pos)
                        last_pos = i
                    #phase shift
                    else:
                        daux = 2*(i + 0.5 - last_pos)
                        last_pos = i + 0.5
                    
                    d.append(int(daux))

        return d

    """
    Detect a sequence of consecutive bright pixels and replace by an average.
    """
    def detectSequence(self, chromosome, phased, i):
        gap = 0
        j = i+1
        last_one = i
        while j < len(chromosome) and gap < 2:
            if chromosome[j] == 0:
                gap += 1
            else:
                gap = 0
                last_one = j
                
            j += 1

        for k in range(i,last_one+1):
            if k < len(chromosome):
                chromosome[k] = 0

        #set l on average        
        d = last_one - i
        i_dh = i + d // 2
        """Nahodna verze
        if d % 2 == 0:
            chromosome[i_dh] = 1
        else:
            r = random.random()        
            if r > 0.5:                
                chromosome[i_dh + 1] = 0
        """
                
        chromosome[i_dh] = 1
        if d % 2 == 0:
            phased[i_dh] = 0
        else:
            phased[i_dh] = 1
                
        return last_one

    """
    Performs averaging procedure over all positions in image.
    """
    def averageRuns(self, pi, phased_i):
        c_id = 0
        for c_id in range(len(pi)):
            #for chromosome in pi:
            chromosome = pi[c_id]
            phased_chr = phased_i[c_id]
            i = 0
            cl = len(chromosome)-1
            while i < cl:
                if chromosome[i] != 0:
                    i = self.detectSequence(chromosome, phased_chr, i)
                    
                i += 1

        return pi

    """
    Compute distances over several offsets.
    """
    def discreteTransform(self, positions):        
        #pixel image for selected offsets
        d = []
        mdiv = int(self.bpp / 1)
        for offset in range(self.bpp):
            if offset % mdiv == 0:
                t0 = time.time()
                pi, pi_list = self.pixelImage(positions, offset)
                t1 = time.time()
                #print("pixel image in: ", (t1-t0))
                phased_i = self.phasedImage(offset)
                t2 = time.time()
                #print("phased image in: ", (t2-t1))
                ave_pi = self.averageRuns(pi, phased_i)
                #self.storePixelDist(ave_pi)
                t3 = time.time()
                #print("average image in: ", (t3-t2))
                d_aux = self.pixelToDist(ave_pi, phased_i)
                #self.storeDist(d_aux)
                #d += self.pixelToDist(ave_pi, phased_i)
                d += d_aux
                t4 = time.time()
                #print("to dist in: ", (t4-t3))
            
        return d

    """
    Triangle distribution interpretation of inter pixel distances.
    """
    def ebppTriangleInterpolation(self):
        dsize = (self.max_d+1)*self.ebpp
        if dsize < 25000:
            dsize = 25000
            
        p0ti = [0]*dsize

        denom = (self.ebpp-1)*(self.ebpp)+self.ebpp

        for i in range(0,self.max_d+1):
            d = 1
            #no phase shift
            if i % 2 == 0:                
                #former distance
                i0 = int(i/2) - 1
                for j in range(i0*self.ebpp + 1, (i0+2)*self.ebpp):
                    if j < (i0+1)*self.ebpp:
                        p0ti[j] += self.p0[i]*d/denom
                        d += 1
                    else:
                        p0ti[j] += self.p0[i]*d/denom
                        d -= 1
            else:
                i0 = int(i/2) - 1
                j_start = i0 * self.ebpp + 1 - int(self.ebpp/2)
                j_end = (i0+2)*self.ebpp - int(self.ebpp/2)
                for j in range(j_start, j_end):
                    if j < (i0+1)*self.ebpp - int(self.ebpp/2):
                        p0ti[j] += self.p0[i]*d/denom
                        d += 1
                    else:
                        p0ti[j] += self.p0[i]*d/denom
                        d -= 1
                
        return p0ti

    """
    Scale by ratio between two bpp. - Not used now.
    """
    def scale(self, p0, bpp0, bpp1):
        p1 = [0]*25001
        r = bpp1/bpp0
        
        for i in range(len(p0)):
            d1 = int(r*i)
            if d1 < 25001:
                p1[d1] += p0[i]

        return p1

    def storeDist(self, ds):
        f = 1000*[0]

        for d in ds:
            if d < 1000:
                f[d] += 1

        o = open("pixel_dist-ptd.csv","w")
        for i in range(0,1000):
            o.write("%d\t%d\n" % (i, f[i]))
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

        