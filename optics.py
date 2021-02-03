import random

from molecule_generator import Molecule
from camera_intensities import *

def LoadCDFRange(start, end):
    cdfs = []
    for _ in range(end+1):
        cdfs.append([])

    for sigma in range(start, end+1):
        cdfs[sigma] = CameraIntensities.CDF(sigma)

    return cdfs

class Optics:
    def __init__(self, molecules, resolution, sizing_coef, pulse_sigma):
        self._mols = molecules
        self._res = resolution
        self._dtr = sizing_coef
        self._psi = pulse_sigma                

        #gaussian pulse
        self.cdf = CameraIntensities.CDF(pulse_sigma)
        self.orig_cdf = self.cdf

    """
    Gaussian diffusion by precomputed cdf(sigma).
    """
    def __precomputedDiffusionCDF(self, pi, p, x, res):        
        bpp = self._res
        
        mid_cdf = len(self.cdf) // 2           

        #intensity in point p
        try:
            pi[x] += (self.cdf[mid_cdf-res+bpp] - self.cdf[mid_cdf-res])
        except:
            print(len(self.cdf), mid_cdf-res+bpp, mid_cdf-res)        

        i = 1        
        while x - i >= 0 and mid_cdf-res-i*bpp >= 0:
            #intensity to the right                   
            pi[x-i] += (self.cdf[mid_cdf-res-(i-1)*bpp] - self.cdf[mid_cdf-res-i*bpp])
            i += 1
        
        i = 1
        while x+i < len(pi) and mid_cdf-res+bpp+i*bpp < len(self.cdf):
            #intensity to the right                        
            pi[x+i] += (self.cdf[mid_cdf-res+bpp+i*bpp] - self.cdf[mid_cdf-res+bpp+(i-1)*bpp])
            i += 1

    """
    Creates pixel image from positions in a sequence.
    """
    def _pixelImage(self, molecule):                
        pi_mol = (molecule.length // self._res + 1)*[0]

        for p in molecule.pos:
            x = p // self._res
            r = p % self._res

            self.__precomputedDiffusionCDF(pi_mol, p, x, r)   

        return pi_mol
   
    """
    Min Max variant of position detection
    """
    def _piToPosMinMax(self, pi):
        res = []
        pos = 1
        s = 0
        while pos < len(pi):
            if pi[pos] > 0:
                l = 0

                while pos + l < len(pi) and pi[pos+l] >= pi[pos + l -1]:
                    s += pi[pos+l]
                    l += 1

                while pos + l < len(pi) and pi[pos+l] > 0 and pi[pos+l] <= pi[pos + l - 1]:
                    s += pi[pos+l]
                    l += 1                
                
                h = (0.5 + (random.random())*self._dtr - self._dtr/2)*s

                cur = 0                
                site = 0
                    
                for i in range(l):
                    if cur + pi[pos+i] >= h:
                        rest = h - cur
                        site = (pos+i)*self._res + rest / pi[pos+i] * self._res
                        break
                    else:
                        cur += pi[pos + i]

                if site > 0:
                    res.append(site)

                pos += l
                if pi[pos-1] > 0:
                    s = pi[pos - 1]
                else:
                    s = 0
            else:
                pos += 1
                s = 0
        
        return res

    def processOne(self):
        for mol in self._mols:
            if len(mol.pos) == 0:
                continue

            #pixel image
            pi = self._pixelImage(mol)                        
            pos = self._piToPosMinMax(pi)            
            mol.pos = pos
            mol.length = len(pi)*375