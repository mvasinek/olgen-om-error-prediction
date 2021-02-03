import math
import time

from scipy.stats import norm

from genome import OMGenome


class CameraIntensities(OMGenome):
    #def __init__(self, pos=None, chrlen=None, resolution=375, lim=0.3, sigma=100):
    #def __init__(self, pos=None, chrlen=None, resolution=375, lim=0.3, cam_a=10, cam_b=30,cam_c=100):
    #def __init__(self, pos=None, chrlen=None, resolution=375, lim=0.3, k=0.01):
    def __init__(self, pos=None, chrlen=None, resolution=375, lim=0.3, k=0.01, l=0.01):
    
        #expected base pairs per pixel
        self.bpp = resolution
        self.lim = lim
        self.chrlen = chrlen      

        #self.cam_a = cam_a
        #self.cam_b = cam_b
        #self.cam_c = cam_c  
        #self.sigma = sigma
        self.cdf = CameraIntensities.CDF(k)
        #self.cdf = SinC2.cdf(k, 375*3)
        self.d = self.discreteTransform(pos)
        
        self.dlen = len(self.d)
        #self.max_d = self.maxDistance(self.d)
        #self.f0 = self.distF(self.d, self.max_d)
        keys, f0 = self.distF2(self.d)
        self.p0ti = self.distP2(keys, f0, self.dlen)
        #self.p0ti = self.distP2(self.f0, self.max_d, self.dlen)        

    """
    Gaussian diffusion.
    """
    def __gaussianDiffusion(self, pi, p, x, res):
        sigma = self.sigma
        bpp = self.bpp

        #intensity in point p
        pi[x] += norm.cdf(p-res+bpp, loc=p, scale=sigma) - norm.cdf(p-res, loc=p, scale=sigma)

        i = 1
        while(x - i >= 0):
            #intensity to the right
            pr = norm.cdf(p-res-(i-1)*bpp, loc=p, scale=sigma) - norm.cdf(p-res-i*bpp, loc=p, scale=sigma)
            pi[x-i] += pr

            if pr < 0.01:
                break

            i += 1
        
        i = 1
        while(x+i < len(pi)):
            #intensity to the right
            pr = norm.cdf(p-res+bpp+i*bpp, loc=p, scale=sigma) - norm.cdf(p-res+bpp+(i-1)*bpp, loc=p, scale=sigma)
            pi[x+i] += pr

            if pr < 0.01:
                break            

            i += 1

    """
    Gaussian diffusion by precomputed cdf(sigma).
    """
    def __precomputedDiffusionCDF(self, pi, p, x, res):        
        bpp = self.bpp
        
        mid_cdf = len(self.cdf) // 2        

        #intensity in point p        
        try:
            pi[x] += self.cdf[mid_cdf-res+bpp] - self.cdf[mid_cdf-res]
        except:
            print(len(self.cdf), mid_cdf-res+bpp, mid_cdf-res)


        i = 1
        #while x - i >= 0 and (mid_cdf-res-i*bpp >= 0) and (mid_cdf-res-(i-1)*bpp < len(self.cdf)):
        while x - i >= 0 and mid_cdf-res-i*bpp >= 0:
            #intensity to the right       
            #try:
            pi[x-i] += self.cdf[mid_cdf-res-(i-1)*bpp] - self.cdf[mid_cdf-res-i*bpp]
            #except:
                #print("x-i",len(pi),x-i,len(self.cdf), mid_cdf-res-(i-1)*bpp, mid_cdf-res-i*bpp)

            i += 1
        
        i = 1
        while x+i < len(pi) and mid_cdf-res+bpp+i*bpp < len(self.cdf):
            #intensity to the right            
            #try:
            pi[x+i] += self.cdf[mid_cdf-res+bpp+i*bpp] - self.cdf[mid_cdf-res+bpp+(i-1)*bpp]
            #except:
            #    print("x+i",len(pi),len(self.cdf), p-res+bpp+i*bpp, p-res+bpp+(i-1)*bpp)

            i += 1

    """
    Creates pixel image from positions in a sequence.
    """
    def _pixelImage(self, pos):
        c_id = 0
        pi = []        
        for chromosome in pos:            
            c_len = self.chrlen[c_id]

            pi_chr = (c_len // self.bpp + 1)*[0]

            for p in chromosome:
                x = p // self.bpp
                r = p % self.bpp

                #self.__gaussianDiffusion(pi_chr, p, x, r)
                #self.__gaussianDiffusionCDF(pi_chr, p, x, r)
                self.__precomputedDiffusionCDF(pi_chr, p, x, r)

            pi.append(pi_chr)

            c_id += 1            

        return pi

    def __addIntensity(self, pi, x, a, b):        
        ar = a/(a+b)
        br = b/(a+b)
        s = (self.cam_a + self.cam_b + self.cam_c)*2
        iss = [self.cam_a/s,self.cam_b/s,self.cam_c/s,self.cam_c/s,self.cam_b/s,self.cam_a/s]
        #iss = [10,20,20,10]

        start_x = x - 3
        for i in range(len(iss)):
            if start_x + i >= 0 and start_x + i < len(pi):
                pi[start_x + i] += br*iss[i]
            if start_x + i + 1 >= 0 and start_x + i + 1 < len(pi):
                pi[start_x + i + 1] += ar*iss[i]
            

    """
    Creates pixel image from positions in a sequence.
    """
    def _pixelImageHandTuned(self, pos):
        c_id = 0
        pi = []        
        for chromosome in pos:            
            c_len = self.chrlen[c_id]

            pi_chr = (c_len // self.bpp + 1)*[0]

            for p in chromosome:
                x = p // self.bpp
                a = p % self.bpp
                b = self.bpp - a

                self.__addIntensity(pi_chr, x, a, b)
                

            pi.append(pi_chr)

            c_id += 1            

        
        return pi

    """
    Applies limit to say what pixel is still bright and what is not
    """
    def _toBright(self, pi, lim):
        pi_bright = []

        for c_pi in pi:
            c_bright = []

            for i in range(len(c_pi)):
                if c_pi[i] > lim:
                    c_bright.append(1)
                else:
                    c_bright.append(0)

            pi_bright.append(c_bright)

        return pi_bright

    """
    Converts to positions deduced from pixel image.
    """
    def _piToPos(self, pi, pi_bright):
        new_pos = []

        for c_id in range(len(pi)):
            new_chr_pos = []
            c_pi = pi[c_id]
            c_br = pi_bright[c_id]
            
            pos = 0
            while pos < len(c_pi):
                if c_br[pos] == 1:
                    start = pos
                    end = pos
                    while end < (len(c_pi) - 1) and c_br[end+1] == 1:
                        end += 1

                    #<start;end> obsahuji svitici pixely + 1px okoli
                    if start > 0:
                        start -= 1
                    if end < len(c_pi)-1:
                        end += 1

                    tot = 0
                    for i in range(start, end+1):
                        tot += c_pi[i]

                    #tot je celkova intenzita napric pixely
                    #uloha najit pixel, kam padne prostredni hodnota                    
                    partial = 0
                    halve = tot/2
                    s_id = start
                    while True:
                        partial += c_pi[s_id]
                        if partial < halve:
                            s_id += 1
                        else:
                            break

                    #pod s_id mame pixel obsahujici prostredni bod
                    #uloha - spocitat res, ktery zbyva na prostredni pixel
                    res = halve - (partial - c_pi[s_id])
                    #dopocitat pozici v prostrednim pixelu
                    mid_pos = self.bpp * res / c_pi[s_id]
                    #uloha - spocitat pozici
                    site_pos = s_id*self.bpp + mid_pos
                    #pridame site pos do pole novych pozic odectenych z obrazku
                    new_chr_pos.append(int(site_pos))

                    pos = end + 1
                else:
                    pos += 1

            new_pos.append(new_chr_pos)

        return new_pos

    """
    Compute distances over several offsets.
    """
    def discreteTransform(self, positions):        
        #pixel image
        t0 = time.time()
        pi = self._pixelImage(positions)        
        #pi = self._pixelImageHandTuned(positions)
        t1 = time.time()
        #print("pixel image in: %.4f" % (t1-t0))
        b_pi = self._toBright(pi, self.lim)        
        t2 = time.time()
        #print("to bright in: %.4f" % (t2-t1))
        c_pos = self._piToPos(pi, b_pi)        
        t3 = time.time()
        #print("pi to pos: %.4f" % (t3-t2))
        d = self.convertToDist(c_pos)        
        t4 = time.time()
        #print("to dist in: %.4f" % (t4-t3))
            
        return d

    @staticmethod
    def CDF(sigma):        
        #three sigma rule
        cdf = []
        x = -3*sigma
        if x > -3*375:
            x = -3*375

        while x <= 0:
            cdf.append(norm.cdf(x, scale=sigma))
            x += 1
                
        cdf_id = len(cdf) - 2 #the point before the middle point located in 0 which is the last in cdf(currently)
        while cdf_id >= 0:
            cdf.append(1 - cdf[cdf_id])
            cdf_id -= 1

        return cdf

#single slit difraction
class SinC2:
    """
    (sin(kx)/x)^2
    """            

    @staticmethod
    def cdf(k, z):
        #first compute total over z points to left and z points to right from zero
        #may contain large integration error !!!
        #first add middle point
        total = SinC2.pdf(-0.5, k)
        for i in range(-z, 0):
            total += (SinC2.pdf(i - 0.5, k) + SinC2.pdf(i+0.5, k))/2            

        for i in range(1, z+1):
            total += (SinC2.pdf(i - 0.5, k) + SinC2.pdf(i+0.5, k))/2

        cdf_a = [0]*(2*z+1)

        scale = 1/total

        y = -1*z

        for i in range(2*z+1): 
            if i == 0:
                cdf_a[i] = scale*(SinC2.pdf(y - 0.5, k) + SinC2.pdf(y+0.5, k))/2
            else:
                cdf_a[i] = cdf_a[i-1] + scale*(SinC2.pdf(y - 0.5, k) + SinC2.pdf(y+0.5, k))/2

            y += 1

        return cdf_a
            
    @staticmethod
    def pdf(x, k):
        return (math.sin(k*x)/(k*x))**2

def sinc(x):
    return math.sin(x) / x

#Double slit difraction
class DoubleSlit:
    @staticmethod
    def cdf(xk, d, r, z):
        total = DoubleSlit.pdf(-0.5, xk, d, r)
        for i in range(-z, 0):
            total += (DoubleSlit.pdf(i - 0.5, xk, d, r) + DoubleSlit.pdf(i + 0.5, xk, d, r))/2

        for i in range(1, z+1):
            total += (DoubleSlit.pdf(i - 0.5, xk, d, r) + DoubleSlit.pdf(i + 0.5, xk, d, r))/2

        cdf_a = [0]*(2*z+1)

        scale = 1/total

        y = -1*z

        for i in range(2*z+1): 
            if i == 0:
                cdf_a[i] = scale*(DoubleSlit.pdf(y - 0.5, xk, d, r) + DoubleSlit.pdf(y + 0.5, xk, d, r))/2
            else:
                cdf_a[i] = cdf_a[i-1] + scale*(DoubleSlit.pdf(y - 0.5, xk, d, r) + DoubleSlit.pdf(y+0.5, xk, d, r))/2

            y += 1

        return cdf_a

    @staticmethod
    def pdf(fx, xk, d, r):
        return (math.cos(d*fx*xk)**2) * (sinc(d*fx/r)**2)
        
        

#print(SinC2.cdf(0.5,10))

#ci = CameraIntensities([[51]], [120], 10, 0.3, 0.5 )

