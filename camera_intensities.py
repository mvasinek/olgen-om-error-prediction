import math
import time

from scipy.stats import norm

from genome import OMGenome


class CameraIntensities(OMGenome):
    
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

class DoubleSlit:
    @staticmethod
    def pdf(d,b,x,k=1):
        return (math.sin(b*x*k)/(b*x*k))**2 * (math.cos(k*d*x))**2