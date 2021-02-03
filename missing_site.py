from genome import OMGenome
from multiprocessing import Process

import math
import os
import random
import sys

sign = lambda x: math.copysign(1, x)

"""
Functions responsible for Bezier function computation.
"""
def b(x,y):
    ts = [i/100 for i in range(101)]
    r = []
    for t in ts:
        px = (1-t)**3*x[0] + 3*t*(1-t)**2*x[1] + 3*(1-t)*t**2*x[2] + t**3*x[3]
        py = (1-t)**3*y[0] + 3*t*(1-t)**2*y[1] + 3*(1-t)*t**2*y[2] + t**3*y[3]

        r.append((px,py))
    return r

def getOnT(x,y,ts):
    r = []
    for t in ts:
        px = (1-t)**3*x[0] + 3*t*(1-t)**2*x[1] + 3*(1-t)*t**2*x[2] + t**3*x[3]
        py = (1-t)**3*y[0] + 3*t*(1-t)**2*y[1] + 3*(1-t)*t**2*y[2] + t**3*y[3]

        r.append((px,py))
    return r

def bWithScale(x,y,xScale):
    ts = [i/100 for i in range(101)]

    r = []
    for t in ts:
        px = ((1-t)**3*x[0] + 3*t*(1-t)**2*x[1] + 3*(1-t)*t**2*x[2] + t**3*x[3])*xScale
        py = (1-t)**3*y[0] + 3*t*(1-t)**2*y[1] + 3*(1-t)*t**2*y[2] + t**3*y[3]

        r.append((px,py))
    return r

def getOnTwithScale(x,y,ts, xScale):
    r = []
    for t in ts:
        px = ((1-t)**3*x[0] + 3*t*(1-t)**2*x[1] + 3*(1-t)*t**2*x[2] + t**3*x[3])*xScale
        py = (1-t)**3*y[0] + 3*t*(1-t)**2*y[1] + 3*(1-t)*t**2*y[2] + t**3*y[3]

        r.append((px,py))
    return r

def bezierCoefficients(a,b,c,d):
    return (-a + 3*b - 3*c + d, 3*a -6*b + 3*c, -3*a + 3*b, a)

def cubicRoots(p1,p2,p3,p4):
    A = p2 / p1
    B = p3 / p1
    C = p4 / p1

    A3 = A/3
    Q = B/3 - A3*A3
    R = -A3*A3*A3 + (A3*B-C)/2
    Q3 = Q*Q*Q
    D = Q3 + R*R

    t = []

    if(D >= 0):
        D2 = math.sqrt(D)
        S = sign(R + D2) * abs(R+D2)**(1/3)
        T = sign(R - D2) * abs(R-D2)**(1/3)
        t.append(-A3 + (S+T))
        
    else:
        th = math.acos( R / math.sqrt(-Q3))
        Q2 = 2* math.sqrt(-Q)

        t.append(Q2 + math.cos(th/3) - A3)
        t.append(Q2 + math.cos((th+2*math.pi)/3) - A3)
        t.append(Q2 + math.cos((th+4*math.pi)/3) - A3)

    return list(filter(lambda x: x >= 0 and x <= 1.0, t))


"""
Intersection of Bezier function with particular fragment length.
"""
def findIntersection(points, line):
    bx = bezierCoefficients(points[0][0], points[1][0], points[2][0], points[3][0])
    by = bezierCoefficients(points[0][1], points[1][1], points[2][1], points[3][1])

    dx = line[0][0] - line[1][0]
    dy = line[1][1] - line[0][1]

    C = line[0][0] * (-1)*dy + line[0][1]*(-1)*dx

    q1 = dy * bx[0] + dx * by[0]
    q2 = dy * bx[1] + dx * by[1]
    q3 = dy * bx[2] + dx * by[2]
    q4 = dy * bx[3] + dx * by[3] + C

    roots = cubicRoots(q1, q2, q3, q4)
    points = []
    for root in roots:
        t = root
        t2 = t*t
        t3 = t*t2

        x = bx[0]*t3 + bx[1]*t2 + bx[2]*t + bx[3]
        y = by[0]*t3 + by[1]*t2 + by[2]*t + by[3]

        if dx != 0.0:
            s = x - line[0][0] / -dx
        else:
            s = y - line[0][1] / dy

        if line[0][0] - line[1][0] != 0:
            s = (x-line[0][0]) / (line[1][0]-line[0][0])
        else:
            s = (y-line[0][1]) / (line[1][1]-line[0][1])

        if t<0 or t>1 or s<0 or s>1:
            pass
        else:
            points.append((x,y))

    return points

"""
Class responsible for introducing missing site error by removing restriction sites by digestion probability.
"""
class MissingSite(OMGenome):
    def __init__(self, gpath, miss_k):
        if gpath != None:
            super().__init__(gpath)

        self.lim = 25000
        self._pd = None

    """
    Sets particular digestion probability.
    Function can be used to enter arbitrary pd as an array of values for all fragment lengths up to self.lim.
    """
    def setPD(self, pd):
        self._pd = pd

    """
    Computatition of PD based on parameterized Bezier function.
    """
    def bezierPD(self, x0, y0, x1, y1, scale):
        x = (0,1,0,1)
        y = (0,0,1,1)

        points = ((0,0),(x0,y0),(x1,y1),(1,1))        
        
        pd = self.lim*[1]
        pd[0] = 0

        for d in range(1,int(scale)):
            df = d/scale
            line = ((df,0),(df,1))
            fi = findIntersection(points, line)
            if len(fi) > 0:
                pd[d] = fi[0][1]            

        return pd

    """
    Digestion probability given as an average of two neighbouring fragmnet lengths.
    """
    def pd(self, d1, d2):
        davg = round((d1+d2)/2)

        if davg >= self.lim:
            return 1        
        
        return self._pd[davg]

    """
    Returns True/False if the site is digested or not
    """
    def digested(self, d1, d2):
        p = self.pd(d1, d2)

        r = random.random()
        if r <= p:
            return True
        return False    

    """
    Processes missing site transformation.
    """
    def processOne(self, queue):
        p1 = []

        chr_id = 0
        for ps in self.pos:
            p1_c = []
            l = len(ps)

            #handle first point
            d_l = ps[0]
            d_r = ps[1]-ps[0]

            if self.digested(d_l, d_r):
                p1_c.append(ps[0])

            for i in range(1,l-1):
                d_l = d_r
                d_r = ps[i+1] - ps[i]

                if self.digested(d_l, d_r):
                    p1_c.append(ps[i])

            #handle last point
            d_l = d_r
            d_r = self.chrlen[chr_id] - ps[l-1]

            if self.digested(d_l, d_r):
                p1_c.append(ps[l-1])

            chr_id += 1
            p1.append(p1_c)

        queue.put(p1)

def LoadPds(wrkdir, iterations, step):   
    if not os.path.exists(wrkdir):
        return None

    pds = {}
    for i in range(1, iterations+1):
        mr = float("%.2f" % (step*i))
        fpath = wrkdir + "/" + "pd-%.2f.csv" % (mr)

        pd = [0]*2001
        f = open(fpath, "r")
        lines = f.readlines()

        for line in lines[1:]:
            l = line.strip()
            if len(l) == 0:
                continue

            v = l.split(';')

            pd[int(v[0])] = float(v[1])

        pds[mr] = pd

    return pds

class EnzymeMissing(MissingSite):
    pds = LoadPds("pd100", 50, 0.01)

    """
    r - enzyme range
    v - missing site ratio, i.e. we expect 10% sites to be missed then v=0.1
    """
    def __init__(self, gpath, r, v):
        if gpath != None:
            super().__init__(gpath, 0)

        self.r = r
        self.v = v

    @staticmethod
    def ReloadPDS():
        EnzymeMissing.pds = LoadPds("pd100", 50, 0.01)

    @staticmethod
    def locate(chr_pos, r):
        if r < chr_pos[0]:
            return (-1,0)
        elif r > chr_pos[-1]:
            return (len(chr_pos)-1, -1)

        left = 0
        right = len(chr_pos) - 1
        
        while left + 1 != right and right > 0 and left < len(chr_pos)-1:
            mid = int((left+right)/2)
            mid_l = mid
            mid_r = mid+1

            if chr_pos[mid_l] <= r and r <= chr_pos[mid_r]:
                return (mid_l, mid_r)

            if chr_pos[mid_l] < r:
                left = mid_r
            elif chr_pos[mid_r] > r:
                right = mid_l

        return (left, right)

    def processOne(self, queue):  
        p1 = []
        for _ in range(len(self.pos)):
            p1.append([])

        chr_id = 0
        for sites in self.pos:
            for i in range(len(sites)):
                davg = int(self.__dAvg(sites, i))
                digest_prob = EnzymeMissing.pds[self.v][davg]
                rn = random.random()

                if rn <= digest_prob:
                    #digested
                    p1[chr_id].append(sites[i])

            chr_id += 1

        return p1
        
    def __dAvg(self, sites, site_id):
        if site_id == 0:
            dl = 2*self.r
        else:
            dl = sites[site_id] - sites[site_id-1]
        if site_id+1 == len(sites):
            dr = 2*self.r
        else:
            dr = sites[site_id+1] - sites[site_id]
                    
        davg = 2*self.r
        if dl <= 2*self.r and dr <= 2*self.r:
            davg = (dl+dr)/2
        elif dl <= 2*self.r:
            davg = (dl+2*self.r)/2
        elif dr <= 2*self.r:
            davg = (dr+2*self.r)/2

        return davg    
    
    def experimentalProcessOne(self):
        tot = self.labelsCount()
        eMiss = tot*self.v
        eDig = int(tot-eMiss)
        digested = []
        for i in range(24):
            digested.append({})

        f_d = [0]*(2*self.r+1)

        i = 0
        dig = 0
        chr_ids = list(range(0,24))
        chr_id_random = random.choices(chr_ids, self.chrlen,k=int(eDig+1))
        while dig < eDig:
            #select chr
            chr_id = chr_id_random[i % (eDig+1)]
            #select binding position
            r = random.randint(1,self.chrlen[chr_id]-1)

            #binary search
            dl_idx, dr_idx = EnzymeMissing.locate(self.pos[chr_id], r)

            if dl_idx == -1:
                if self.pos[chr_id][0] - self.r > r and not 0 in digested[chr_id]:
                    digested[chr_id][0] = True
                    dig += 1       
            elif dr_idx == -1:
                if self.pos[chr_id][-1] + self.r > r and not len(self.pos[chr_id])-1 in digested[chr_id]:
                    digested[chr_id][len(self.pos)-1] = True
                    dig += 1                    
            else:
                dl = r - self.pos[chr_id][dl_idx]
                dr = self.pos[chr_id][dr_idx] - r

                if dl < dr:
                    #odstran site na pozici j-1
                    if dl < self.r and not dl_idx in digested[chr_id]:
                        digested[chr_id][dl_idx] = True
                        davg = self.__dAvg(self.pos[chr_id], dl_idx)
                        f_d[int(davg)] += 1                        
                        dig += 1                    
                else:
                    if dr < self.r and not dr_idx in digested[chr_id]:
                        digested[chr_id][dr_idx] = True
                        davg = self.__dAvg(self.pos[chr_id], dr_idx)                        
                        f_d[int(davg)] += 1
                        dig += 1     

            i += 1

        #print('digestion completed')
        f_all = [0]*(2*self.r+1)

        for sites in self.pos:
            for site_id in range(1,len(sites)-1):
                davg = int(self.__dAvg(sites, site_id))                
                f_all[davg] += 1

        pd = [0]*(2*self.r+1)
        for i in range(2*self.r+1):            
            if f_all[i] == 0:
                pd[i] = 0
            else:
                pd[i] = f_d[i] / f_all[i]

        return pd

    @staticmethod
    def NoLabelingPortion(sites, chrlen, erange):
        no_lab = 0

        chr_id = 0
        for chr_pos in sites:
            cur_pos = 0

            for pos in chr_pos:
                if pos - erange > cur_pos:
                    no_lab += (pos - erange) - cur_pos

                cur_pos = pos + erange

            #add the rest of chromosome
            if cur_pos < chrlen[chr_id]:
                no_lab += chrlen[chr_id] - cur_pos

            chr_id += 1            

        gsize = 0
        for l in chrlen:
            gsize += l

        return (no_lab, no_lab/gsize)

    @staticmethod
    def GeneratePDS(genome, iterations):
        try:
            os.mkdir("pd100")
        except:
            print("unable to create pd100 directory")
            exit(1)

        tot = 50
        cur = 0
        for i in range(1,51):
            mr = 0.01*i
            EnzymeMissing.GeneratePD(genome, mr, 1000, iterations)
            cur += 1            
            print("Generated %d/%d" % (cur,tot))

    @staticmethod
    def GeneratePD(genome, mr, e, iterations):
        pd = [0]*(2*e+1)
        for _ in range(iterations):
            em = EnzymeMissing(None, 1000, mr)
            em.pos = genome.pos
            em.chrlen = genome.chrlen

            one_pd = em.experimentalProcessOne()

            for i in range(len(one_pd)):
                pd[i] += one_pd[i]

        for i in range(len(pd)):
            pd[i] /= iterations

        #store pd
        fpath = "pd100/pd-%.2f.csv" % (mr)
        f = open(fpath,"w")
        f.write("Average Distance; Digestion Probability\n")
        for i in range(len(pd)):
            f.write("%d;%.5f\n" % (i, pd[i]))
        f.close()

        print("PD generation finished: %.2f" % mr)


"""
omg = OMGenome(sys.argv[1])
iterations = int(sys.argv[2])

ps = []
for i in range(1,51):
    mr = 0.01*i

    p = Process(target=EnzymeMissing.GeneratePD, args=(omg, mr, 1000, iterations))

    ps.append(p)

for p in ps:
    p.start()

for p in ps:
    p.join()
"""

"""
f = open("no-label-part-genome.csv", "w")
for e in range(1,201):
    erange = e*50

    f.write("%d;%3f\n" % (erange, EnzymeMissing.NoLabelingPortion(omg.pos, omg.chrlen, erange)[1]))

f.close()
"""



#em = EnzymeMissing(sys.argv[1], 1000, 0.1)
#em.storePD()