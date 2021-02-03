from genome import OMGenome

import math
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

class EnzymeMissing(MissingSite):
    """
    r - enzyme range
    v - missing site ratio, i.e. we expect 10% sites to be missed then v=0.1
    """
    def __init__(self, gpath, r, v):
        if gpath != None:
            super().__init__(gpath, 0)

        self.r = r
        self.v = v
        self.v_scaled = 0

        #self.v_scaled = self._correction()
        #print(self.v_scaled)

    @staticmethod
    def locate(chr_pos, r):
        if r < chr_pos[0]:
            return (-1,0)
        elif r > chr_pos[-1]:
            return (len(chr_pos), -1)

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
        tot = self.labelsCount()
        eMiss = tot*self.v

        i = 0
        chr_ids = list(range(0,24))

        while i < eMiss:
            #select chr
            chr_id = random.choices(chr_ids, self.chrlen,k=1)[0]
            
            r = random.randint(1,self.chrlen[chr_id]-1)

            #binary search
            dl_idx, dr_idx = EnzymeMissing.locate(self.pos[chr_id], r)

            if dl_idx == -1:
                #self.pos[chr_id].remove(self.pos[chr_id][0])
                del self.pos[chr_id][0]
                #self.pos[chr_id].pop(0)
            elif dr_idx == -1:
                #self.pos[chr_id].remove(self.pos[chr_id][-1])
                del self.pos[chr_id][-1]
                #self.pos[chr_id].pop(-1)
            else:
                dl = r - self.pos[chr_id][dl_idx]
                dr = self.pos[chr_id][dr_idx] - r

                if dl < dr:
                    #odstran site na pozici j-1
                    #self.pos[chr_id].remove(self.pos[chr_id][dl_idx])
                    del self.pos[chr_id][dl_idx]
                    #self.pos[chr_id].pop(dl_idx)
                else:
                    #self.pos[chr_id].remove(self.pos[chr_id][dr_idx])
                    del self.pos[chr_id][dr_idx]
                    #self.pos[chr_id].pop(dr_idx)

            """
            Working but slow
            #najdi nejblizsi a odstran jej
            for j in range(len(self.pos[chr_id])):
                if self.pos[chr_id][j] > r:
                    #j-1 pozice prvniho mensiho
                    #j pozice prniho vetsiho
                    dl = 0
                    if j-1 >= 0:
                        dl = r - self.pos[chr_id][j-1]

                    dr = self.pos[chr_id][-1]
                    if j < self.chrlen[chr_id]:
                        dr = self.pos[chr_id][j]

                    if dl < dr:
                        #odstran site na pozici j-1
                        self.pos[chr_id].remove(self.pos[chr_id][j-1])
                    else:
                        self.pos[chr_id].remove(self.pos[chr_id][j])

                    break
            """

            i += 1

        queue.put(self.pos)

    def correction(self):
        #compute p(average distance) distribution
        max_d = 4*self.r
        p_ave_d = [0]*(4*self.r + 1)

        tot = 0

        for chr_p in self.pos:
            tot += len(chr_p)
            if len(chr_p) == 0:
                break

            if len(chr_p) == 1:
                p_ave_d[max_d] += 1
                break

            for i in range(0, len(chr_p)):
                ave_dist = 0

                if i == 0:
                    ave_dist = chr_p[i+1]-chr_p[i] + 2*self.r
                elif i == len(chr_p) - 1:
                    ave_dist = chr_p[i]-chr_p[i-1] + 2*self.r
                else:
                    ave_dist = chr_p[i]-chr_p[i-1] + chr_p[i+1]-chr_p[i]

                if ave_dist > 4*self.r:
                    ave_dist = 4*self.r

                p_ave_d[int(ave_dist)] += 1

        #compute scale factor
        scale = 0
        #f = open("p_ave_d.csv","w")
        for i in range(4*self.r+1):
            p_ave_d[i] /= tot
            #f.write("%d\t%.8f\t%.8f\t%.8f\n" % (i, p_ave_d[i], i/(4*self.r), p_ave_d[i]*(i/(4*self.r))))
            scale += p_ave_d[i]*(i/(4*self.r))
        #f.close()

        #print("scale: ", self.v / scale)

        self.v_scaled = self.v / scale


    def storePD(self):
        f = open("pd.csv","w")

        for digested_part in range(4*self.r):
            missed = 4*self.r - digested_part

            p_miss = self.v_scaled * missed / (4*self.r) #p_digested    
            p_digested = 1-p_miss

            f.write("%d\t%.8f\n" % (int(digested_part/2), p_digested))

        f.close()

    def digested(self, dl, dr):
        if dl > 2*self.r:
            dl = 2*self.r
        if dr > 2*self.r:
            dr = 2*self.r

        sumd = (4*self.r) - (dl + dr)/2

        p_miss = self.v_scaled * sumd / (4*self.r) #p_digested
        p_d = 1 - p_miss
        rd = random.random()

        if rd < p_d:
            return True

        return False

#em = EnzymeMissing(sys.argv[1], 1000, 0.1)
#em.storePD()
