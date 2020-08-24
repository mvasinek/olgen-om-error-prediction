from genome import OMGenome

import math
import random

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
