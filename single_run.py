import math
import os
import sys
import time

from multiprocessing import Queue, Process

from genome import OMGenome

from missing_site import MissingSite
from insertion_site import InsertionSite
from phased_resolution import PhasedResolution


"""
Class responsible for performing a single simulation.
"""
class SingleRun(OMGenome):
    def __init__(self, gpath, args, irate, iterations, x0, y0, x1, y1, scale, resolution=450, miss_konst=1.0, ebpp=500):
        self.genome_path = gpath
        self.args = args
        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        self.irate = irate
        self.iterations = iterations

        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.scale = scale
        self.resolution = resolution
        self.miss_k = miss_konst
        self.ebpp = ebpp

        self.p0ti = self.processFull()

    """
    Simulation is performed self.iterations times, each by separate process.
    """
    def processFull(self):
        pt = None

        queue = Queue()

        procs = [Process(target=self.processOne, args=(queue,)) for i in range(self.iterations)]

        for p in procs:
            p.start()

        results = [queue.get() for p in procs]

        for p in procs:
            p.join()

        for ti in results:
            if pt == None:
                pt = ti                
            else:
                for i in range(25000):
                    pt[i] += ti[i]

        for i in range(25000):
            pt[i] = pt[i] / self.iterations

        return pt

    """
    Consecutively applies sequence of three transformations on the reference genome distributions.
    """
    def processOne(self, pqueue):
        #proved nad pozicemi missing site
        mi = MissingSite(None, self.miss_k)
        mi.pos = self.pos
        mi.chrlen = self.chrlen
        mi.setPD(mi.bezierPD(self.x0, self.y0, self.x1, self.y1, self.scale))
        
        queue = Queue()
        p = Process(target = mi.processOne, args=(queue,))
        p.start()        
        
        pos_mi = queue.get()        
        p.join()

        #proved insertion
        ins = InsertionSite(None,self.chrlen, pos_mi)
        ins.pos = pos_mi
        ins.expectedInsertions(self.irate)
        
        p = Process(target = ins.processOne, args=(queue,))
        p.start()
        
        pos_ins = queue.get()        
        p.join()
        
        #proved resolution
        gg = PhasedResolution(self.genome_path, pos_ins, self.chrlen, int(self.resolution), int(self.ebpp))
        pqueue.put(gg.p0ti)

