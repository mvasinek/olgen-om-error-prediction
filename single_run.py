import math
import os
import sys
import time

from multiprocessing import Queue, Process

from genome import OMGenome

from missing_site import MissingSite, EnzymeMissing
from insertion_site import InsertionSite
from sizing_error import SizingError, SizingErrorLi, SizingErrorLaplace
from phased_resolution import PhasedResolution
from camera import Camera
from camera_intensities import CameraIntensities
from stretch_factor import StretchFactor


"""
Class responsible for performing a single simulation.
"""
class SingleRun(OMGenome):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, resolution=450, miss_konst=1.0, ebpp=500, sigma=100):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0,sigma=100):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, res, miss_konst=1.0,sigma=100):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0, sigma=100, cam_lim=0.3, cam_sigma=100):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0, sigma=100, res=375, cam_lim=0.3, cam_a=10,cam_b=30,cam_c=100):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0, sigma=100, res=375, cam_lim=0.3, k=0.01):
    #def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0, sigma=100, cam_lim=0.3):
    def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, e_range, v_miss, sigma=100, res=375, cam_lim=0.3, k=0.01, sf=1.0):
        self.genome_path = gpath
        self.ag_genome = ag_genome
        self.emols = emols
        self.args = args
        self.pos = pos
        self.chrlen = chr_len
        #self.chrlen = 24*[0]
        #self.pos = self.parsePositions()
        self.res = int(res)
        self.irate = irate
        self.iterations = iterations
        
        self.detection_limit = cam_lim
        #self.detection_sigma = cam_sigma
        #self.cam_a = cam_a
        #self.cam_b = cam_b
        #self.cam_c = cam_c
        self.k = k

        #self.x0 = x0
        #self.y0 = y0
        #self.x1 = x1
        #self.y1 = y1

        self.e_range = int(e_range)
        self.v_miss = v_miss

        #self.scale = scale
        #self.resolution = resolution
        #self.miss_k = miss_konst
        #self.ebpp = ebpp
        self.sigma = sigma
        self.sf = sf

        #print(self.detection_limit, self.detection_sigma)

        self.p0ti = self.processFull()

    """
    Simulation is performed self.iterations times, each by separate process.
    """
    def processFull(self):
        pt = None

        queue = Queue()

        procs = [Process(target=self.processOne, args=(queue,)) for i in range(self.iterations)]

        """
        for p in procs:
            p.start()

        results = [queue.get() for p in procs]

        for p in procs:
            p.join()
        """

        results = []
        processed = 0
        p_id = 0
        num_pids = self.args.procs
        pid_start = p_id
        pid_end = p_id + num_pids

        while processed < self.iterations:
            for i in range(pid_start, pid_end):                
                procs[i].start()

            for i in range(pid_start, pid_end):                
                results.append(queue.get())

            for i in range(pid_start, pid_end):                
                procs[i].join()

            processed += pid_end - pid_start

            pid_start = pid_end
            pid_end += num_pids

            if pid_end > self.iterations:
                pid_end = self.iterations

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
        #mi = MissingSite(None, self.miss_k)
        #mi.pos = self.pos
        #mi.chrlen = self.chrlen
        #mi.setPD(mi.bezierPD(self.x0, self.y0, self.x1, self.y1, self.scale))        

        mi = EnzymeMissing(None, self.e_range, self.v_miss)
        mi.pos = self.pos
        mi.chrlen = self.chrlen
        #mi.correction()        
        
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

        #apply molecule stretch factor variation
        stf = StretchFactor(None, self.chrlen, pos_ins, self.sf)
        stf.pos = pos_ins
        
        p = Process(target = stf.processOne, args=(queue,))
        p.start()
        pos_stf = queue.get()        
        p.join()


        #apply sizing error        
        if self.args.use_laplace_li:
            siz = SizingErrorLi(None, 0, self.chrlen, pos_ins)
        else:
            siz = SizingError(None, self.sigma, self.chrlen, pos_stf)

        siz.pos = pos_stf
        p = Process(target = siz.processOne, args=(queue,))
        p.start()        
        pos_siz = queue.get()        
        p.join()
        
        
        if self.args.disable_camera:
            #print('camera disabled')
            pqueue.put(siz.distFromPositions(pos_siz))
        else:
            #apply camera transformation
            #print('camera enabled')
            #gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, 375, 375)            
            #gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, 450, 450)
            #gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, self.res, self.res)
            #gg = PhasedResolution(self.genome_path, self.pos, self.chrlen, self.res, self.res)
            #gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, self.res, self.res)
            #gg = Camera(self.genome_path, pos_siz, self.chrlen, 375, 375)
            #gg = Camera(self.genome_path, pos_siz, self.chrlen, self.res, self.res)
            #gg = CameraIntensities(pos_siz, self.chrlen, 375, self.detection_limit, self.detection_sigma)
            #gg = CameraIntensities(pos_siz, self.chrlen, self.res, self.detection_limit, self.cam_a, self.cam_b, self.cam_c)
            gg = CameraIntensities(pos_siz, self.chrlen, self.res, self.detection_limit, self.k)
            #gg = CameraIntensities(pos_siz, self.chrlen, 375, self.detection_limit, 0)

            #if self.args.autonoised_genome == 'none':
            #    pqueue.put(self.correctionFromDistribution(gg.p0ti, self.emols, self))
            #else:
            #pti = self.correctionFromDistribution(gg.p0ti, self.emols, self.ag_genome)
            pqueue.put(gg.p0ti)
        
        
        #proved resolution
        #gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, int(self.res), int(self.res))        
        #gg = PhasedResolution(self.genome_path, self.pos, self.chrlen, int(self.res), int(self.res))        
        #pqueue.put(gg.p0ti)

class NoCameraSingleRun(SingleRun):
    def __init__(self, gpath, bnxpath, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0, sigma=100):
        self.genome_path = gpath
        self.bnx_path = bnxpath
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
        
        self.miss_k = miss_konst
        
        self.sigma = sigma

        self.p0ti = self.processFull()

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

        #apply sizing error        
        if self.args.use_laplace_li:
            siz = SizingErrorLi(None, 0, self.chrlen, pos_ins)
        else:
            siz = SizingError(None,self.sigma, self.chrlen, pos_ins)

        siz.pos = pos_ins        
        p = Process(target = siz.processOne, args=(queue,))
        p.start()        
        pos_siz = queue.get()        
        p.join()

        pqueue.put(siz.distFromPositions(pos_siz, True, self.bnx_path, self.genome_path))

class OnlyCorrectionSingleRun(SingleRun):
    def __init__(self, gpath, bnxpath, iterations):
        self.genome_path = gpath
        self.bnx_path = bnxpath

        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        self.iterations = iterations

        self.p0ti = self.processFull()

    def processOne(self, pqueue):
        pqueue.put(self.distFromPositions(self.pos, True, self.bnx_path, self.genome_path))        

class InsertionCorrectionSingleRun(SingleRun):
    def __init__(self, gpath, bnxpath, iterations,irate):
        self.genome_path = gpath
        self.bnx_path = bnxpath

        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        self.iterations = iterations
        self.irate = irate

        self.p0ti = self.processFull()

    def processOne(self, pqueue):
        #proved insertion
        queue = Queue()
        
        ins = InsertionSite(None,self.chrlen, self.pos)
        ins.pos = self.pos
        ins.expectedInsertions(self.irate)
        
        p = Process(target = ins.processOne, args=(queue,))
        p.start()
        pos_ins = queue.get()        
        p.join()

        pqueue.put(self.distFromPositions(pos_ins, True, self.bnx_path, self.genome_path))

class IndelsCorrectionSingleRun(SingleRun):
    def __init__(self, gpath, bnxpath, args, irate, iterations, miss_konst=1.0):
        self.genome_path = gpath
        self.bnx_path = bnxpath
        self.args = args
        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        self.irate = irate
        self.iterations = iterations
        
        self.miss_k = miss_konst

        self.p0ti = self.processFull()

    """
    Consecutively applies sequence of three transformations on the reference genome distributions.
    """
    def processOne(self, pqueue):
        #proved nad pozicemi missing site
        mi = MissingSite(None, self.miss_k)
        mi.pos = self.pos
        mi.chrlen = self.chrlen
        mi.setPD(25000*[self.miss_k])
        
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

        pqueue.put(self.distFromPositions(pos_ins, True, self.bnx_path, self.genome_path))

class BezierCorrectionSingleRun(SingleRun):
    def __init__(self, gpath, bnxpath, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0):
        self.genome_path = gpath
        self.bnx_path = bnxpath
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
        
        self.miss_k = miss_konst

        self.p0ti = self.processFull()

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

        pqueue.put(self.distFromPositions(pos_ins, True, self.bnx_path, self.genome_path))

class SizingCorrectionSingleRun(SingleRun):
    def __init__(self, gpath, bnxpath, args, irate, iterations, miss_konst=1.0, sigma=100):
        self.genome_path = gpath
        self.bnx_path = bnxpath
        self.args = args
        self.chrlen = 24*[0]
        self.pos = self.parsePositions()
        self.irate = irate
        self.iterations = iterations
        
        self.miss_k = miss_konst
        self.sigma = sigma

        self.p0ti = self.processFull()

    """
    Consecutively applies sequence of three transformations on the reference genome distributions.
    """
    def processOne(self, pqueue):
        #proved nad pozicemi missing site
        mi = MissingSite(None, self.miss_k)
        mi.pos = self.pos
        mi.chrlen = self.chrlen
        mi.setPD(25000*[self.miss_k])
        
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

        siz = SizingError(None,self.sigma, self.chrlen, pos_ins)

        siz.pos = pos_ins        
        p = Process(target = siz.processOne, args=(queue,))
        p.start()        
        pos_siz = queue.get()        
        p.join()

        pqueue.put(self.distFromPositions(pos_siz, True, self.bnx_path, self.genome_path))

class LaplaceSingleRun(SingleRun):
    def __init__(self, gpath, ag_genome, pos, chr_len, emols, args, irate, iterations, x0, y0, x1, y1, scale, miss_konst=1.0, mb=None, levels=None):
        self.genome_path = gpath
        self.ag_genome = ag_genome
        self.args = args
        self.chrlen = chr_len
        self.pos = pos
        self.irate = irate
        self.iterations = iterations
        self.emols = emols

        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.scale = scale        
        self.miss_k = miss_konst              
        self.mb = self.__linearToPair(mb)
        self.levels = levels

        self.p0ti = self.processFull()

    """
    Converts linear array of value into pairs.
    [1,2,3,4,5,6] -> [(1,2),(3,4),(5,6)]
    """
    def __linearToPair(self, a):
        p = []

        for i in range(int(len(a)/2)):
            p.append((a[i], a[i+1]))

        return p

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

        #apply sizing error                
        siz = SizingErrorLaplace(None, self.levels, self.mb, self.chrlen, pos_ins)
        siz.pos = pos_ins        
        
        p = Process(target = siz.processOne, args=(queue,))
        p.start()        
        pos_siz = queue.get()        
        p.join()
        
        if self.args.disable_camera:
            #print('camera disabled')
            pqueue.put(siz.distFromPositions(pos_siz))
        else:
            #apply camera transformation
            #print('camera enabled')
            gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, 375, 375)
            if self.args.autonoised_genome == 'none':
                pqueue.put(self.correctionFromDistribution(gg.p0ti, self.emols, self.ag_genome))
            else:
                pti = self.correctionFromDistribution(gg.p0ti, self.emols, self.ag_genome)
                pqueue.put(pti)

            #gg = PhasedResolution(self.genome_path, pos_siz, self.chrlen, int(self.resolution), int(self.ebpp))
            #pqueue.put(gg.p0ti)