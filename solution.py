import math
import os
import sys
import time
import copy

from multiprocessing import Queue, Process, freeze_support

from genome import OMGenome
from bnx import BNX, BNXData

from single_run import SingleRun, LaplaceSingleRun, NoCameraSingleRun, OnlyCorrectionSingleRun, InsertionCorrectionSingleRun, IndelsCorrectionSingleRun, BezierCorrectionSingleRun, SizingCorrectionSingleRun

from scipy.optimize import rosen, differential_evolution

"""
mean squared error between two lists of values a,b
lim - variable limiting the index in a,b
"""
def mse(a,b,lim):
    mse = 0
    for i in range(lim):
        mse += (a[i]-b[i])**2
    return mse/lim

"""
mean absolute error between two lists of values a,b
lim - variable limiting the index in a,b
"""
def mae(a,b,lim):
    mae = 0
    for i in range(lim):
        mae += abs(a[i]-b[i])
    return mae

def maeBins(a, b, lim):
    mae = 0

    for i in range(int(lim/25)):
        p_sum_a = 0
        p_sum_b = 0
        for j in range(25):
            p_sum_a += a[i*25 + j]
            p_sum_b += b[i*25 + j]

        mae += 25*abs(p_sum_a - p_sum_b)

    return mae

"""
Class responsible for start of simulation using differential evolution algorithm.
"""
class Solution:
    def __init__(self, genome_path, bnx_path, args):
        self._g_path = genome_path
        self._bnx_path = bnx_path
        self._args = args

        #directory to store output
        self._folder, self._bnx_name = self.__parsePath(bnx_path)
        
        self._limit = args.dist_max_len
        self._min_mae = 1000

        self._opt_x = None

        self._orig_genome = OMGenome(self._g_path)
        self._orig_pos = self._orig_genome.parsePositions()
        self._orig_chrlen = self._orig_genome.chrlen
        self._orig_ilen = self._orig_genome.intervalLength()

        if self._args.autonoised_genome == 'none':
            self._ag_genom = self._orig_genome
        else:
            self._ag_genom = OMGenome(self._args.autonoised_genome)

        self._bnxf = None

        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            self._bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            self._bnxf = BNX(self._bnx_path)

        self.bnxp = self._bnxf.p0
        self._bnx_f_emols = BNX.ExpectedMolsPerGenome(self._ag_genom.intervalLength(), self._bnxf.mol_dist)
        print("emols: %.4f" % self._bnx_f_emols)

    def bnx(self):
        return self._bnx_name

    """
    Returns tupple: (folder of input bnx, bnx file name)
    """
    def __parsePath(self, bnx):
        if "\\" in bnx:
            #windows path
            v = bnx.split("\\")
            return ("\\".join(v[:-1]) + "\\", v[-1])
        else:
            #linux
            v = bnx.split("/")
            return ("/".join(v[:-1]) + "/", v[-1])

    """
    Method first collects all arguments from argparse objects and initializes differential evoluion.
    Returns three objects with distributions of fragment lengths>
    - orig -  from reference genome
    - bnxf - fom bnx file
    - onea - from simulation based on optimal fit by differential evolution
    """
    def _collect(self):
        #original bnx
        #orig = OMGenome(self._g_path)
        orig = self._orig_genome

        bnxf = None        

        a = self._args

        x0 = (a.min_x0, a.max_x0)
        y0 = (a.min_y0, a.max_y0)

        x1 = (a.min_x1, a.max_x1)
        y1 = (a.min_y1, a.max_y1)

        v_miss = (0, 0.4)
        e_range = (1000, 1000)

        b_lim = (a.min_bez_lim, a.max_bez_lim)
        ir = (a.min_irate, a.max_irate)

        ebpp = (a.min_ebpp, a.max_ebpp)
        ud = (a.min_ud, a.max_ud)

        sigma = (a.min_sigma, a.max_sigma)

        det_lim = (a.min_det_thr, a.max_det_thr)
        det_sigma = (a.min_det_sigma, a.max_det_sigma)

        #if self._args.de_obpp:
            #obpp = (a.min_obpp, a.max_obpp)
        #    bounds = (x0,y0,x1,y1,b_lim,ir, ebpp, ud, obpp, sigma)
        #else:
        #    bounds = (x0,y0,x1,y1,b_lim,ir, ebpp, ud, sigma)

        cam_a = (0,100)
        cam_b = (0,100)
        cam_c = (0,100)
        #cam_k = (0,0.01)
        cam_k = (300,380)
        sf = (0.8,1.2)

        #bounds = (x0,y0,x1,y1,b_lim,ir,ud,sigma,det_lim,det_sigma)
        #bounds = (x0,y0,x1,y1,b_lim,ir,ud,sigma,ebpp,det_lim,cam_a, cam_b, cam_c)
        #bounds = (x0,y0,x1,y1,b_lim,ir,ud,sigma,det_lim,cam_k)
        #bounds = (x0,y0,x1,y1,b_lim,ir,ud,sigma,det_lim)
        bounds = (v_miss, ir,sigma,det_lim,cam_k, sf)
        #bounds = (v_miss, ir)
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        pos = copy.deepcopy(self._orig_pos)

        #if self._args.de_obpp:
        #    onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7],x[8],sigma=x[9])
        #else:
        #    onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7],sigma=x[8])
            
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],sigma=x[7])
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7],sigma=x[8])
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7],x[8],x[9])
        onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[1], self._args.iterations, 1000,x[0],x[2],375,x[3],x[4],x[5])
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[1], self._args.iterations, 1250,x[0],358,375,0.216,327.9)
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[1], self._args.iterations, 1250,x[0],358,375,0.269,341)
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, 0.05, self._args.iterations, 1250,0.1,358,375,0.269,341)
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, 0.1, self._args.iterations, 1000,0.3,373,375,0.287,329.939)
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7],x[8])
        #X0=[0.90943215,0.01514273]
#X1=[0.3950192,0.16554716]
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, 0.004, self._args.iterations, 0.90943215,0.01514273,0.39570192,0.16554716,1602,600,0.954,sigma=0.394)
        #onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, 0.039, self._args.iterations, 4991,0.025,255,375,0.226,0.003224)
        pos = None

        return (orig, self._bnxf, onea)

    """
    Fitness function for differential evolution. Fitness is evaluated by comparison of bnx fragment lengths distribution 
    with simulated one feeded by parameters given by differential evolution.

    Returns mean absolute error
    """
    def _fitness(self, x):
        #x0 = x[0]
        #y0 = x[1]
        #x1 = x[2]
        #y1 = x[3]
        #scale = x[4]
        #e_range = x[0]
        #v_miss = x[1]
        #fpr = x[2]              
        #miss_k = x[6]

        #miss_k = x[6]

        pos = copy.deepcopy(self._orig_pos)        

        #if self._args.de_obpp:
            #res_ebpp = x[8]
            #sigma = x[9]
            #onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols, self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,res,miss_k,res_ebpp, sigma)
        #else:
        #sigma = x[3]
        #res = x[8]
        #det_lim = x[4]
        #det_sigma = x[9]
        #cam_k = x[5]
        #cam_a = x[10]
        #cam_b = x[11]
        #cam_c = x[12]

        #miss_k = x[0]
        #fpr = x[1]

            #onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols,  self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,res,miss_k,500, sigma)
        #onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols,  self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,miss_k, sigma, det_lim, det_sigma)
        #onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols,  self._args, fpr, self._args.iterations, e_range, v_miss, sigma, 375, det_lim, cam_k)
        #onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols,  self._args, fpr, self._args.iterations, 1250, miss_k, 358, 375, 0.269,341)
        onea = SingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[1], self._args.iterations, 1000,x[0],x[2],375,x[3],x[4],x[5])
        #onea = SingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols,  self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,miss_k, sigma, det_lim)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:
            #if self._args.de_obpp:
            #    print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_path, x[4],x[5],x[6],x[7],x[9], mae_r, self._min_mae))
            #else:
            #    print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_path, x[4],x[5],x[6],x[7],x[8], mae_r, self._min_mae))
            print("%s\t %.3f\t %.3f\t %.3f\t %.3f\t%.3f\t%.6f\t%.3f\t%.3f" % (self._bnx_path, x[0],x[1],x[2],x[3],x[4],x[5], mae_r, self._min_mae))
            #print("%s\t%.3f\t%.3f\t%.3f\t%.3f" % (self._bnx_path, x[0],x[1], mae_r, self._min_mae))
            #print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_path, x[4],x[5],x[6],x[7],x[8], mae_r, self._min_mae))
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    """
    Stores parameters in a separate file.
    """
    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)

        f.write("Enzyme distance: %.3f" % self._opt_x[0])
        f.write("Missing rate: %.3f" % self._opt_x[1])
        #f.write("Insertion rate: %.3f" % self._opt_x[2])
        #f.write("Sizing error: %.3f" % self._opt_x[3])
        #f.write("Camera threshold: %.3f" % self._opt_x[4])
        #f.write("Sinc factor: %.6f" % self._opt_x[5])

        #f.write("Parameters of Bezier curve X0=[x0,y0] and X1=[x1,y1]\n")
        #f.write("X0=[%.8f,%.8f]\n" % (self._opt_x[0], self._opt_x[1]))
        #f.write("X1=[%.8f,%.8f]\n\n" % (self._opt_x[2], self._opt_x[3]))

        #f.write("Fragment lengths described by Bezier curve are smaller than: %.0f\n" % self._opt_x[4])
        #f.write("Insertion rate: %.4f\n" % self._opt_x[5])        
        #f.write("Digestion probability for fragment lengths greater than Bezier limit: %.4f\n" % self._opt_x[6])
        
        #f.write("Sigma = %.0f\n\n" % self._opt_x[7])

        #f.write("Light intensity threshold = %.4f\n\n" % self._opt_x[8])
        #f.write("Light intensity spread sigma = %.4f\n\n" % self._opt_x[9])
        
        f.write("Mean absolute error: %.3f" % err)
        
        f.close()


    """
    Stores fragment lengths distribution.
    """
    def _storeDistribution(self, results):
        f = open(self._folder + self._args.output + ".csv", "w")

        f.write("Fragment Length\tGenome\tBNX\tPrediction\n")

        orig = results[0].p0
        bnx = results[1].p0
        tipo = results[2].p0ti
        
        for i in range(self._limit):
            f.write("%d\t%.8f\t%.8f\t%.8f\n" % ((i+1), orig[i], bnx[i], tipo[i]))
            
        f.close()
        
    """
    Public method to call storage of distribution(optional) and parameters(mandatory)
    """    
    def store(self):
        results = self._collect()

        if self._args.store_dist:
            self._storeDistribution(results)

        #compute mean absolute error
        err = mae(results[2].p0ti, self.bnxp, self._limit)
        self._storeParameters(err)


class NoCameraSolution(Solution):
    def _collect(self):
        #original bnx
        orig = OMGenome(self._g_path)
        bnxf = None
        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        a = self._args

        x0 = (a.min_x0, a.max_x0)
        y0 = (a.min_y0, a.max_y0)

        x1 = (a.min_x1, a.max_x1)
        y1 = (a.min_y1, a.max_y1)

        b_lim = (a.min_bez_lim, a.max_bez_lim)
        ir = (a.min_irate, a.max_irate)
        
        ud = (a.min_ud, a.max_ud)

        sigma = (a.min_sigma, a.max_sigma)

        bounds = (x0,y0,x1,y1,b_lim,ir, ud, sigma)
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        onea = NoCameraSingleRun(self._g_path, self._bnx_path, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],sigma=x[7])
            
        return (orig, bnxf, onea)

    def _fitness(self, x):
        x0 = x[0]
        y0 = x[1]
        x1 = x[2]
        y1 = x[3]
        scale = x[4]
        fpr = x[5]      
        
        miss_k = x[6]

        
        sigma = x[7]
        onea = NoCameraSingleRun(self._g_path, self._bnx_path, self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,miss_k, sigma)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:            
            print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_name, x[4],x[5],x[6],x[7], mae_r, self._min_mae))
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)

        f.write("Parameters of Bezier curve X0=[x0,y0] and X1=[x1,y1]\n")
        f.write("X0=[%.8f,%.8f]\n" % (self._opt_x[0], self._opt_x[1]))
        f.write("X1=[%.8f,%.8f]\n\n" % (self._opt_x[2], self._opt_x[3]))

        f.write("Fragment lengths described by Bezier curve are smaller than: %.0f\n" % self._opt_x[4])
        f.write("Insertion rate: %.4f\n" % self._opt_x[5])
        
        f.write("Digestion probability for fragment lengths greater than Bezier limit: %.4f\n" % self._opt_x[6])
        
        f.write("Sigma = %.0f\n\n" % self._opt_x[7])

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class OnlyCorrectionSolution(Solution):
    def __init__(self, genome_path, bnx_path, args):
        super().__init__(genome_path, bnx_path,args)

    def _collect(self):
        #original bnx
        orig = OMGenome(self._g_path)
        bnxf = None
        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        onea = OnlyCorrectionSingleRun(self._g_path, self._bnx_path, self._args.iterations)
            
        return (orig, bnxf, onea)

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Only correction solution\n\n")
        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class InsertionCorrectionSolution(Solution):
    def __init__(self, genome_path, bnx_path, args):
        super().__init__(genome_path, bnx_path,args)

    def _collect(self):
        #original bnx
        orig = OMGenome(self._g_path)
        bnxf = None
        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        ir = (self._args.min_irate, self._args.max_irate)
                
        bounds = (ir, (0,1))
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        onea = InsertionCorrectionSingleRun(self._g_path, self._bnx_path, self._args.iterations, irate=x[0])
            
        return (orig, bnxf, onea)

    def _fitness(self, x):        
        fpr = x[0]      
        
        onea = InsertionCorrectionSingleRun(self._g_path, self._bnx_path, self._args.iterations, irate=fpr)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:            
            print("%s\t %.3f\t %.3f\t %.3f" % (self._bnx_name, fpr, mae_r, self._min_mae))
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Only correction solution\n\n")
        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)

        f.write("Insertion rate: %.3f\n\n" % self._opt_x[0])

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class IndelsCorrectionSolution(Solution):
    def _collect(self):
        #original bnx
        orig = OMGenome(self._g_path)
        bnxf = None
        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        a = self._args
      
        ir = (a.min_irate, a.max_irate)
        
        ud = (a.min_ud, a.max_ud)

        bounds = (ir, ud)
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        onea = IndelsCorrectionSingleRun(self._g_path, self._bnx_path, self._args, x[0], self._args.iterations, x[1])
            
        return (orig, bnxf, onea)

    def _fitness(self, x):        
        fpr = x[0]      
        
        miss_k = x[1]
        
        onea = IndelsCorrectionSingleRun(self._g_path, self._bnx_path, self._args, fpr, self._args.iterations, miss_k)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:            
            print("%s\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_name, x[0], x[1], mae_r, self._min_mae))
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)
        
        f.write("Insertion rate: %.4f\n" % self._opt_x[0])
        
        f.write("Digestion probability: %.4f\n" % self._opt_x[1])

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class BezierCorrectionSolution(Solution):
    def _collect(self):
        #original bnx
        orig = OMGenome(self._g_path)
        bnxf = None
        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        a = self._args

        x0 = (a.min_x0, a.max_x0)
        y0 = (a.min_y0, a.max_y0)

        x1 = (a.min_x1, a.max_x1)
        y1 = (a.min_y1, a.max_y1)

        b_lim = (a.min_bez_lim, a.max_bez_lim)
        ir = (a.min_irate, a.max_irate)
        
        ud = (a.min_ud, a.max_ud)

        bounds = (x0,y0,x1,y1,b_lim,ir, ud)
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        onea = BezierCorrectionSingleRun(self._g_path, self._bnx_path, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6])
            
        return (orig, bnxf, onea)

    def _fitness(self, x):
        x0 = x[0]
        y0 = x[1]
        x1 = x[2]
        y1 = x[3]
        scale = x[4]
        fpr = x[5]      
        
        miss_k = x[6]
        
        onea = BezierCorrectionSingleRun(self._g_path, self._bnx_path, self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,miss_k)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:            
            print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_name, x[4],x[5],x[6], mae_r, self._min_mae))
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)

        f.write("Parameters of Bezier curve X0=[x0,y0] and X1=[x1,y1]\n")
        f.write("X0=[%.8f,%.8f]\n" % (self._opt_x[0], self._opt_x[1]))
        f.write("X1=[%.8f,%.8f]\n\n" % (self._opt_x[2], self._opt_x[3]))

        f.write("Fragment lengths described by Bezier curve are smaller than: %.0f\n" % self._opt_x[4])
        f.write("Insertion rate: %.4f\n" % self._opt_x[5])
        
        f.write("Digestion probability for fragment lengths greater than Bezier limit: %.4f\n" % self._opt_x[6])
        
        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class IndelsSizingCorrectionSolution(Solution):
    def _collect(self):
        #original bnx
        orig = OMGenome(self._g_path)
        bnxf = None
        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx
            bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        a = self._args
      
        ir = (a.min_irate, a.max_irate)
        
        ud = (a.min_ud, a.max_ud)

        sigma = (a.min_sigma, a.max_sigma)

        bounds = (ir, ud, sigma)
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        onea = SizingCorrectionSingleRun(self._g_path, self._bnx_path, self._args, x[0], self._args.iterations, x[1], x[2])
            
        return (orig, bnxf, onea)

    def _fitness(self, x):        
        fpr = x[0]      
        
        miss_k = x[1]

        sigma = x[2]
        
        onea = SizingCorrectionSingleRun(self._g_path, self._bnx_path, self._args, fpr, self._args.iterations, miss_k, sigma)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:            
            print("%s\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_name, x[0], x[1], x[2], mae_r, self._min_mae))
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)
        
        f.write("Insertion rate: %.4f\n" % self._opt_x[0])
        
        f.write("Digestion probability: %.4f\n" % self._opt_x[1])

        f.write("Sigma: %.4f\n" % self._opt_x[2])

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class LaplaceSolution(Solution):
    def __init__(self, genome_path, bnx_path, args):
        super().__init__(genome_path, bnx_path,args)

        self.levels = None

    def _collect(self):
        #original bnx
        orig = self._orig_genome
        bnxf = None

        a = self._args

        x0 = (a.min_x0, a.max_x0)
        y0 = (a.min_y0, a.max_y0)

        x1 = (a.min_x1, a.max_x1)
        y1 = (a.min_y1, a.max_y1)

        b_lim = (a.min_bez_lim, a.max_bez_lim)
        ir = (a.min_irate, a.max_irate)

        ud = (a.min_ud, a.max_ud)

        #we must extend bounds by Laplace(m,b) coefficients - each level means 2 more variables
        #bounds for m 0.8;1.1
        #bounds for b 0.01;0.2

        upper_levels = []
        for l in range(self._args.laplace_levels - 1):
            upper_levels.append(l*self._args.laplace_step + self._args.laplace_step)

        self.levels = upper_levels

        bounds = (x0,y0,x1,y1,b_lim,ir, ud)

        #l levels yields 2*l new variables 
        for l in range(self._args.laplace_levels):
            bounds += ((0.8,1.1),(0.01,0.2))
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        pos = copy.deepcopy(self._orig_pos)

        onea = LaplaceSingleRun(self._g_path, self._ag_genom,pos, self._orig_chrlen, self._bnx_f_emols, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],mb=x[8:8+2*self._args.laplace_levels], levels=self.levels)
            
        return (orig, bnxf, onea)

    def _fitness(self, x):
        x0 = x[0]
        y0 = x[1]
        x1 = x[2]
        y1 = x[3]
        scale = x[4]
        fpr = x[5]              
        miss_k = x[6]        

        pos = copy.deepcopy(self._orig_pos)
    
        mb = x[7:7+2*self._args.laplace_levels]
        onea = LaplaceSingleRun(self._g_path, self._ag_genom, pos, self._orig_chrlen, self._bnx_f_emols, self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,miss_k, mb=mb, levels=self.levels)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:            
            print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_name, x[4],x[5],x[6],x[7], mae_r, self._min_mae))

            print("Laplace(mb): ", mb)
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r

    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)

        f.write("Parameters of Bezier curve X0=[x0,y0] and X1=[x1,y1]\n")
        f.write("X0=[%.8f,%.8f]\n" % (self._opt_x[0], self._opt_x[1]))
        f.write("X1=[%.8f,%.8f]\n\n" % (self._opt_x[2], self._opt_x[3]))

        f.write("Fragment lengths described by Bezier curve are smaller than: %.0f\n" % self._opt_x[4])
        f.write("Insertion rate: %.4f\n" % self._opt_x[5])        
        f.write("Digestion probability for fragment lengths greater than Bezier limit: %.4f\n" % self._opt_x[6])

        mb = self._opt_x[7:7+2*self._args.laplace_levels] 
        
        f.write("Laplace(m,b):")
        f.write(str(mb) + "\n\n")

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()

class LiSolution(Solution):
    def _collect(self):
        pass

    def _fitness(self, x):
        pass

    def _storeParameters(self, err):
        pass