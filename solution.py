import math
import os
import sys
import time

from multiprocessing import Queue, Process, freeze_support

from genome import OMGenome
from bnx import BNX

from single_run import SingleRun

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
        orig = OMGenome(self._g_path)
        bnxf = BNX(self._bnx_path)

        self.bnxp = bnxf.p0

        a = self._args

        x0 = (a.min_x0, a.max_x0)
        y0 = (a.min_y0, a.max_y0)

        x1 = (a.min_x1, a.max_x1)
        y1 = (a.min_y1, a.max_y1)

        b_lim = (a.min_bez_lim, a.max_bez_lim)
        ir = (a.min_irate, a.max_irate)

        ebpp = (a.min_ebpp, a.max_ebpp)
        ud = (a.min_ud, a.max_ud)

        if self._args.de_obpp:
            obpp = (a.min_obpp, a.max_obpp)
            bounds = (x0,y0,x1,y1,b_lim,ir, ebpp, ud, obpp)
        else:
            bounds = (x0,y0,x1,y1,b_lim,ir, ebpp, ud)

        print(bounds)
        
        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        if self._args.de_obpp:
            onea = SingleRun(self._g_path, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7],x[8])
        else:
            onea = SingleRun(self._g_path, self._args, x[5], self._args.iterations, x[0],x[1],x[2],x[3],x[4],x[6],x[7])
            
        return (orig, bnxf, onea)

    """
    Fitness function for differential evolution. Fitness is evaluated by comparison of bnx fragment lengths distribution 
    with simulated one feeded by parameters given by differential evolution.

    Returns mean absolute error
    """
    def _fitness(self, x):
        x0 = x[0]
        y0 = x[1]
        x1 = x[2]
        y1 = x[3]
        scale = x[4]
        fpr = x[5]      
        res = x[6]
        miss_k = x[7]

        if self._args.de_obpp:
            res_ebpp = x[8]
            onea = SingleRun(self._g_path, self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,res,miss_k,res_ebpp)
        else:
            onea = SingleRun(self._g_path, self._args, fpr, self._args.iterations, x0,y0,x1,y1,scale,res,miss_k,500)

        #triangle
        mae_r = mae(onea.p0ti, self.bnxp, self._limit)

        if self._args.verbose:
            print("%s\t %.0f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f" % (self._bnx_name, x[4],x[5],x[6],x[7], mae_r, self._min_mae))
        
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

        f.write("Parameter of Bezier curve X0=[x0,y0] and X1=[x1,y1]\n")
        f.write("X0=[%.8f,%.8f]\n" % (self._opt_x[0], self._opt_x[1]))
        f.write("X1=[%.8f,%.8f]\n\n" % (self._opt_x[2], self._opt_x[3]))

        f.write("Fragment lengths described by Bezier curve are smaller than: %.0f\n" % self._opt_x[4])
        f.write("Insertion rate: %.4f\n" % self._opt_x[5])
        f.write("EBPP = %.0f\n" % self._opt_x[6])
        f.write("Digestion probability for fragment lengths greater than Bezier limit: %.4f\n" % self._opt_x[7])
        f.write("OBPP = %.0f\n\n" % self._opt_x[8])

        f.write("Mean absolute error: %.3f" % err)
        
        f.close()


    """
    Stores fragment lengths distribution.
    """
    def _storeDistribution(self, results):
        f = open(self._folder + self._bnx_name + ".csv", "w")

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

        
        

