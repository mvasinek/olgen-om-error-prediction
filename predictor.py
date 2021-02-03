import copy
import os

from bnx import BNX, BNXData
from genome import OMGenome
from single_sim import SingleSim
from missing_site import EnzymeMissing
from molecule_generator import MoleculeGenerator

from scipy.optimize import rosen, differential_evolution

"""
mean absolute error between two lists of values a,b
lim - variable limiting the index in a,b
"""
def mae(a,b,lim):
    mae = 0
    #for i in range(1500, lim):
    for i in range(lim):
        mae += abs(a[i]-b[i])
    return mae

class Predictor:
    def __init__(self, bnx_path, args):
        self._g_path = args.genome
        self._bnx_path = bnx_path
        self._args = args

        #directory to store output
        self._folder, self._bnx_name = self.__parsePath(bnx_path)

        self._limit = args.dist_max_len
        self._min_mae = 1000

        self._opt_x = None
             
        self._orig_genome = OMGenome(self._g_path)
        self._orig_pos = copy.deepcopy(self._orig_genome.pos)

        if not os.path.exists("pd100"):
            print("Missing site table has not been generated yet.")
            print("Generation of missing site table has started. It may take a while to produce. This action is performed only once.")

            EnzymeMissing.GeneratePDS(self._orig_genome, self._args.missing_iterations)
            EnzymeMissing.ReloadPDS()

        self._bnxf = None

        if self._bnx_path.endswith(".data"):
            #assume simplified bnx data
            self._bnxf = BNXData(self._bnx_path)
        else:
            #standard bnx            
            self._bnxf = BNX(self._bnx_path)

        self.bnxp = self._bnxf.p0

        self.bnx_sites_dist = MoleculeGenerator.SitesDistBNX(BNX.ParseSitesNo(self._bnx_path))

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

    def _collect(self):
        orig = self._orig_genome
        
        a = self._args

        v_miss = (a.min_mrate, a.max_mrate)
        ir = (a.min_irate, a.max_irate)
        sf = (a.min_stretch, a.max_stretch)
        pulse = (a.min_light_sigma, a.max_light_sigma)        
        sizing_coef = (a.min_sigma, a.max_sigma)
        
        bounds = (v_miss, ir, sf, sizing_coef, pulse)        

        result = differential_evolution(self._fitness, bounds, maxiter=self._args.de_max_iter, popsize=self._args.de_pop_size)

        x = result.x

        onea = SingleSim(self._args, self._bnx_path, self._orig_genome, x[1], 1000, x[0], x[2], x[3], x[4])

        return (orig, self._bnxf, onea)

    def _fitness(self, x):
        miss_rate = x[0]
        irate = x[1]        
        sf = x[2]
        sizing = x[3]
        det_sigma = x[4]     

        self._orig_genome.pos = copy.deepcopy(self._orig_pos)
        
        onea = SingleSim(self._args, self._bnx_path, self._orig_genome, irate, 1000, miss_rate, sf, sizing, det_sigma)        

        #triangle
        if self._args.eval_type == "distance":
            mae_r = mae(onea.p_dist, self.bnxp, self._limit)
        elif self._args.eval_type == "sites":
            mae_r = mae(onea.p_sites, self.bnx_sites_dist, self._args.sites_max_len)
        elif self._args.eval_type == "both":
            mae_r = mae(onea.p_dist, self.bnxp, self._limit) + mae(onea.p_sites, self.bnx_sites_dist, self._args.sites_max_len)

        if self._args.verbose:            
            print("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f" % (self._bnx_path, x[0],x[1],x[2],x[3],x[4], mae_r, self._min_mae))            

            f = open(self._folder + self._args.output + ".out", "a")            
            f.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (self._bnx_path, x[0],x[1],x[2],x[3],x[4], mae_r, self._min_mae))            
            f.close()
        
        if mae_r < self._min_mae:
            self._min_mae = mae_r
            self._opt_x = x
        
        return mae_r
    
    """
    Public method to call storage of distribution(optional) and parameters(mandatory)
    """    
    def store(self):
        results = self._collect()
        
        if self._args.store_dist:
            self._storeLengthDistribution(results,1)
            self._storeLengthDistribution(results,20)
            self._storeLengthDistribution(results,100)
            self._storeSitesDistribution(results)
        
        #compute mean absolute error
        err = mae(results[2].p_dist, self.bnxp, self._limit)
        self._storeParameters(err)   

    """
    Stores parameters in a separate file.
    """
    def _storeParameters(self, err):
        f = open(self._folder + self._args.output, "w")

        f.write("Input file: %s\n" % self._bnx_path)
        f.write("Genome file: %s\n\n" % self._g_path)
        
        f.write("Missing rate: %.3f\n" % self._opt_x[0])
        f.write("Insertion rate: %.3f\n" % self._opt_x[1])
        f.write("Stretch factor: %.3f\n" % self._opt_x[2])
        f.write("Sizing error: %.3f\n" % self._opt_x[3])
        f.write("Light intensity spread sigma = %.4f\n\n" % self._opt_x[4])        
        
        f.write("Mean absolute error: %.8f" % err)
        
        f.close()     

    """
    Stores fragment lengths distribution.
    """
    def _storeLengthDistribution(self, results, bin_size=100):
        f = open(self._folder + self._args.output + "_" + str(bin_size) + "_dist.csv", "w")

        f.write("Fragment Length\tGenome\tBNX\tPrediction\n")

        orig = results[0].p0
        bnx = results[1].p0
        tipo = results[2].p_dist

        bin_size = bin_size
        bins_no = self._limit // bin_size

        for i in range(bins_no):
            s_o = 0
            s_b = 0
            s_t = 0
            for j in range(bin_size):
                idx = i*bin_size + j
                s_o += orig[idx]
                s_b += bnx[idx]
                s_t += tipo[idx]

            f.write("%d\t%.8f\t%.8f\t%.8f\n" % ((i+1)*bin_size, s_o/bin_size, s_b/bin_size, s_t/bin_size))
        
        f.close()

    """
    Stores fragment lengths distribution.
    """
    def _storeSitesDistribution(self, results):
        f = open(self._folder + self._args.output + "_sites.csv", "w")

        f.write("Sites Number\tBNX\tPrediction\n")

        sites_dist = results[2].p_sites

        for i in range(len(sites_dist)):
            f.write("%d\t%.8f\t%.8f\n" % (i,self.bnx_sites_dist[i],sites_dist[i]))

        f.close()

class SingleSimPredictor(Predictor):
    def _collect(self):
        orig = self._orig_genome        
        a = self._args
        onea = SingleSim(self._args, self._bnx_path, self._orig_genome, a.min_irate, 1000, a.min_mrate, a.min_stretch, a.min_sigma, a.min_light_sigma)
        self._opt_x = [a.min_mrate, a.min_irate, a.min_stretch, a.min_sigma, a.min_light_sigma]
        return (orig, self._bnxf, onea)