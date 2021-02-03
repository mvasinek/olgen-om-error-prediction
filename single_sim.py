from multiprocessing import Queue, Process
import time

from bnx import BNX
from false_sites import FalseSites
from molecule_generator import MoleculeGenerator, Molecule
from optical_variation import OpticalVariation
from optics import Optics

class SingleSim:
    def __init__(self, args, bnx_path, genome, irate, erange, vmiss, stretch_factor, sizing, pulse_sigma):
        self.args = args
        self.bnx_path = bnx_path
        self.genome = genome
        self.iterations = args.iterations        
        self.m_generator = MoleculeGenerator(bnx_path, self.genome, 1)

        self.irate = irate
        self.erange = erange
        self.vmiss = float("%.2f" % vmiss)
        self.sf = stretch_factor
        self.sizing = sizing
        self.p_sigma = pulse_sigma
        
        self.p_dist,self.p_sites = self.processFull()

    """
    Simulation is performed self.iterations times, each by separate process.
    """
    def processFull(self):
        pt = None
        ps = None

        queue = Queue()

        procs = [Process(target=self.processOne, args=(queue,)) for i in range(self.iterations)]

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

        for value_pair in results:
            ti = value_pair[0]
            if pt == None:
                pt = ti                
            else:
                for i in range(self.args.dist_max_len):
                    pt[i] += ti[i]

        for value_pair in results:
            si = value_pair[1]
            if ps == None:
                ps = si
            else:
                for i in range(self.args.sites_max_len):
                    ps[i] += si[i]        

        m = []
        for _, _, mols in results:
            m += mols

        if self.args.generate_bnx:
            BNX.ToBNX(self.bnx_path.replace(".bnx","_generated.bnx"), m)            

        for i in range(self.args.dist_max_len):
            pt[i] = pt[i] / self.iterations

        for i in range(self.args.sites_max_len):
            ps[i] = ps[i] / self.iterations

        return (pt,ps)
    
    def processOne(self, queue):        
        #bionano dataset specific setting
        if "GM09888" in self.bnx_path or "GM11428" in self.bnx_path or "GM24143" in self.bnx_path:
            molecules = self.m_generator.generate(self.erange, self.vmiss, "female")
        else:
            molecules = self.m_generator.generate(self.erange, self.vmiss, "male")

        ins = FalseSites(self.irate, self.genome)
        ins.processOne(molecules)
        
        optics = Optics(molecules, 375, self.sizing, self.p_sigma)
        optics.processOne()

        sf = OpticalVariation(molecules, self.sf)
        p_result = sf.processOne()

        #bionano dataset specific setting
        if "GM24149" in self.bnx_path or "GM24143" in self.bnx_path or "GM24385" in self.bnx_path:
            queue.put((p_result, self.m_generator.sitesDist(molecules, False), molecules))
        else:
            queue.put((p_result, self.m_generator.sitesDist(molecules, True), molecules))