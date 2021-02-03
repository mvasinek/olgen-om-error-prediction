import random

from molecule_generator import Molecule

from multiprocessing import Queue

class FalseSites:
    def __init__(self, false_rate, genome):
        self._false_rate = false_rate
        self._genome = genome

    def processOne(self, molecules, queue=None):
        ir = self._genome.insertionRate(self._false_rate)

        for molecule in molecules:
            m_len = molecule.length
            
            #expected false positive
            E_fp = ir*m_len
            
            res = E_fp - int(E_fp)
            rd = random.random()

            if rd < res or E_fp == 0:
                E_fp = int(E_fp)            
            else:
                E_fp = int(E_fp) + 1

            for _ in range(E_fp):
                rd_pos = random.randint(0, m_len)
                molecule.pos.append(rd_pos)
            
            molecule.pos.sort()

        if queue != None:
            queue.put(molecules)