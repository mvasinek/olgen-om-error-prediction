import math
import os
import sys
import time

from multiprocessing import Queue, Process, freeze_support

from solution import Solution

"""
Class responsible for simulation start over list of input files.
"""
class Analysis:
    """
    args - collection of input arguments by argparse
    """
    def __init__(self, args):
        self._args = args
        self._solutions = self.__collectSolutions()
    
    """
    returns a list of solution objects
    """
    def __collectSolutions(self):
        solutions = []
        for bnx in self._args.files:
            solutions.append(Solution(self._args.genome, bnx, self._args))
        return solutions

    """
    Starts simulation for each input file as a separate process.
    """
    def store(self):
        procs = []
        for solution in self._solutions:
            p = Process(target=solution.store)
            procs.append((p,solution))

        for p,s in procs:
            p.start()
            if self._args.verbose:
                print("Solution %s started." % s.bnx())
            
        for p,s in procs:
            p.join()
            if self._args.verbose:
                print("Solution %s finished." % s.bnx())

        print("Analysis finished!")
