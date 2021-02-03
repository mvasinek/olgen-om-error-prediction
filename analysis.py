import math
import os
import sys
import time

from multiprocessing import Queue, Process, freeze_support

from solution import Solution, LaplaceSolution, NoCameraSolution, OnlyCorrectionSolution, InsertionCorrectionSolution, IndelsCorrectionSolution, BezierCorrectionSolution, IndelsSizingCorrectionSolution

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
            if self._args.use_laplace:
                solutions.append(LaplaceSolution(self._args.genome, bnx, self._args))
            elif self._args.disable_camera:
                solutions.append(NoCameraSolution(self._args.genome, bnx, self._args))
            elif self._args.fragment_correction:
                solutions.append(OnlyCorrectionSolution(self._args.genome, bnx, self._args))
            elif self._args.insertion_correction:
                solutions.append(InsertionCorrectionSolution(self._args.genome, bnx, self._args))
            elif self._args.indels_correction:
                solutions.append(IndelsCorrectionSolution(self._args.genome, bnx, self._args))
            elif self._args.bezier_correction:
                solutions.append(BezierCorrectionSolution(self._args.genome, bnx, self._args))
            elif self._args.indels_sizing_correction:
                solutions.append(IndelsSizingCorrectionSolution(self._args.genome, bnx, self._args))
            else:
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
