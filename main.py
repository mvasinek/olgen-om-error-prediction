import argparse
import os
from multiprocessing import freeze_support

"""
Copyright
Michal Vasinek
michal.vasinek@vsb.cz
31.1.2021 Ostrava, Czech Republic
"""

from analysis import Analysis

if __name__ == "__main__":
    freeze_support()

    parser = argparse.ArgumentParser(description="Program for optical mapping error prediction.", prog="OM Error Predictor")

    parser.add_argument('files', nargs='*', help='list of input bnx files')
    parser.add_argument('-p', '--procs', type=int, default=1, help='maximum number of processes to be used (default: 1) to process one file')
    parser.add_argument('-o', '--output', default="results.txt", help='file to store results into')

    parser.add_argument('--store_dist', action='store_true', default=0, help='should the program store the final distribution, the distribution is stored within the directory of input bnx')
    parser.add_argument('--generate_bnx', action='store_true', default=0, help='should the program eventually generate bnx')
    parser.add_argument('--single_sim', action='store_true', default=0, help='should the program run only single simulation with prefeined parameters')
    
    parser.add_argument('-g', '--genome', help='genome of interest in CMAP', required=True)
    parser.add_argument('-i', '--iterations', type=int, default=10, help='number of iterations to obtain average distribution for one run of differential evolution')
        
    parser.add_argument('--min_irate', type=float, default=0, help='DE: minimum insertion rate')
    parser.add_argument('--max_irate', type=float, default=0.3, help='DE: maximum insertion rate')

    parser.add_argument('--min_mrate', type=float, default=0, help='DE: minimum miss rate')
    parser.add_argument('--max_mrate', type=float, default=0.5, help='DE: maximum miss rate')
    
    parser.add_argument('--min_sigma', type=float, default=0, help='DE: minimum standard deviation used in determination of sizing error')
    parser.add_argument('--max_sigma', type=float, default=1.0, help='DE: maximum standard deviation used in determination of sizing error')

    parser.add_argument('--min_light_sigma', type=float, default=100, help='DE: minimum standard deviation used in light intensity spread over pixels')
    parser.add_argument('--max_light_sigma', type=float, default=700.0, help='DE: maximum standard deviation used in light intensity spread over pixels')

    parser.add_argument('--min_stretch', type=float, default=0.5, help='DE: minimum stretch factor')
    parser.add_argument('--max_stretch', type=float, default=1.5, help='DE: maximum stretch factor')
    
    parser.add_argument('--dist_max_len', type=int, default=25000, help='Maximum fragment length to be considered for distribution comparison and output')
    parser.add_argument('--sites_max_len', type=int, default=100, help='Maximum number of sites per molecule')

    parser.add_argument('--de_max_iter', type=int, default=10, help='Differential evolution: the maximum number of generations over which the entire population is evolved')
    parser.add_argument('--de_pop_size', type=int, default=15, help='A multiplier for setting the total population size. See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html')

    parser.add_argument('--missing_iterations', type=int, default=100, help='The number of iterations to be used for generation of digestion probability, higher number means better precision. To regenerate digestion probability delete pd100 directory.')

    parser.add_argument('-t', '--eval_type', default="distance", help='Evaluation type: distance, sites or both')
    
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')    

    args = parser.parse_args()

    for fpath in args.files:
        if not os.path.exists(fpath):
            print("Input file: %s does not exists." % fpath)
            exit(1)

    if not os.path.exists(args.genome):
        print("Genome: %s does not exists." % args.genome)
        exit(1)
        
    a = Analysis(args)    
    a.store()

    exit(0)
    
    
    
