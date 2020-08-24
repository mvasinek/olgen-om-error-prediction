import argparse
import os

"""
Copyright
Michal Vasinek
michal.vasinek@vsb.cz
19.8.2020 Ostrava, Czech Republic
"""

from analysis import Analysis

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Program for optical mapping error prediction.", prog="OM Error Predictor")

    parser.add_argument('files', nargs='*', help='list of input bnx files')
    #parser.add_argument('-p', '--procs', type=int, default=1, help='maximum number of processes to be used (default: 1)')
    parser.add_argument('-o', '--output', default="results.txt", help='file to store results into')

    parser.add_argument('--store_dist', action='store_true', default=0, help='should the program store the final distribution, the distribution is stored within the directory of input bnx')
    
    parser.add_argument('-g', '--genome', help='genome of interest in CMAP', required=True)
    parser.add_argument('-i', '--iterations', type=int, default=10, help='number of iterations to obtain average distribution for one run of differential evolution')
    
    parser.add_argument('--min_x0', type=float, default=0, help='DE: minimum x0 Bezier parameter')
    parser.add_argument('--min_x1', type=float, default=0, help='DE: minimum x1 Bezier parameter')
    parser.add_argument('--min_y0', type=float, default=0, help='DE: minimum y0 Bezier parameter')
    parser.add_argument('--min_y1', type=float, default=0, help='DE: minimum y1 Bezier parameter')

    parser.add_argument('--max_x0', type=float, default=1, help='DE: maximum x0 Bezier parameter')
    parser.add_argument('--max_x1', type=float, default=1, help='DE: maximum x1 Bezier parameter')
    parser.add_argument('--max_y0', type=float, default=1, help='DE: maximum y0 Bezier parameter')
    parser.add_argument('--max_y1', type=float, default=1, help='DE: maximum y1 Bezier parameter')

    parser.add_argument('--min_bez_lim', type=int, default=1000, help='DE: minimum fragment length described by Bezier')
    parser.add_argument('--max_bez_lim', type=int, default=2200, help='DE: maximum fragment length described by Bezier')
    
    parser.add_argument('--min_irate', type=float, default=0, help='DE: minimum insertion rate')
    parser.add_argument('--max_irate', type=float, default=0.3, help='DE: maximum insertion rate')

    parser.add_argument('--min_ebpp', type=int, default=400, help='DE: minimum expected bpp')
    parser.add_argument('--max_ebpp', type=int, default=520, help='DE: maximum expected bpp')

    parser.add_argument('--min_ud', type=float, default=0.9, help='DE: minimum digestion probability greater than Bezier fragment length')
    parser.add_argument('--max_ud', type=float, default=1.0, help='DE: maximum digestion probability greater than Bezier fragment length')

    parser.add_argument('--de_obpp', action='store_true', default=0, help='use perturbation of OBPP in differential evolution, otherwise set to default=500')
    parser.add_argument('--min_obpp', type=int, default=470, help='DE: minimum operating bpp')
    parser.add_argument('--max_obpp', type=int, default=550, help='DE: maximum operating bpp')

    parser.add_argument('--dist_max_len', type=int, default=25000, help='Maximum fragment length to be considered for distribution comparison and output')

    parser.add_argument('--de_max_iter', type=int, default=10, help='Differenetial evolution: the maximum number of generations over which the entire population is evolved')
    parser.add_argument('--de_pop_size', type=int, default=15, help='A multiplier for setting the total population size. See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html')
    
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')    

    args = parser.parse_args()

    for fpath in args.files:
        if not os.path.exists(fpath):
            print("Input file: %s does not exists." % fpath)
            exit(0)

    if not os.path.exists(args.genome):
        print("Genome: %s does not exists." % args.genome)
        exit(0)
        
    
    a = Analysis(args)
    a.store()

    exit(0)
    
    
    
