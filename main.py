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
    parser.add_argument('-p', '--procs', type=int, default=1, help='maximum number of processes to be used (default: 1) to process one file')
    parser.add_argument('-o', '--output', default="results.txt", help='file to store results into')

    parser.add_argument('--store_dist', action='store_true', default=0, help='should the program store the final distribution, the distribution is stored within the directory of input bnx')
    
    parser.add_argument('-g', '--genome', help='genome of interest in CMAP', required=True)
    parser.add_argument('-ag', '--autonoised_genome', help='genome of interest in CMAP for use in autonoise procedure', default='none')
    parser.add_argument('-i', '--iterations', type=int, default=10, help='number of iterations to obtain average distribution for one run of differential evolution')
    
    parser.add_argument('--min_x0', type=float, default=0, help='DE: minimum x0 Bezier parameter')
    parser.add_argument('--min_x1', type=float, default=0, help='DE: minimum x1 Bezier parameter')
    parser.add_argument('--min_y0', type=float, default=0, help='DE: minimum y0 Bezier parameter')
    parser.add_argument('--min_y1', type=float, default=0, help='DE: minimum y1 Bezier parameter')

    parser.add_argument('--max_x0', type=float, default=1, help='DE: maximum x0 Bezier parameter')
    parser.add_argument('--max_x1', type=float, default=1, help='DE: maximum x1 Bezier parameter')
    parser.add_argument('--max_y0', type=float, default=1, help='DE: maximum y0 Bezier parameter')
    parser.add_argument('--max_y1', type=float, default=1, help='DE: maximum y1 Bezier parameter')

    parser.add_argument('--min_bez_lim', type=int, default=0, help='DE: minimum fragment length described by Bezier')
    parser.add_argument('--max_bez_lim', type=int, default=5000, help='DE: maximum fragment length described by Bezier')
    
    parser.add_argument('--min_irate', type=float, default=0, help='DE: minimum insertion rate')
    parser.add_argument('--max_irate', type=float, default=0.3, help='DE: maximum insertion rate')

    parser.add_argument('--min_ebpp', type=int, default=350, help='DE: minimum expected bpp')
    parser.add_argument('--max_ebpp', type=int, default=450, help='DE: maximum expected bpp')

    parser.add_argument('--min_ud', type=float, default=0.8, help='DE: minimum digestion probability greater than Bezier fragment length')
    parser.add_argument('--max_ud', type=float, default=1.0, help='DE: maximum digestion probability greater than Bezier fragment length')

    parser.add_argument('--min_sigma', type=float, default=0, help='DE: minimum standard deviation used in determination of sizing error')
    parser.add_argument('--max_sigma', type=float, default=500.0, help='DE: maximum standard deviation used in determination of sizing error')

    parser.add_argument('--min_det_sigma', type=float, default=0, help='DE: minimum standard deviation used in light intensity spread over pixels')
    parser.add_argument('--max_det_sigma', type=float, default=400.0, help='DE: maximum standard deviation used in light intensity spread over pixels')

    parser.add_argument('--min_det_thr', type=float, default=0.1, help='DE: minimum light intensity threshold to consider pixel a site')
    parser.add_argument('--max_det_thr', type=float, default=0.5, help='DE: maximum light intensity threshold to consider pixel a site')

    parser.add_argument('--use_laplace', action='store_true', default=0, help='Sizing error computed with Laplace distribution')
    parser.add_argument('--use_laplace_li', action='store_true', default=0, help='Sizing error computed with Laplace distribution based on Li et al.')
    parser.add_argument('--laplace_levels', type=int, default=5, help='Number of levels of Laplace distribution')
    parser.add_argument('--laplace_step', type=float, default=500.0, help='Step determines a size of the level, for instance step 1000 and levels 3 yields intervals: <0;1000>;(1000;2000>;(2000;infty)')

    parser.add_argument('--de_obpp', action='store_true', default=0, help='Use perturbation of OBPP in differential evolution, otherwise set to default=500')
    parser.add_argument('--min_obpp', type=int, default=470, help='DE: minimum operating bpp')
    parser.add_argument('--max_obpp', type=int, default=550, help='DE: maximum operating bpp')

    parser.add_argument('--dist_max_len', type=int, default=25000, help='Maximum fragment length to be considered for distribution comparison and output')

    parser.add_argument('--de_max_iter', type=int, default=10, help='Differential evolution: the maximum number of generations over which the entire population is evolved')
    parser.add_argument('--de_pop_size', type=int, default=15, help='A multiplier for setting the total population size. See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.differential_evolution.html')

    parser.add_argument('--disable_camera', action='store_true', default=0, help='Disables camera transformation.')
    parser.add_argument('--fragment_correction', action='store_true', default=0, help='Only fragment length correction due to separation into molecules applied.')
    parser.add_argument('--insertion_correction', action='store_true', default=0, help='Only insertion and fragment length correction.')
    parser.add_argument('--indels_correction', action='store_true', default=0, help='Only missing, insertion and fragment length correction.')
    parser.add_argument('--bezier_correction', action='store_true', default=0, help='Missing with bezier coefficients, insertion and fragment length correction.')
    parser.add_argument('--indels_sizing_correction', action='store_true', default=0, help='Missing, insertion, sizing error and fragment length correction.')

    
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
    
    
    
