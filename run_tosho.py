# Generates solutions for a given N and target distribution
# Usage: python run_tosho.py [input filename] [output directory] [number of threads]
# Example usage: python run_tosho.py input.txt tosho_solutions 4
#

from __future__ import division
import sys
import string

from solution import *
from greedy_first_better_neighbor import *
#from archive.greedy_best_neighbor import *
from merge_outs import *

from multiprocessing import Process, Value

x_range = None
y_range = None
z_range = None
poresN = None
solutionsN = None
target = None
solutions_count = None
output_dir = None

def getNtosho(poresN):
    return poresN ** 3 + (poresN - 1) ** 3

def genHash(length=10):
    return ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(length))

def infinite_greedy(thd):
    while solutions_count.value < solutionsN:
        solution = greedy(target, x_range, y_range, z_range, "Thread #%d" % thd)
        prefix = "%s/" % output_dir
        suffix = "_N_%d_cost_%d_%s" % (poresN, solution.cost, genHash())
        coords_file, neigh_file, neigh_WTF_file = solution.exportForTosho(prefix, suffix)
        solutions_count.value += 1
        merge_outs(coords_file, neigh_file, 'conductance/sheeeiiiiit.csv')
        merge_outs(coords_file, neigh_WTF_file, 'conductance/fuuuckk.csv')

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    threadsN = int(sys.argv[3])
    assert threadsN > 0 and threadsN < 50

    print "Input filename = ", input_file
    print "Output dir = ", output_dir
    print "# of threads = ", threadsN

    with open(input_file, "r") as f:
        # Read # of solutions to generate
        solutionsN = int(f.readline())
        print "Input # solutions = ", solutionsN

        # Read # of pores
        #
        poresN = int(f.readline())
        print "Input # pores = ", poresN
        x_range = [0, poresN * 2 - 2]
        y_range = x_range
        z_range = x_range
        assert getN(x_range, y_range, z_range) == getNtosho(poresN)

        # Read distribution
        # Assumes it's whitespace-separated
        #
        distribution = [float(x) for x in f.readline().split()]
        print "Input raw distribution = ", distribution
        target = truncateAndScale(distribution, x_range, y_range, z_range)
        print "Normalized = ", target

    solutions_count = Value('i', 0)
    procs = [Process(target=infinite_greedy, args=(i,)) for i in range(threadsN)]
    [p.start() for p in procs]
    [p.join() for p in procs]
