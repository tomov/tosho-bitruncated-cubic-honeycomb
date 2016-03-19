from __future__ import division

from solution import *
#from genetic import *
#from greedy_best_of_first_k import *
from greedy_first_better_neighbor import *

from multiprocessing import Process

# assumes the corner points are at even coordinates
# and that center points are at odd coordinates
# see https://en.wikipedia.org/wiki/Cubic_crystal_system

x_range = [0, 200]
y_range = [0, 200]
z_range = [0, 200]

ASA_Sandstone_1 = [0.1590810577, 0.1239705245, 0.1352405722, 0.1404421326, 0.1005635024, 0.07195491981, 0.05938448201, 0.0455136541, 0.04421326398, 0.03077589944, 0.02080624187, 0.01517121803, 0.0125704378, 0.00736887733, 0.00736887733, 0.004768097096, 0.004768097096, 0.002600780234, 0.003034243606, 0.001733853489, 0.002167316862, 0.001300390117, 0.0008669267447, 0.0004334633723, 0.0008669267447, 0.001300390117, 0.0004334633723, 0.0004334633723, 0, 0, 0, 0, 0, 0.0004334633723, 0, 0, 0, 0, 0.0004334633723]


def infinite_greedy(thd):
    best = None
    while True:
        new_solution = greedy(target, x_range, y_range, z_range, str(thd))
        if not best or best.cost > new_solution.cost:
            best = new_solution
            best.serialize("best_greedy_" + str(thd) + ".json")
        break

if __name__ == '__main__':
    print truncate(ASA_Sandstone_1)
    target = truncateAndScale(ASA_Sandstone_1, x_range, y_range, z_range)
    print target

    procs = [Process(target=infinite_greedy, args=(i,)) for i in range(1)]
    [p.start() for p in procs]
    [p.join() for p in procs]

    #genetic(target, x_range, y_range, z_range)
