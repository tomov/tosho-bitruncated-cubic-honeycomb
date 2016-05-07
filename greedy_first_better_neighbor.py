# Gradient descent w/ heuristic --
# At each step, go to first neighbor that is better
# Keep iterating over neighbors in round-robin fashion across steps
#

from __future__ import division
import random
import operator
import copy
import datetime
import numpy as np
import math
import datetime

from solution import *

def greedy(target, x_range, y_range, z_range, name):

    print '\n', name, ': ----------- New solution ---------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN(x_range, y_range, z_range)
    prob = throatsN/maxThroatsN

    solution = Solution(target, x_range, y_range, z_range, randomize=False, prob=prob)
    print name, ': initial cost = ' + str(solution.cost)

    # randomize the pores
    #
    points = solution.throats.keys()
    all_idx = range(len(points))
    random.shuffle(all_idx)

    # randomize the neighbors of each pore i.e. the throats
    #
    adj_idxs = []
    for i in range(len(all_idx)):
        r = range(14)
        random.shuffle(r)
        adj_idxs.append(r)

    then = datetime.datetime.now()

    cur_idx = 0
    it = 0
    # while we keep finding better solutions
    #
    while True:
        neighbor = None
        start_idx = cur_idx
        # keep going round-robin across all pores
        #
        while True:
            idx = all_idx[cur_idx]
            point = points[idx]
            bits = solution.throats[point]
            # for each pore, try toggling all throats coming out of it
            #
            for b in adj_idxs[idx]:
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 1 - bits[b])
                # see if the "neighboring" solution is better
                #
                if new_cost < solution.cost:
                    neighbor = (new_cost, point, b, 1 - bits[b])
                    break
            # stop if we found a better "neighboring" solution
            #
            if neighbor:
                break
            cur_idx = (cur_idx + 1) % len(points)
            # if we went all the way round and found nothing => call it quits
            #
            if cur_idx == start_idx:
                break

        # we went all the way and didn't find a better neighbor => we're at a local minimum
        #
        if not neighbor:
            break
        cost = neighbor[0]
        point = neighbor[1]
        b = neighbor[2]
        value = neighbor[3]

        #sanity = copy.deepcopy(solution) # uncomment for sanity

        solution.setAndRecalc(point, b, value)
        assert solution.cost == cost

        #sanity.set(point, b, value) # uncomment for sanity
        #sanity.recalc() # uncomment for sanity
        #assert solution.isEqual(sanity) # uncomment for sanity

        it += 1

        if it % 1000 == 0:
            print '\n', name, 'iter = ', it, ' cost = ', solution.cost, ' secs per iter = %.10lf' % ((datetime.datetime.now() - then).total_seconds() / (it + 1))
            print name, 'target = ', target
            print name, 'solution = ', solution.hist


    print '\n', name, 'FINAL: iter = ', it, ' cost = ', solution.cost, ' secs per iter = %.10lf' % ((datetime.datetime.now() - then).total_seconds() / (it + 1))
    print name, ': target = ', target
    print name, ': solution = ', solution.hist
    print name, ': Total secs = ', (datetime.datetime.now() - then).total_seconds()
    solution.sanity()

    return solution
