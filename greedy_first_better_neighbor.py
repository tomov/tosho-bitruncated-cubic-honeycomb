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

    solution = Solution(target, x_range, y_range, z_range, randomize=True, prob=prob)
    print name, ': initial cost = ' + str(solution.cost)

    then = datetime.datetime.now()
    points = solution.throats.keys()
    cur_idx = 0
    it = 0

    while True:
        neighbor = None
        start_idx = cur_idx
        while True:
            point = points[cur_idx]
            bits = solution.throats[point]
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 1 - bits[b])
                if new_cost < solution.cost:
                    neighbor = (new_cost, point, b, 1 - bits[b])
                    break
            if neighbor: # we found a better neighbor
                break
            cur_idx = (cur_idx + 1) % len(points)
            if cur_idx == start_idx: # we went all the way round
                break

        if not neighbor:
            break # we're at a local minimum
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

        if it % 10000 == 0:
            print '\n', name, 'iter = ', it, ' cost = ', solution.cost, ' secs per iter = %.10lf' % ((datetime.datetime.now() - then).total_seconds() / (it + 1))
            print name, 'target = ', target
            print name, 'solution = ', solution.hist


    print '\n', name, 'FINAL: iter = ', it, ' cost = ', solution.cost, ' secs per iter = %.10lf' % ((datetime.datetime.now() - then).total_seconds() / (it + 1))
    print name, ': target = ', target
    print name, ': solution = ', solution.hist
    print name, ': Total secs = ', (datetime.datetime.now() - then).total_seconds()
    solution.sanity()

    return solution
