# Gradient descent to steepest among K random neighbors.
# Much better than looking at all neighbors but still meh
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

    print '\n -------------------- Gradient descent! ---------------------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN(x_range, y_range, z_range)
    prob = throatsN/maxThroatsN

    solution = Solution(target, x_range, y_range, z_range, randomize=True, prob=prob)
    print 'cost = ' + str(solution.cost)

    costs = []

    then = datetime.datetime.now()

    points = solution.throats.keys()
    all_idx = range(len(points))
    random.shuffle(all_idx)

    it = 0
    while it < 1000000000:
        neighbor = None
        idx = np.random.choice(all_idx, 10, replace=False) # TODO param
        for i in idx:
            point = points[i]
            bits = solution.throats[point]
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 1 - bits[b])
                if not neighbor or neighbor[0] > new_cost:
                    neighbor = (new_cost, point, b, 1 - bits[b])
                    break

        cost = neighbor[0]
        point = neighbor[1]
        b = neighbor[2]
        value = neighbor[3]
        #sanity = copy.deepcopy(solution) # uncomment for sanity
        #sanity.set(point, b, value) # uncomment for sanity
        #sanity.recalc() # uncomment for sanity
        solution.setAndRecalc(point, b, value)
        assert solution.cost == cost
        #assert solution.isEqual(sanity) # uncomment for sanity

        print name, 'iter = ', it, ' cost = ', solution.cost, ' time per iter = ', (datetime.datetime.now() - then).total_seconds() / (it + 1)
        print name, 'target = ', target
        print name, 'solution = ', solution.hist

        costs.append(solution.cost)
        if len(costs) > 40 and cost in costs[-40:-20]:
            break

        it += 1

    print 'Total time = ', (datetime.datetime.now() - then).total_seconds()
    solution.sanity()

    return solution
