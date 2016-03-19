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

    cur_idx = 0
    for it in range(100000):
        neighbor = None
        #idx = np.random.choice(all_idx, 10, replace=False) # TODO param
        start_idx = cur_idx
        terminate = False
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
            if neighbor:
                break
            cur_idx = (cur_idx + 1) % len(points)
            if cur_idx == start_idx:
                terminate = True # we've reached a local minimum
                break

        if terminate:
            break
        best = neighbor

        cost = best[0]
        point = best[1]
        b = best[2]
        value = best[3]
        #sanity = copy.deepcopy(solution)
        #sanity.set(point, b, value)
        #sanity.recalc()
        solution.setAndRecalc(point, b, value)
        assert solution.cost == cost
        #assert solution.isEqual(sanity)

        print name, 'iter = ', it, ' cost = ', solution.cost, ' time per iter = ', (datetime.datetime.now() - then).total_seconds() / (it + 1)
        print name, 'target = ', target
        print name, 'solution = ', solution.hist

        costs.append(solution.cost)
        if len(costs) > 40 and cost in costs[-40:-20]:
            break

    print 'Total time = ', (datetime.datetime.now() - then).total_seconds()

    return solution
