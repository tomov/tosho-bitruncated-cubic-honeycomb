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

    for it in range(100000):
        neighbor = None
        for point, bits in solution.throats.iteritems():
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 1 - bits[b])
                if not neighbor or new_cost < neighbor[0]:
                    neighbor = (new_cost, point, b, 1 - bits[b])
        best = neighbor

        cost = best[0]
        point = best[1]
        b = best[2]
        value = best[3]
        solution.set(point, b, value)
        solution.recalc()
        assert solution.cost == cost
        print name, 'iter = ', it, ' cost = ', solution.cost, ' time per iter = ', (datetime.datetime.now() - then).total_seconds() / (it + 1)
        print name, 'target = ', target
        print name, 'solution = ', solution.hist

        costs.append(solution.cost)
        if len(costs) > 40 and cost in costs[-40:-20]:
            break

    return solution
