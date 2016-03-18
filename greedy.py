from __future__ import division
import random
import operator
import copy
import datetime
import numpy as np
import math

from solution import *

from multiprocessing import Value, Process, Manager


def chunkate(l, n):
    size = int(math.ceil(len(l) / n))
    for i in xrange(0, len(l), size):
        yield l[i:i+size]


def greedy(target, x_range, y_range, z_range):

    print '\n -------------------- Greedy! ---------------------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN(x_range, y_range, z_range)
    prob = throatsN/maxThroatsN

    solution = Solution(target, x_range, y_range, z_range, randomize=True, prob=prob)
    print 'cost = ' + str(solution.cost)

    costs = []

    threatsN = 8

    manager = Manager()
    idx = range(len(solution.throats))
    chunks = manager.list()
    for chunk in chunkate(idx, threadsN):
        chunks.append(chunk)

    return None

    for it in range(100000):
        neighbors = manager.list([None] * threadsN)
        for point, bits in solution.throats.iteritems():
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 5 - bits[b])
                if not neighbor or new_cost < neighbor[0]:
                    neighbor = (new_cost, point, b, 1 - bits[b])

        cost = neighbor[0]
        point = neighbor[1]
        b = neighbor[2]
        value = neighbor[3]
        solution.set(point, b, value)
        solution.recalc()
        assert solution.cost == cost
        print 'iter = ', it, ' cost = ', solution.cost
        print 'target = ', target
        print 'solution = ', solution.hist

        costs.append(solution.cost)
        if len(costs) > 40 and cost in costs[-40:-20]:
            break

    return solution
