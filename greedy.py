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


def bestInRange(solution, neighbors, points, chunks, i, x_range, y_range, z_range):
    for idx in chunks[i]:
        point = points[idx]
        bits = solution.throats[point]
        for b in range(14):
            adj = getAdj(point, b)
            if not isIn(adj, x_range, y_range, z_range):
                continue
            new_cost = solution.costIfSet(point, b, 1 - bits[b])
            if not neighbors[i] or new_cost < neighbors[i][0]:
                neighbors[i] = (new_cost, point, b, 1 - bits[b])

def greedy(target, x_range, y_range, z_range):

    print '\n -------------------- Gradient descent! ---------------------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN(x_range, y_range, z_range)
    prob = throatsN/maxThroatsN

    solution = Solution(target, x_range, y_range, z_range, randomize=True, prob=prob)
    print 'cost = ' + str(solution.cost)

    costs = []

    threadsN = 8

    manager = Manager()
    idx = range(len(solution.throats))
    chunks = manager.list()
    for chunk in chunkate(idx, threadsN):
        chunks.append(chunk)
    points = solution.throats.keys()

    for it in range(100000):
        neighbors = manager.list([None] * threadsN)
        procs = [Process(target=bestInRange, args=(solution, neighbors, points, chunks, i, x_range, y_range, z_range)) for i in range(threadsN)]
        [p.start() for p in procs]
        [p.join() for p in procs]

        print 'neighbors = ', neighbors
        best = None
        for neighbor in neighbors:
            if not best or neighbor[0] < best[0]:
                best = neighbor

        neighbor = None
        for point, bits in solution.throats.iteritems():
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 1 - bits[b])
                if not neighbor or new_cost < neighbor[0]:
                    neighbor = (new_cost, point, b, 1 - bits[b])

        print 'best =', best
        print 'neighbor = ', neighbor
        assert neighbor == best

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
