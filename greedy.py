from __future__ import division
import random
import operator
import copy
import datetime
import numpy as np

from solution import *


def greedy(target, x_range, y_range, z_range):

    print '\n -------------------- Greedy! ---------------------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN(x_range, y_range, z_range)
    throatIdxs = getRandomShuffleThroatIdxs(x_range, y_range, z_range) # randomly assign indexes to the throats
    prob = throatsN/maxThroatsN

    solution = Solution(target, x_range, y_range, z_range, randomize=True, prob=prob)
    print 'cost = ' + str(solution.cost)

    for it in range(1000):
        neighbors = []
        for point, bits in solution.throats.iteritems():
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_cost = solution.costIfSet(point, b, 1 - bits[b], target)
                neighbors.append((new_cost, point, b, 1 - bits[b]))

        neighbors.sort()
        cost = neighbors[0][0]
        point = neighbors[0][1]
        b = neighbors[0][2]
        value = neighbors[0][3]
        solution.set(point, b, value)
        solution.recalc(target)
        assert solution.cost == cost
        print 'cost = ', solution.cost, ' cost = ', solution.cost
        print 'target = ', target
        print 'solution = ', solution.hist

