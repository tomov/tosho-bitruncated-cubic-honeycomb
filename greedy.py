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
    print 'fit = ' + str(solution.fit)

    for it in range(1):
        neighbors = []
        for point, bits in solution.throats.iteritems():
            print 'for ', point
            for b in range(14):
                adj = getAdj(point, b)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                new_fit = solution.fitIfSet(point, b, 1-bits[b], target)

                neighbor = copy.deepcopy(solution)
                neighbor.set(point, b, 1 - bits[b])
                neighbor.recalc(target)
              
                print '    neighbor ', b, ' = ', neighbor.fit, ' vs. ', new_fit
                assert neighbor.fit == new_fit

                #neighbors.append(neighbor)

