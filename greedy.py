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

    solution = Solution(target, x_range, y_range, z_range)

    buckets = [14, 0, 13, 1, 12, 2, 11, 3, 10, 4, 9, 5, 8, 6, 7]
    points = solution.throats.keys()
    for bucket in buckets:
        print 'at ' + str(bucket)
        which = np.random.choice(solution.throats, points)
        for p in which:
            bs = np.random.choice(solution.throats, :)
            points.remove(p)


    populationN = 50 # TODO param
    gensN = 30 # TODO param
    mutationProb = 1 / 10000 # TODO param
    crossoverProb = 0.7 # TODO param
    population = [Solution(target, x_range, y_range, z_range, randomize=True, prob=prob) for _ in range(populationN)]

    fitTotals = []
    fitBests = []

    best = copy.deepcopy(population[0])

    for gen in range(gensN):
        fitTotal = sum([solution.fit for solution in population])
        fitTotals.append(fitTotal)
        bestIdx, fitBest = max(enumerate([solution.fit for solution in population]), key=operator.itemgetter(1))
        fitBests.append(fitBest)
        #new_population = [population[bestIdx]] # always pass on the best TODO param
        new_population = []
        if best.fit < population[bestIdx].fit:
            best = copy.deepcopy(population[bestIdx])

        print datetime.datetime.now(), ': At gen ', gen, ': fitTotal = ', fitTotal, ' fitBest = ', fitBest
        #print 'fitness = ', sorted([solution.fit for solution in population])

        # roulette wheel
        #
        while len(new_population) < populationN:
            r = random.random()
            cumsum = 0
            for solution in population:
                cumsum += solution.fit / fitTotal
                if cumsum > r:
                    new_population.append(copy.deepcopy(solution))
                    break

        # crossover
        #
        for i in range(len(population) // 2):
            mom = population[i * 2]
            dad = population[i * 2 + 1]
            if random.random() < crossoverProb:
                switchPos = int(random.random() * maxThroatsN)
                for point in mom.throats.keys():
                    for b in range(14):
                        if b < revAdjIdx[b]:
                            adj = getAdj(point, b)
                            if not isIn(adj, x_range, y_range, z_range):
                                continue
                            if throatIdxs[point][b] < switchPos:
                                continue
                            bitMom = mom.throats[point][b]
                            bitDad = dad.throats[point][b]
                            mom.set(point, b, bitDad)
                            dad.set(point, b, bitMom)
            mom.recalc(target)
            dad.recalc(target)

        # mutate
        #
        for solution in population:
            for point, bits in solution.throats.iteritems():
                for b in range(14):
                    if b < revAdjIdx[b]:
                        adj = getAdj(point, b)
                        if not isIn(adj, x_range, y_range, z_range):
                            continue
                        if random.random() < mutationProb:
                            # print 'flip!'
                            solution.set(point, b, 1 - bits[b])
            solution.recalc(target)

        population = new_population

    print 'BEST :  (target = ', target, ')'
    best.output()
