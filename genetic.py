from __future__ import division
import random
import operator
import copy
import datetime

from solution import *


def genetic(target, x_range, y_range, z_range):

    print '\n -------------------- Genetic! ---------------------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN(x_range, y_range, z_range)
    throatIdxs = getRandomShuffleThroatIdxs(x_range, y_range, z_range) # randomly assign indexes to the throats for the genetic algorithm

    print throatsN, ' / ', maxThroatsN, ' = ', throatsN / maxThroatsN

    populationN = 50 # TODO param
    gensN = 30 # TODO param
    prob = throatsN/maxThroatsN # TODO try with prob = 0.5
    mutationProb = 1 / 10000 # TODO param
    crossoverProb = 0.7 # TODO param
    population = [Solution(target, x_range, y_range, z_range, randomize=True, prob=prob) for _ in range(populationN)]

    fitAvgs = []
    fitBests = []

    best = copy.deepcopy(population[0])

    for gen in range(gensN):
        fitTotal = sum([solution.fit for solution in population])
        fitAvg = fitTotal / populationN
        fitAvgs.append(fitAvg)
        bestIdx, fitBest = max(enumerate([solution.fit for solution in population]), key=operator.itemgetter(1))
        fitBests.append(fitBest)
        #new_population = [population[bestIdx]] # always pass on the best TODO param
        new_population = []
        if best.fit < population[bestIdx].fit:
            best = copy.deepcopy(population[bestIdx])

        print datetime.datetime.now(), ': At gen ', gen, ': fitAvg = ', fitAvg, ' fitBest = ', fitBest, ' costBest = ', 1 / fitBest
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
            mom.recalc()
            dad.recalc()

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
            solution.recalc()

        population = new_population

    print 'BEST :  (target = ', target, ')'
    best.output()
