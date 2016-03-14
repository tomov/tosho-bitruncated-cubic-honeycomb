from __future__ import division
import random
import operator
import copy
import datetime

# assumes the corner points are at even coordinates
# and that center points are at odd coordinates
# see https://en.wikipedia.org/wiki/Cubic_crystal_system

x_range = [0, 10]
y_range = [0, 10]
z_range = [0, 10]

ASA_Sandstone_1 = [0.1590810577, 0.1239705245, 0.1352405722, 0.1404421326, 0.1005635024, 0.07195491981, 0.05938448201, 0.0455136541, 0.04421326398, 0.03077589944, 0.02080624187, 0.01517121803, 0.0125704378, 0.00736887733, 0.00736887733, 0.004768097096, 0.004768097096, 0.002600780234, 0.003034243606, 0.001733853489, 0.002167316862, 0.001300390117, 0.0008669267447, 0.0004334633723, 0.0008669267447, 0.001300390117, 0.0004334633723, 0.0004334633723, 0, 0, 0, 0, 0, 0.0004334633723, 0, 0, 0, 0, 0.0004334633723]

# Total number of lattice points in range
#
def getN(x_range, y_range, z_range):
    xs = x_range[1] - x_range[0]
    ys = y_range[1] - y_range[0]
    zs = z_range[1] - z_range[0]
    ret = (xs // 2 + 1) * (ys // 2 + 1) * (zs // 2 + 1)
    ret += ((xs - 1) // 2 + 1) * ((ys - 1) // 2 + 1) * ((zs - 1) // 2 + 1)
    return ret

# truncate / extend a histogram to 14 buckets
#
def truncate(hist):
    if len(hist) < 15:
        return hist.extend(15 - len(hist))
    hist[14] = sum(hist[14:])
    del hist[15:]
    return hist

# truncate histogram and scale so that sum is N
#
def truncateAndScale(hist):
    hist = truncate(hist)
    n = getN(x_range, y_range, z_range)
    hist = [int(round(h * n)) for h in hist]
    diff = n - sum(hist)
    hist[14] += diff
    assert hist[14] >= 0
    return hist

# Is a lattice point in the range
#
def isIn(point):
    x, y, z = point[0], point[1], point[2]
    if x < x_range[0] or x > x_range[1]:
        return False
    if y < y_range[0] or y > y_range[1]:
        return False
    if z < z_range[0] or z > z_range[1]:
        return False
    return True

# List all neighbors of lattice point
#
def getAllAdj(point):
    x, y, z = point[0], point[1], point[2]
    ret = []
    # corners
    #
    for new_x in [x - 1, x + 1]:
        for new_y in [y - 1, y + 1]:
            for new_z in [z - 1, z + 1]:
                ret.append((new_x, new_y, new_z))
    # faces
    #
    for new_x in [x - 2, x + 2]:
        ret.append((new_x, y, z))
    for new_y in [y - 2, y + 2]:
        ret.append((x, new_y, z))
    for new_z in [z - 2, z + 2]:
        ret.append((x, y, new_z))
    return ret

# Index of point in corresponding neighbor's adjacency list, for all neighbors.
# The corresponding neigbors are as returned by getAllAdj().
# Notice this is the same for all points b/c of symmetry
#
def getRevAdjIdx():
    adj = getAllAdj((0, 0, 0))
    ret = []
    for point in adj:
        adj_adj = getAllAdj(point)
        ret.append(adj_adj.index((0, 0, 0)))
    return ret

revAdjIdx = getRevAdjIdx()
adjDeltas = getAllAdj((0, 0, 0))

def getAdj(point, i):
    adj = (adjDeltas[i][0] + point[0], adjDeltas[i][1] + point[1], adjDeltas[i][2] + point[2])
    return adj

def getInitThroats():
    throats = dict()
    for x in range(x_range[0], x_range[1] + 1, 2):
        for y in range(y_range[0], y_range[1] + 1, 2):
            for z in range(z_range[0], z_range[1] + 1, 2):
                throats[(x, y, z)] = [0] * 14
    for x in range(x_range[0] + 1, x_range[1] + 1, 2):
        for y in range(y_range[0] + 1, y_range[1] + 1, 2):
            for z in range(z_range[0] + 1, z_range[1] + 1, 2):
                throats[(x, y, z)] = [0] * 14
    return throats

def getThroatsN():
    throats = getInitThroats()
    ret = 0
    for point, bits in throats.iteritems():
        for i in range(14):
            if i < revAdjIdx[i]:
                adj = getAdj(point, i)
                if not isIn(adj):
                    continue
                ret += 1
    return ret

def getRandomShuffleThroatIdxs():
    throatIdxs = getInitThroats()
    perm = range(getThroatsN())
    random.shuffle(perm)
    p = 0
    for point, permIdx in throatIdxs.iteritems():
        for i in range(14):
            if i < revAdjIdx[i]:
                adj = getAdj(point, i)
                if not isIn(adj):
                    continue
                permIdx[i] = perm[p]
                throatIdxs[adj][revAdjIdx[i]] = perm[p]
                p += 1
    assert p == len(perm)
    return throatIdxs

throatIdxs = getRandomShuffleThroatIdxs()

def calcThroatCns(throats):
    throatCns = dict()
    for point, bits in throats.iteritems():
        throatCns[point] = sum(bits)
    return throatCns

def calcHist(throatCns):
    hist = [0] * 15
    for point, cn in throatCns.iteritems():
        hist[throatCns[point]] += 1
    return hist

def calcCost(hist, target):
    diffs = [hist[i] - target[i] for i in range(15)]
    return sum([d * d for d in diffs])
    
class Solution():
    throats = dict() # mapping point -> bitvector of throats for each neighbor, as ordered by getAllAdj
    throatCns = dict() # mapping point -> total # of throats
    hist = None # histogram of throatCns
    cost = None # cost of histogram = sum of square difference btwn hist and target
    fit = None # 1 / cost

    def output(self):
        print 'Solution -- cost = ', self.cost, ', hist = ', self.hist, ' fit = ', self.fit
        for point, bits in self.throats.iteritems():
            print point, ' -> (', self.throatCns[point], ' total) ',  bits

    def recalc(self, target):
        self.throatCns = calcThroatCns(self.throats)
        self.hist = calcHist(self.throatCns)
        self.cost = calcCost(self.hist, target)
        self.fit = 1 / self.cost
        self.sanity(target) # TODO enable every now and then

    def __init__(self, target, randomize=False, prob=0):
        self.throats = getInitThroats()
        if randomize:
            for point, bits in self.throats.iteritems():
                for i in range(14):
                    if i < revAdjIdx[i] and random.random() < prob:
                        adj = getAdj(point, i)
                        if not isIn(adj):
                            continue
                        bits[i] = 1
                        self.throats[adj][revAdjIdx[i]] = 1
        self.recalc(target)

    def sanity(self, target):
        assert len(self.throats) == getN(x_range, y_range, z_range)
        assert len(self.throats) == len(self.throatCns)
        for point, bits in self.throats.iteritems():
            assert sum(bits) == self.throatCns[point]
            assert len(bits) == 14
            for i in range(14):
                adj = getAdj(point, i)
                assert bits[i] == 0 or bits[i] == 1
                if not isIn(adj):
                    assert bits[i] == 0
                    assert adj not in self.throats
                    continue
                assert bits[i] == self.throats[adj][revAdjIdx[i]]
            assert self.throatCns[point] <= 14
            assert self.throatCns[point] >= 0
        assert self.hist == calcHist(self.throatCns)
        assert self.cost == calcCost(self.hist, target)
        assert self.fit == 1 / self.cost


def genetic(target):

    print '\n -------------------- Genetic! ---------------------------'

    throatsN = sum([i * target[i] for i in range(15)]) / 2 # approximate b/c of faces but still good
    maxThroatsN = getThroatsN()
    print throatsN, ' / ', maxThroatsN, ' = ', throatsN / maxThroatsN

    populationN = 500 # TODO param
    gensN = 30 # TODO param
    prob = throatsN/maxThroatsN # TODO try with prob = 0.5
    mutationProb = 1 / 10000 # TODO param
    crossoverProb = 0.7 # TODO param
    population = [Solution(target, randomize=True, prob=prob) for _ in range(populationN)]

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
                            if not isIn(adj):
                                continue
                            if throatIdxs[point][b] < switchPos:
                                continue
                            bitMom = mom.throats[point][b]
                            bitDad = dad.throats[point][b]
                            mom.throats[point][b] = bitDad
                            mom.throats[adj][revAdjIdx[b]] = bitDad
                            dad.throats[point][b] = bitMom
                            dad.throats[adj][revAdjIdx[b]] = bitMom
            mom.recalc(target)
            dad.recalc(target)

        # mutate
        #
        for solution in population:
            for point, bits in solution.throats.iteritems():
                for b in range(14):
                    if b < revAdjIdx[b]:
                        adj = getAdj(point, b)
                        if not isIn(adj):
                            continue
                        if random.random() < mutationProb:
                           # print 'flip!'
                            bits[b] = 1 - bits[b]
                            solution.throats[adj][revAdjIdx[b]] = bits[b]
            solution.recalc(target)

        population = new_population

    print 'BEST :  (target = ', target, ')'
    best.output()

if __name__ == '__main__':
    print getN([0, 15 * 2], [0, 15 * 2], [0, 15 * 2])

    adj = getAllAdj((1, 1, 1))
    for point in adj:
        print point
    print range(14)
    print revAdjIdx

    print truncate(ASA_Sandstone_1)
    target = truncateAndScale(ASA_Sandstone_1)
    print target

    print getN(x_range, y_range, z_range)
    solution = Solution(target=target, randomize=True, prob = 0.5)
    print len(solution.throats)
    solution.output()

    #for point, idxs in throatIdxs.iteritems():
    #    print point, ' -> ', idxs

    genetic(target)
