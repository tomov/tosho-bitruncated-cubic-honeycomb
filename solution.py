from __future__ import division
import random

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
def truncateAndScale(hist, x_range, y_range, z_range):
    hist = truncate(hist)
    n = getN(x_range, y_range, z_range)
    hist = [int(round(h * n)) for h in hist]
    diff = n - sum(hist)
    hist[14] += diff
    assert hist[14] >= 0
    return hist

# Is a lattice point in the range
#
def isIn(point, x_range, y_range, z_range):
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

def getInitThroats(x_range, y_range, z_range):
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

def getThroatsN(x_range, y_range, z_range):
    throats = getInitThroats(x_range, y_range, z_range)
    ret = 0
    for point, bits in throats.iteritems():
        for i in range(14):
            if i < revAdjIdx[i]:
                adj = getAdj(point, i)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                ret += 1
    return ret

def getRandomShuffleThroatIdxs(x_range, y_range, z_range):
    throatIdxs = getInitThroats(x_range, y_range, z_range)
    perm = range(getThroatsN(x_range, y_range, z_range))
    random.shuffle(perm)
    p = 0
    for point, permIdx in throatIdxs.iteritems():
        for i in range(14):
            if i < revAdjIdx[i]:
                adj = getAdj(point, i)
                if not isIn(adj, x_range, y_range, z_range):
                    continue
                permIdx[i] = perm[p]
                throatIdxs[adj][revAdjIdx[i]] = perm[p]
                p += 1
    assert p == len(perm)
    return throatIdxs

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
    x_range = None
    y_range = None
    z_range = None

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

    def __init__(self, target, x_range, y_range, z_range, randomize=False, prob=0):
        self.throats = getInitThroats(x_range, y_range, z_range)
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        if randomize:
            for point, bits in self.throats.iteritems():
                for i in range(14):
                    if i < revAdjIdx[i] and random.random() < prob:
                        adj = getAdj(point, i)
                        if not isIn(adj, x_range, y_range, z_range):
                            continue
                        bits[i] = 1
                        self.throats[adj][revAdjIdx[i]] = 1
        self.recalc(target)

    def set(self, point, b, value):
        assert isIn(point, self.x_range, self.y_range. self.z_range)
        adj = getAdj(point, b)
        assert isIn(adj, self.x_range, self.y_range. self.z_range)
        self.throats[point][b] = value
        self.throats[adj][revAdjIdx[b]] = value

    def sanity(self, target):
        assert len(self.throats) == getN(self.x_range, self.y_range, self.z_range)
        assert len(self.throats) == len(self.throatCns)
        for point, bits in self.throats.iteritems():
            assert sum(bits) == self.throatCns[point]
            assert len(bits) == 14
            for i in range(14):
                adj = getAdj(point, i)
                assert bits[i] == 0 or bits[i] == 1
                if not isIn(adj, self.x_range, self.y_range, self.z_range):
                    assert bits[i] == 0
                    assert adj not in self.throats
                    continue
                assert bits[i] == self.throats[adj][revAdjIdx[i]]
            assert self.throatCns[point] <= 14
            assert self.throatCns[point] >= 0
        assert self.hist == calcHist(self.throatCns)
        assert self.cost == calcCost(self.hist, target)
        assert self.fit == 1 / self.cost
