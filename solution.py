# Represens a single solution to the distribution
# assumes the corner points are at even coordinates
# and that center points are at odd coordinates
# see https://en.wikipedia.org/wiki/Cubic_crystal_system

from __future__ import division
import random
import json
import math
import numpy as np


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
# assumes sum is 1
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

    #random.shuffle(ret, lambda: 0.2342) # TODO ARGHHH WTF
    #ret.reverse()
    #print 'HERE BITCH'
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
        idx = adj_adj.index((0, 0, 0))
        ret.append(idx)
    for i in range(len(ret)):
        assert ret[ret[i]] == i
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
            adj = getAdj(point, i)
            if point > adj:
                continue
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
    target = None

    throats = dict() # mapping point -> bitvector of throats for each neighbor, as ordered by getAllAdj
    throatCns = dict() # mapping point -> total # of throats
    hist = None # histogram of throatCns
    cost = None # cost of histogram = sum of square difference btwn hist and target
    fit = None # 1 / cost

    def isEqual(self, other):
        if self.x_range != other.x_range:
            return False
        if self.y_range != other.z_range:
            return False
        if self.z_range != other.y_range:
            return False
        if self.target != other.target:
            return False
        if self.throats != other.throats:
            return False
        if self.throatCns != other.throatCns:
            return False
        if self.hist != other.hist:
            return False
        if self.cost != other.cost:
            return False
        if self.fit != other.fit:
            return False
        return True

    def serialize(self, filename):
        self.sanity()
        with open(filename, "w") as f:
            d = {
                'x_range': self.x_range,
                'y_range': self.y_range,
                'z_range': self.z_range,
                'target': self.target,
                'throats_keys': self.throats.keys(),
                'throats_values': self.throats.values()
            }
            f.write(json.dumps(d))

    def deserialize(self, filename):
        with open(filename, "r") as f:
            s = f.read()
            d = json.loads(s)
            self.x_range = d['x_range']
            self.y_range = d['y_range']
            self.z_range = d['z_range']
            self.target = d['target']
            throats_keys = d['throats_keys']
            throats_values = d['throats_values']
            throats_keys = [tuple(k) for k in throats_keys]
            self.throats = dict(zip(throats_keys, throats_values))
            self.recalc()

    def output(self):
        print 'Solution -- cost = ', self.cost, ', hist = ', self.hist, ' fit = ', self.fit
        for point, bits in self.throats.iteritems():
            print point, ' -> (', self.throatCns[point], ' total) ',  bits

    def exportForTosho(self, prefix = "", suffix = ""):
        self.sanity()
        coords = []
        neigh = []
        neigh_WTF = []
        pointIdx = dict()
        for point, bits in self.throats.iteritems():
            pointIdx[point] = len(coords)
            coords.append([point[0], point[1], point[2]])
        WTF = dict()
        for point, bits in self.throats.iteritems():
            for i in range(14):
                if i < revAdjIdx[i]:
                    adj = getAdj(point, i)
                   # if point > adj:
                   #     continue
                    if not isIn(adj, self.x_range, self.y_range, self.z_range):
                        continue
                    if bits[i]:
                        idxA = pointIdx[point]
                        idxB = pointIdx[adj]
                        if idxA < 0 or idxB < 0:
                            print 'SHIT'
                            print idxA, idxB
                            assert False
                        if (idxA, idxB) in WTF:
                            continue
                        WTF[(idxA, idxB)] = 1
                        WTF[(idxB, idxA)] = 1
                        neigh.append([idxA, idxB, 1])
        for point, bits in self.throats.iteritems():
            for i in range(14):
                    adj = getAdj(point, i)
                    #if point > adj:
                    #    continue
                    if not isIn(adj, self.x_range, self.y_range, self.z_range):
                        continue
                    if bits[i]:
                        idxA = pointIdx[point]
                        idxB = pointIdx[adj]
                        assert WTF[(idxA, idxB)]
                        assert WTF[(idxB, idxA)]
        used = dict()
        for shit, _ in WTF.iteritems():
            if (shit[0], shit[1]) not in used:
                neigh_WTF.append([shit[0], shit[1], 1])
                used[(shit[0], shit[1])] = True
                used[(shit[1], shit[0])] = True
        print len(neigh)
        print len(neigh_WTF)
        print sum([i * self.hist[i] for i in range(15)]) / 2
        assert len(neigh) == sum([i * self.hist[i] for i in range(15)]) / 2
        assert len(neigh_WTF) == sum([i * self.hist[i] for i in range(15)]) / 2
        for point in coords: # my coords are * 2 and from 0
            point[0] = point[0] / 2 + 100
            point[1] = point[1] / 2 + 100
            point[2] = point[2] / 2 + 100
        for n in neigh: # matlab starts from 1
            n[0] += 1
            n[1] += 1
        for n in neigh_WTF: # matlab starts from 1
            n[0] += 1
            n[1] += 1
        coords_file = "%scoords%s.csv" % (prefix, suffix)
        with open(coords_file, "w") as f:
            for point in coords:
                f.write("%f,%f,%f\n" % (point[0], point[1], point[2]))
        neigh_file = "%sneigh%s.csv" % (prefix, suffix)
        with open(neigh_file, "w") as f:
            for n in neigh:
                f.write("%d,%d,%d\n" % (n[0], n[1], n[2]))
        neigh_WTF_file = "%sneigh%s_WTF.csv" % (prefix, suffix)
        with open(neigh_WTF_file, "w") as f:
            for n in neigh_WTF:
                f.write("%d,%d,%d\n" % (n[0], n[1], n[2]))
        return coords_file, neigh_file, neigh_WTF_file

    # requires that ranges and target are initialized accordingly
    #
    def importFromTosho(self, prefix = "", suffix = ""):
        self.throats = getInitThroats(self.x_range, self.y_range, self.z_range)
        pointIdx = dict()
        points = []
        with open("%scoords%s.csv" % (prefix, suffix), "r") as f:
            idx = 0
            for line in f.readlines():
                c = [int((float(x) - 100) * 2) for x in line.split(",")]
                assert len(c) == 3
                point = (c[0], c[1], c[2])
                assert point not in pointIdx
                assert isIn(point, self.x_range, self.y_range, self.z_range)
                pointIdx[point] = idx
                points.append(point)
                idx += 1
        with open("%sneigh%s.csv" % (prefix, suffix), "r") as f:
            for line in f.readlines():
                n = [int(x) for x in line.split(",")]
                pointA = points[n[0] - 1]
                pointB = points[n[1] - 1]
                assert pointA < pointB
                for i in range(14):
                    adj = getAdj(pointA, i)
                    if adj == pointB:
                        self.throats[pointA][i] = 1
                        self.throats[adj][revAdjIdx[i]] = 1
                        break
        self.recalc()

    def recalc(self):
        self.throatCns = calcThroatCns(self.throats)
        self.hist = calcHist(self.throatCns)
        self.cost = calcCost(self.hist, self.target)
        self.fit = 1 / self.cost
        self.sanity() # TODO enable every now and then

    def __init__(self, target, x_range, y_range, z_range, randomize=False, prob=0):
        self.throats = getInitThroats(x_range, y_range, z_range)
        self.x_range = x_range
        self.y_range = y_range
        self.z_range = z_range
        self.target = target
        if randomize:
            points = self.throats.keys()
            all_idx = range(len(points))
            random.shuffle(all_idx)
            for idx in all_idx:
                point = points[idx]
                bits = self.throats[point]
                bit_idxs = range(14)
                random.shuffle(bit_idxs)
                # WTF TODO WHYYYY THIS DOESN'T WORK
                rand_bits = np.random.choice(range(14), 4, replace=False)
                print getAllAdj((0, 0, 0))
                for i in range(14):
                    if i < revAdjIdx[i] and random.random() < prob:
                        adj = getAdj(point, i)
                        #if point > adj:
                        #    continue
                        print point, ' --> ', adj
                        if not isIn(adj, x_range, y_range, z_range):
                            continue
                        bits[i] = 1
                        self.throats[adj][revAdjIdx[i]] = 1
        self.recalc()

    def set(self, point, b, value):
        assert value == 0 or value == 1
        assert isIn(point, self.x_range, self.y_range, self.z_range)
        adj = getAdj(point, b)
        assert isIn(adj, self.x_range, self.y_range, self.z_range)
        self.throats[point][b] = value
        self.throats[adj][revAdjIdx[b]] = value

    def __hist_deltas(self, point, b, value, adj):
        deltas = dict()
        for p in [point, adj]:
            old_bucket = self.throatCns[p]
            if old_bucket in deltas:
                deltas[old_bucket] -= 1
            else:
                deltas[old_bucket] = -1
            new_bucket = old_bucket + (1 if value else -1)
            assert new_bucket >= 0 and new_bucket <= 14
            if new_bucket in deltas:
                deltas[new_bucket] += 1
            else:
                deltas[new_bucket] = 1
        return deltas


    def setAndRecalc(self, point, b, value):
        old_value = self.throats[point][b]
        if old_value == value:
            return
        self.set(point, b, value)
        adj = getAdj(point, b)
        assert isIn(adj, self.x_range, self.y_range, self.z_range)

        # get hist deltas
        deltas = self.__hist_deltas(point, b, value, adj)

        # update throatCns AFTER getting hist deltas
        if old_value:
            self.throatCns[point] -= 1
            self.throatCns[adj] -= 1
        else:
            self.throatCns[point] += 1
            self.throatCns[adj] += 1

        # update cost
        for bucket, delta in deltas.iteritems():
            self.cost -= (self.hist[bucket] - self.target[bucket]) ** 2
            self.cost += (self.hist[bucket] + delta - self.target[bucket]) ** 2
        self.fit = 1 / self.cost

        # update hists AFTER update cost
        for bucket, delta in deltas.iteritems():
            self.hist[bucket] += delta

        #self.sanity()


    def costIfSet(self, point, b, value):
        assert isIn(point, self.x_range, self.y_range, self.z_range)
        if self.throats[point][b] == value:
            return self.cost
        adj = getAdj(point, b)
        assert isIn(adj, self.x_range, self.y_range, self.z_range)

        deltas = self.__hist_deltas(point, b, value, adj)

        new_cost = self.cost
        for bucket, delta in deltas.iteritems():
            new_cost -= (self.hist[bucket] - self.target[bucket]) ** 2
            new_cost += (self.hist[bucket] + delta - self.target[bucket]) ** 2
        return new_cost

    def sanity(self):
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
        assert self.cost == calcCost(self.hist, self.target)
        assert self.fit == 1 / self.cost

