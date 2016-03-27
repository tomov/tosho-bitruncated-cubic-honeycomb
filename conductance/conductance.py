import math
import heapq as heap

coords = None # array of pores = (pore X, pore Y, pore Z, pore radius in m)
neigh = None  # array of throats = (throat pore 1 index in coords, throat pore 2, throat radius in m)

def ThTotLen(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    return math.sqrt((pore1[0] - pore2[0]) ** 2 + (pore1[1] - pore2[1]) ** 2 + (pore1[2] - pore2[2]) ** 2)

def ThLen(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    ret = ThTotLen(throat, coords) - pore1[3] - pore2[3]
    return ret

def G(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    viscosity = 0.001 # Pa * s
    def addend(pore):
        return pore[3] / ((math.pi * (pore[3] * 2/3) ** 4) / (8 * viscosity))

    addend_throat = ThLen(throat, coords) / ((math.pi * throat[2] ** 4) / (8 * viscosity))
    denominator = addend(pore1) + addend(pore2) + addend_throat
    return ThTotLen(throat, coords) / denominator

def readCSV(filename):
    coords = []
    neigh = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(',')
            if line[0] != '':
                pore = (float(line[0]) * 1e-6, float(line[1]) * 1e-6, float(line[2]) * 1e-6, float(line[6]) * 1e-6)
                coords.append(pore)
            throat = (int(line[3]) - 1, int(line[4]) - 1, float(line[5]) * 1e-6)
            neigh.append(throat)
    return coords, neigh

def dijkstra(coords, neigh, starting_pore_idxs, ending_pore_idxs):
    # Build adjacency lists
    #
    adj = dict()
    for throat in neigh:
        idx1 = throat[0]
        idx2 = throat[1]
        g = G(throat, coords)
        l = ThTotLen(throat, coords)
        print 'throat ', idx1, idx2, l/g, l, g
        if idx1 in adj:
            adj[idx1].append((idx2, l/g, l))
        else:
            adj[idx1] = [(idx2, l/g, l)]
        if idx2 in adj:
            adj[idx2].append((idx1, l/g, l))
        else:
            adj[idx2] = [(idx1, l/g, l)]

    # Dijkstra
    # minimizing sum of l/g, keeping track of total l for finding total g in the end
    #
    pq = []
    [heap.heappush(pq, (0, 0, idx)) for idx in starting_pore_idxs]

    print 'starting: '
    for idx in starting_pore_idxs:
        print idx, ' = ', coords[idx]

    visited = set()
    prev = dict()
    dist = dict()
    for idx in starting_pore_idxs:
        dist[idx] = (0, 0)

    while len(pq) > 0:
        top = heap.heappop(pq)
        lg, l, idx = top
        if idx in visited:
            continue
        visited.add(idx)
        #print 'in ', idx, ' with l/g = ', lg, ', l = ', l, ', g = ', l/lg if lg != 0 else 'inf'
        assert dist[idx] == (lg, l)
        if idx in adj:
            for neighbor, delta_lg, delta_l in adj[idx]:
                new_lg = lg + delta_lg
                new_l = l + delta_l
         #       print '        to ', neighbor, ' with l/g = ', new_lg, ', l = ', new_l, ', g = ', new_l/new_lg
                if neighbor not in dist or dist[neighbor][0] > new_lg:
                    dist[neighbor] = (new_lg, new_l)
                    prev[neighbor] = idx
                    heap.heappush(pq, (new_lg, new_l, neighbor))
          #          print '                  add to heap! ', idx, ' --> ', neighbor, ' ', (new_lg, new_l, neighbor)
    
    gs = []
    for idx in ending:
        if idx not in dist:
            continue
        print 'ending ', idx, ' -> ', dist[idx]
        lg = dist[idx][0]
        l = dist[idx][1]
        g = l / lg
        gs.append(g)
        print '               g = ', g

    return max(gs)

def WTFxsWhy():
    xs = dict()
    for pore in coords:
        x = pore[0]
        if x in xs:
            xs[x] += 1
        else:
            xs[x] = 1

    xs = [(x, c) for x, c in xs.iteritems()]
    xs = sorted(xs)
    prev = 0
    for x, c in xs:
        print x, c, '           ', x - prev
        prev = x

if __name__ == '__main__':
    #coords, neigh = readCSV('Sand1Ori.csv')
    coords, neigh = readCSV('Sand1N10Ind9.csv')

    starting = []
    ending = []
    min_x = min([pore[0] for pore in coords])
    max_x = max([pore[0] for pore in coords])
    for i in range(len(coords)):
        pore = coords[i]
        if pore[0] == min_x:
            starting.append(i)
        if pore[0] == max_x:
            ending.append(i)
    max_g = dijkstra(coords, neigh, starting, ending)
    print 'MAX g = ', max_g
