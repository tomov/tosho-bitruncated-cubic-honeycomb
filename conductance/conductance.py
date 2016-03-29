# Find the path with maximum conductance from left to right boundary
# usage: python conductance.py [input file/dir] [output file]
# example #1: python conductance.py input.csv output.csv
# example #2: python conductance.py datadir output.csv
# If the input filename ends with '.csv', it is read as the only input file
# Otherwise, it is treated as a directory and all .csv files from that directory and its subdirectories (!) are used as input files
# The output for all files is appended (!) to the output file.
#

import math
import sys
import os
import heapq as heap

"""
coords # array of pores = (pore X, pore Y, pore Z, pore radius in m)
neigh  # array of throats = (throat pore 1 index in coords, throat pore 2, throat radius in m)
left # idxs of pores on left boundary
right # idxs of pores on right boundary
ucs # size of single cell (???)
n # number of pores along axis
permeability # as computed by the complicated formulas
"""

do_print = False

def ThTotLen(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    return math.sqrt((pore1[0] - pore2[0]) ** 2 + (pore1[1] - pore2[1]) ** 2 + (pore1[2] - pore2[2]) ** 2)

def ThLen(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    ret = ThTotLen(throat, coords) - pore1[3] - pore2[3]
    if ret < 0:
        return 1e-6
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
    # Row format:
    #  0 1 2  3     4     5   6  7  8  9  10  11 
    # [X,Y,Z,Pore1,Pore2,ThR,Pr,LB,RB,UCS,N,Perm] 
    #
    coords = []
    neigh = []
    left = []
    right = []
    ucs = None
    n = None
    permeability = None
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(',')
            line = [l.strip() for l in line]
            if line[0] != '' and line[0] != '0':
                pore = (float(line[0]) * 1e-6, float(line[1]) * 1e-6, float(line[2]) * 1e-6, float(line[6]) * 1e-6)
                coords.append(pore)
            if line[3] != '' and line[3] != '0':
                throat = (int(line[3]) - 1, int(line[4]) - 1, float(line[5]) * 1e-6)
                neigh.append(throat)
            if line[7] != '' and line[7] != '0':
                idx = int(line[7]) - 1
                left.append(idx)
            if line[8] != '' and line[8] != '0':
                idx = int(line[8]) - 1
                right.append(idx)
            if line[9] != '' and line[9] != '0':
                ucs = float(line[9])
            if line[10] != '' and line[10] != '0':
                n = int(line[10])
            if line[11] != '' and line[11] != '0':
                permeability = float(line[11])

    return coords, neigh, left, right, ucs, n, permeability

def dijkstra(coords, neigh, starting_pore_idxs, ending_pore_idxs):
    # Build adjacency lists
    #
    adj = dict()
    for throat in neigh:
        idx1 = throat[0]
        idx2 = throat[1]
        g = G(throat, coords)
        l = ThTotLen(throat, coords)
        if do_print:
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

    if do_print:
        print 'starting: '
        for idx in starting_pore_idxs:
            print idx, ' = ', coords[idx]

    visited = set()
    prev = dict()
    dist = dict()
    for idx in starting_pore_idxs:
        dist[idx] = (0, 0)
        prev[idx] = None

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
   
    # Find max g
    #
    gs = []
    for idx in ending_pore_idxs:
        if idx not in dist:
            continue
        if do_print:
            print 'ending ', idx, ' -> ', dist[idx], ' prev = ', prev[idx]
        lg = dist[idx][0]
        l = dist[idx][1]
        g = l / lg
        gs.append((g, idx))
        if do_print:
            print '               g = ', g

    max_g = max(gs)[0]

    # Find its path
    #
    idx = max(gs)[1]
    path = []
    while idx is not None:
        path.append(idx)
        idx = prev[idx]
    path = list(reversed(path))

    return max_g, path

# In case they're not given
#
def getBoundaries(coords):
    starting = []
    ending = []
    min_x = min([pore[0] for pore in coords])
    max_x = max([pore[0] for pore in coords])
    left_x = [pore[0] - pore[3] for pore in coords]
    right_x = [pore[0] + pore[3] for pore in coords]
    for i in range(len(coords)):
        if left_x[i] <= min_x:
            starting.append(i)
        if right_x[i] >= max_x:
            ending.append(i)
    return starting, ending

def solve(infile, outfile):
    coords, neigh, left, right, ucs, n, permeability = readCSV(infile)
    print '---------- solving', infile, '--------------'
    print 'N = ', n
    print 'UCS = ', ucs
    print 'Perm = ', permeability
    print '# pores = ', len(coords)
    print '# throats = ', len(neigh)
    print '# pores on left boundary = ', len(left)
    print '# pores on right boundary = ', len(right)

    print 'Running dijkstra....'
    max_g, path = dijkstra(coords, neigh, left, right)
    print 'MAX g = ', max_g
    print 'Path = ', path
    assert path[0] in left
    assert path[-1] in right

    # TODO sanity path calc make sure it's same as max g

    with open(outfile, 'a') as f:
        res = [infile, permeability, max_g, max_g / (n * ucs), ' '.join([str(p + 1) for p in path])]
        f.write(','.join([str(x) for x in res]) + '\n')

if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]

    if infile.lower().endswith('.csv'):
        # single file
        #
        solve(infile, outfile)
    else:
        # directory of files
        #
        dirname = infile 
        for (path, dirs, files) in os.walk(dirname):
            print '\n============ EXPLORING DIRECTORY', path, '=================\n'
            for filename in files:
                if filename.lower().endswith('.csv'):
                    solve(os.path.join(path, filename), outfile)
