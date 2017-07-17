# Find the path with maximum conductance from left to right boundary
# usage: python conductance.py [input file/dir] [output file] [direction]
# example #1: python conductance.py input.csv output.csv 0
# example #2: python conductance.py datadir output.csv 1
# If the input filename ends with '.csv', it is read as the only input file
# Otherwise, it is treated as a directory and all .csv files from that directory and its subdirectories (!) are used as input files
# The output for all files is appended (!) to the output file.
#

import math
import sys
import os
import heapq as heap
import datetime

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
        return pore[3] / ((math.pi * (pore[3] * 0.58) ** 4) / (8 * viscosity))

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
    ucs = 0
    n = 0
    permeability = 0
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


def writeCSV(coords, neigh, left, right, ucs, n, permeability, filename):
    lines = max(len(coords), len(neigh))
    with open(filename, "w") as f:
        for l in range(lines):
            if l >= len(coords):
                f.write("0,0,0,")
            else:
                f.write("%.6lf,%.6lf,%.6lf," % (coords[l][0] * 1e6, coords[l][1] * 1e6, coords[l][2] * 1e6))

            if l >= len(neigh):
                f.write("0,0,0,")
            else:
                f.write("%d,%d,%.6lf," % (neigh[l][0] + 1, neigh[l][1] + 1, neigh[l][2] * 1e6))

            if l >= len(coords):
                f.write("0,")
            else:
                f.write("%.6lf," % (coords[l][3] * 1e6))

            if l >= len(left):
                f.write("0,")
            else:
                f.write("%d," % left[l])

            if l >= len(right):
                f.write("0,")
            else:
                f.write("%d," % right[l])

            if l > 0:
                f.write("0,0,0\n")
            else:
                f.write("%e,%d,%e\n" % (ucs, n, permeability))


def dijkstra(coords, neigh, starting_pore_idxs, ending_pore_idxs):
    # Build adjacency lists
    #
    adj = dict()
    for throat in neigh:
        u = throat[0]
        v = throat[1]
        g = G(throat, coords)
        l = ThTotLen(throat, coords)
        if do_print:
            print 'throat ', u, v, l/g, l, g
        if u in adj:
            adj[u].append((v, l/g, l))
        else:
            adj[u] = [(v, l/g, l)]
        if v in adj:
            adj[v].append((u, l/g, l))
        else:
            adj[v] = [(u, l/g, l)]

    # Dijkstra
    # minimizing sum of l/g, keeping track of total l for finding total g in the end
    #
    pq = []
    [heap.heappush(pq, (0, 0, u)) for u in starting_pore_idxs]

    if do_print:
        print 'starting: '
        for u in starting_pore_idxs:
            print u, ' = ', coords[u]

    visited = set()
    prev = dict()
    dist = dict()
    for u in starting_pore_idxs:
        dist[u] = (0, 0)
        prev[u] = None

    while len(pq) > 0:
        top = heap.heappop(pq)
        lg, l, u = top
        if u in visited:
            continue
        visited.add(u)
        #print 'in ', u, ' with l/g = ', lg, ', l = ', l, ', g = ', l/lg if lg != 0 else 'inf'
        assert dist[u] == (lg, l)
        if u in adj:
            for v, delta_lg, delta_l in adj[u]:
                new_lg = lg + delta_lg
                new_l = l + delta_l
         #       print '        to ', v, ' with l/g = ', new_lg, ', l = ', new_l, ', g = ', new_l/new_lg
                if v not in dist or dist[v][0] > new_lg:
                    dist[v] = (new_lg, new_l)
                    prev[v] = u 
                    heap.heappush(pq, (new_lg, new_l, v))
          #          print '                  add to heap! ', idx, ' --> ', v, ' ', (new_lg, new_l, v)
   
    # Find max g
    #
    gs = []
    for u in ending_pore_idxs:
        if u not in dist:
            continue
        if do_print:
            print 'ending ', u, ' -> ', dist[u], ' prev = ', prev[u]
        lg = dist[u][0]
        l = dist[u][1]
        g = l / lg
        gs.append((g, u))
        if do_print:
            print '               g = ', g

    max_g = max(gs)[0]

    # Find its path
    #
    u = max(gs)[1]
    path = []
    while u is not None:
        path.append(u)
        u = prev[u]
    path = list(reversed(path))

    return max_g, path

def edmondsKarp(coords, neigh, starting_pore_idxs, ending_pore_idxs, doubleVertices=False):
    vertices = coords[:]
    # double edges since it's a non-directred graph
    # note that this is different from adding reverse edges!
    # so in the end, there are 4 edges between each pair of vertices:
    # A -> B with capacity 1 and it's reverse with capacity 0, and
    # B -> A with capacity 1 and it's reverse with capacity 0
    #
    edges = neigh[:] + [(e[1], e[0], e[2]) for e in neigh]
   
    # add source & sink
    #
    source = len(vertices)
    vertices.append((0, 0, 0, 0))
    sink = len(vertices)
    vertices.append((0, 0, 0, 0))
    if do_print:
        print 'souce = ', source,' sink = ', sink
    
    for idx in starting_pore_idxs:
        edges.append((source, idx, 0))
    for idx in ending_pore_idxs:
        edges.append((idx, sink, 0))

    # add rev edges -- edge[i + edges_n] is rev of edge[i]
    #
    edges_n = len(edges)
    cap = [1] * len(edges)
    flow = [0] * len(edges)
    rev = [i + edges_n for i in range(edges_n)] # idx of rev edge
    for idx in range(edges_n):
        u = edges[idx][0]
        v = edges[idx][1]
        edges.append((v, u, 0))
        flow.append(0)
        cap.append(0)
        rev.append(idx)

    # Make capacities to/from source/sink infinite
    #
    for idx in range(len(edges)):
        u = edges[idx][0]
        v = edges[idx][1]
        if u == sink or u == source or v == sink or v == source:
            cap[idx] = 100000000

    # Double vertices -- for each vertex, create in-vertex and out-vertex
    # and have all incoming edges to former and outgoing edges from latter
    # with a capacity 1 edge between in-vertex and out-vertex
    # that way you find the min cut of pores / vertices
    # for convenience, in-vertex is vertices[i] and out-vertex is vertices[i + vertices_n]
    #
    if doubleVertices:
        vertices_n = len(vertices)
        for i in range(vertices_n):
            vertices.append(vertices[i])

        for idx in range(edges_n): # make existing edges u_out -> v_in
            e = edges[idx]
            new_edge = (e[0] + vertices_n, e[1], e[2])
            edges[idx] = new_edge

            e = edges[idx + edges_n] # reverse edge
            assert e[1] + vertices_n == new_edge[0]
            assert e[0] == new_edge[1]
            new_edge = (e[0], e[1] + vertices_n, e[2])
            edges[idx + edges_n] = new_edge

        source += vertices_n # source is out-vertex

        for u in range(vertices_n): # add in-out edges
            v = u + vertices_n

            edges.append((u, v, 0))
            flow.append(0)
            cap.append(1)
            rev.append(len(edges))

            edges.append((v, u, 0)) # and reverse edge
            flow.append(0)
            cap.append(0)
            rev.append(len(edges) - 2)

    # Build adjacency lists -- adj[v] = list of outgoing edges from v
    #
    adj = dict()
    for idx in range(len(edges)):
        u = edges[idx][0]
        if u in adj:
            adj[u].append(idx)
        else:
            adj[u] = [idx]

    # BFS
    #
    def BFS():
        foundPath = False
        prev = [None] * len(vertices)
        prev_edge = [None] * len(vertices)
        visited = [0] * len(vertices)
        visited[source] = True
        pq = [source]
        head = 0
        while head < len(pq) and not foundPath:
            u = pq[head]
            head += 1
            if do_print:
                print '   in ', u
            for idx in adj[u]:
                e = edges[idx]
                v = e[1]
                cf = cap[idx] - flow[idx]
                if do_print:
                    print '         to ', v, ' flow = ', flow[idx], ' cap = ', cap[idx], ' cf = ', cf, ' visited ', visited[v]
                if cf > 0 and not visited[v]:
                    if do_print:
                        print '                       go!'
                    pq.append(v)
                    visited[v] = True
                    prev[v] = u
                    prev_edge[v] = idx
                    if v == sink:
                        if do_print:
                            print '                           SINK!'
                        foundPath = True
                        break
        return foundPath, prev, prev_edge, visited

    # Keep finding augmenting paths
    #
    while True:
        if do_print:
            print '----------- iteration ----------------'

        foundPath, prev, prev_edge, _ = BFS()
        if not foundPath:
            # no more augmeting paths
            #
            break

        # Reconstruct path
        #
        v = sink
        path = []
        path_edges = []
        while v != source:
            path.append(v)
            path_edges.append(prev_edge[v])
            v = prev[v]
        path.append(source)
        path = list(reversed(path))
        path_edges = list(reversed(path_edges))
        if do_print:
            print '    augmenting path = ', path
            print '    augmenting path edges = ', path_edges

        # Find min augmented capacity
        #
        min_cf = None
        for idx in path_edges:
            u = edges[idx][0]
            v = edges[idx][1]
            cf = cap[idx] - flow[idx]
            assert cf > 0
            if do_print:
                print '              min_cf :  edge ', idx, ' = ', u, v, ' cf = ', cf
            if min_cf is None or min_cf > cf:
                min_cf = cf
        if do_print:
            print '     Min cf = ', min_cf

        # Augment path
        #
        for idx in path_edges:
            u = edges[idx][0]
            v = edges[idx][1]
            ridx = rev[idx]
            flow[idx] += min_cf
            flow[ridx] -= min_cf
            assert edges[ridx][0] == v
            assert edges[ridx][1] == u
            assert flow[idx] <= cap[idx]
            assert flow[ridx] <= cap[ridx]
            if do_print:
                print '              augment :  edge ', idx, ' = ', u, v, ' flow = ', flow[idx], ' r flow (for ', ridx, ') = ', flow[ridx]

    # calculate the max flow
    #
    max_flow_source = 0
    max_flow_sink = 0
    for idx in range(len(edges)):
        if edges[idx][0] == source:
            if do_print:
                print '    max flow source ', idx, ' (', source, edges[idx][1], ') += ', flow[idx]
            max_flow_source += flow[idx]
        if edges[idx][1] == sink:
            if do_print:
                print '    max flow sink ', idx, ' (', edges[idx][1], sink, ') += ', flow[idx]
            max_flow_sink += flow[idx]

    # find the critical throats
    #
    foundPath, _, _, visited = BFS()
    critical = []
    assert not foundPath
    for idx in range(len(edges)):
        u = edges[idx][0]
        v = edges[idx][1]
        cf = cap[idx] - flow[idx]
        if (u, v) in critical or (v, u) in critical:
            continue
        if doubleVertices and v != u + vertices_n:
            # we're looking at critical pores => only consider the in-out edges
            #
            continue
        if cf == 0 and visited[u] != visited[v]:
            critical.append((u, v))
            if do_print:
                print '                     critical edge = ', u, v
    print 'Total = ', len(critical), ' critical'

    print 'MAX FLOW source = ', max_flow_source, ' sink (should be same) = ', max_flow_sink
    assert max_flow_source == max_flow_sink
    assert len(critical) == max_flow_source
    return max_flow_source, critical

# In case they're not given
# direction = 0 => X, 1 => Y, 2 => Z
#
def getBoundaries(coords, direction):
    starting = []
    ending = []
    min_b = min([pore[direction] for pore in coords])
    max_b = max([pore[direction] for pore in coords])
    left_b = [pore[direction] - pore[3] for pore in coords]
    right_b = [pore[direction] + pore[3] for pore in coords]
    for i in range(len(coords)):
        if left_b[i] <= min_b:
            starting.append(i)
        if right_b[i] >= max_b:
            ending.append(i)
    return starting, ending

def solve(infile, outfile, direction):
    coords, neigh, left, right, ucs, n, permeability = readCSV(infile)
    print '\n\n--------------------- solving', infile, '---------------------------\n\n'
    print 'N = ', n
    print 'UCS = ', ucs
    print 'Perm = ', permeability
    print '# pores = ', len(coords)
    print '# throats = ', len(neigh)
    print '# pores on left boundary = ', len(left)
    print '# pores on right boundary = ', len(right)

    if len(left) == 0 or len(right) == 0:
        left, right = getBoundaries(coords, direction=direction)
        print 'Computing boundaries: # left = ', len(left), ', # right = ', len(right)

    # Find the critical throats
    #
    print '\nFinding critical throats...'
    then = datetime.datetime.now()
    max_flow_throats, critical_throats = edmondsKarp(coords, neigh, left, right, doubleVertices=False)
    print 'Edmonds Karp took ', (datetime.datetime.now() - then).total_seconds(), ' sec'
    assert len(critical_throats) == max_flow_throats
    print 'max flow (throats) = ', max_flow_throats
    print 'critical edges (throats) = ', critical_throats

    # Find the critical pores
    #
    print '\nFinding critical pores...'
    max_flow_pores, critical_pores = edmondsKarp(coords, neigh, left, right, doubleVertices=True)
    assert len(critical_pores) == max_flow_pores
    print 'max flow (pores) = ', max_flow_pores
    print 'critical edges (pores) = ', critical_pores

    # Find the path of max conductance... or preferential flow path... TODO figure out which one
    #
    print '\nRunning dijkstra....'
    max_g, path = dijkstra(coords, neigh, left, right)
    print 'MAX g = ', max_g
    print 'Path = ', path
    assert path[0] in left
    assert path[-1] in right

    # TODO sanity path calc make sure it's same as max g

    with open(outfile, 'a') as f:
        res = [infile, permeability, max_g, max_g / (n * ucs), ' '.join([str(v + 1) for v in path])]
        res.append(len(critical_throats))
        res.append(' '.join(['%d-%d' % (e[0] + 1, e[1] + 1) for e in critical_throats]))
        res.append(len(critical_pores))
        res.append(' '.join([str(e[0] + 1) for e in critical_pores]))

        f.write(','.join([str(x) for x in res]) + '\n')

if __name__ == '__main__':
    infile = sys.argv[1]
    outfile = sys.argv[2]
    direction = int(sys.argv[3])

    if infile.lower().endswith('.csv'):
        # single file
        #
        solve(infile, outfile, direction)
    else:
        # directory of files
        #
        dirname = infile 
        for (path, dirs, files) in os.walk(dirname):
            print '\n============ EXPLORING DIRECTORY', path, '=================\n'
            for filename in files:
                if filename.lower().endswith('.csv'):
                    solve(os.path.join(path, filename), outfile, direction)
