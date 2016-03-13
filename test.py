from __future__ import division

# assumes the corner points are at even coordinates
# and that center points are at odd coordinates
# see https://en.wikipedia.org/wiki/Cubic_crystal_system

x_range = [0, 5]
y_range = [0, 5]
z_range = [0, 5]

# Total number of lattice points in range
#
def getN(x_range, y_range, z_range):
    xs = x_range[1] - x_range[0]
    ys = y_range[1] - y_range[0]
    zs = z_range[1] - z_range[0]
    ret = (xs // 2 + 1) * (ys // 2 + 1) * (zs // 2 + 1)
    ret += ((xs - 1) // 2 + 1) * ((ys - 1) // 2 + 1) * ((zs - 1) // 2 + 1)
    return ret

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

if __name__ == '__main__':
    print 'wtf'
    print getN([0, 4], [0, 4], [0, 4])

    adj = getAllAdj((1, 1, 1))
    for point in adj:
        print point
    print getRevAdjIdx()
