# Takes coords and neigh file (as output by run_tosho.py)
# and merges them in a csv file in the format accepted by conductance.py
# Usage: python merge_outs.py coords.csv neigh.csv conductance/shit.csv
#
# CSV format: [X,Y,Z,Pore1,Pore2,ThR,Pr,LB,RB,UCS,N,Perm]

import sys

def merge_outs(coords_file, neigh_file, out_file):
    with open(coords_file, 'r') as f:
        coords = f.readlines()

    with open(neigh_file, 'r') as f:
        neigh = f.readlines()

    lines = max(len(coords), len(neigh))
    with open(out_file, 'w') as f:
        for l in range(lines):
            z = []
            if l >= len(coords):
                z.extend([0, 0, 0])
            else:
                z.extend([x.strip() for x in coords[l].split(',')])
            if l >= len(neigh):
                z.extend([0, 0, 0])
            else:
                z.extend([x.strip() for x in neigh[l].split(',')])
            if l == 0:
                z.extend([1, 0, 0, 1, 1, 1])
            else:
                z.extend([1, 0, 0, 0, 0, 0])
            f.write(','.join([str(x) for x in z]) + '\n')

if __name__ == "__main__":
    coords_file = sys.argv[1]
    neigh_file = sys.argv[2]
    out_file = sys.argv[3]

    merge_outs(coords_file, neigh_file, out_file)
