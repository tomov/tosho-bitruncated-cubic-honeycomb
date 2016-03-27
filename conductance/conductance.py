import math

def ThTotLen(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    return math.sqrt((pore1[0] - pore2[0]) ** 2 + (pore1[1] - pore2[1]) ** 2 + (pore1[2] - pore2[2]) ** 2)

def ThLen(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    ret = ThTotLen(throat, coords) - pore1[3] - pore2[3]
    assert ret >= 0
    return ret

def G(throat, coords):
    pore1 = coords[throat[0]]
    pore2 = coords[throat[1]]
    ret = ThTotLen(throat, coords) / 

    ThTotLen(i)
       /
       ( 
          PR(neigh(i,1)) 
          / (
              (
                pi * (
                       PR(neigh(i,1)) * 2/3
                     )^4
              )
              /
              ( 8 * viscosity )
            )

          +

          PR(neigh(i,2))
          / (
              (
                pi * (
                       PR(neigh(i,2)) * 2/3
                     )^4
              )
              /
              ( 8 * viscosity )
            )
          
          +
          
          ThLen(i)
          /
          (
            ( pi * neigh(i,3)^4 )
            / ( 8*viscosity ) 
          )
       )

def readCSV(filename):
    coords = []
    neigh = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.split(',')
            if line[0] != '':
                pore = (float(line[0]), float(line[1]), float(line[2]), float(line[6]))
                coords.append(pore)
            throat = (int(line[3]) - 1, int(line[4]) - 1, float(line[5]))
            neigh.append(throat)
    return coords, neigh

if __name__ == '__main__':
    coords, neigh = readCSV('Sand1Ori.csv')

    print ThTotLen(neigh[0], coords)
