from solution import *

from test import ASA_Sandstone_1 as target

if __name__ == '__main__':

    solution = Solution(target, [0, 1], [0, 1], [0, 1])
    solution.deserialize("best_greedy_0.json")

    print 'cost = ', solution.cost
    print 'hist = ', solution.hist
    print 'target = ', solution.target

    solution.exportForTosho()

    sanity = Solution(solution.target, solution.x_range, solution.y_range, solution.z_range)
    sanity.importFromTosho()

    assert sanity.isEqual(solution)
