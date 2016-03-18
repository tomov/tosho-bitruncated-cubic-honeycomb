from solution import *

solution = Solution([0] * 15, [0, 1], [0, 1], [0, 1])

solution.deserialize('best_greedy.json')

print 'cost = ', solution.cost
print 'target = ', solution.target
print 'solution = ', solution.hist
print 'x = ', solution.x_range
print 'y = ', solution.y_range
print 'z = ', solution.z_range
