################################
#
# Avi Schwarzschild
# Fall 2015 - reaserch
# misc.py for useful functions
# 
###############################
import numpy as np

def cell_values(grid, interpolants, integrals):
        values = []
        for i in xrange(len(grid)-1):
            j = 0
            x0found = x1found = False
            while not(x0found and x1found):
                ints = interpolants[j]
                if grid[i] >= ints[0]: 
                    i0 = j
                    x0found = True
                if grid[i+1] <= ints[1]: 
                    i1 = j
                    x1found = True
                j += 1
            sum = 0
            for k in xrange(i0, i1+1):
                if interpolants[k][0] >= grid[i]: 
                    a = interpolants[k][0]
                else: 
                    a = grid[i]
                if interpolants[k][1] <= grid[i+1]:
                    b = interpolants[k][1] 
                else:
                    b = grid[i+1]
                p = integrals[k](b) - integrals[k](a)
                sum += p
            values.append(sum/(grid[i+1]-grid[i]))
        return np.array(values)