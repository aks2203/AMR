################################
#
# Avi Schwarzschild
# Fall 2015 - reaserch
# bathymetry functions
# 
###############################

import random
from scipy.interpolate import lagrange
import numpy as np
import matplotlib.pyplot as plt
from misc import cell_values


class Bathymetry:
    
    def __init__(self, bath_data):
        self.bath_data = bath_data
        self.interpolants = []
        self.integrals = []
        for j in xrange(len(bath_data)-1):
            li = lagrange([bath_data[i,0] for i in (j,j+1)], [bath_data[i,1] for i in (j,j+1)])
            self.interpolants.append((bath_data[j, 0], bath_data[j+1, 0], li))
        for i in self.interpolants:
            self.integrals.append(np.polyint(i[2]))
        
    def bath_cell_values(self, grid):
        return cell_values(grid, self.interpolants, self.integrals)