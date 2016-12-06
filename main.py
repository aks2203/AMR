################################
#
# Avi Schwarzschild
# Fall 2016 - reaserch 
# main2.py
# 
###############################


from Grid import Grid
from BathymetryClass import Bathymetry
import random
from scipy.interpolate import lagrange
import numpy as np
import matplotlib.pyplot as plt
import os
import pdb


#inputs:
# bathymetry_data = np.array([(0,-2), (2,-1), (4, 0), (6, 3)])
bathymetry_data = np.array([(0,-2.0), (1,-1.8), (2, -1.4), (3, -1.6), (4, -1.3), (5, -1.2), (6, -1.1)])

grid_edges = np.array([0, 2, 4, 6])

if not os.path.exists('_output'):
    os.makedirs('_output')
else:
    os.system('rm -rf _output/*')
    print 'rm -rf _output/*'


# initialize the grid object:
my_grid = Grid(grid_edges, bathymetry_data)
my_grid.fill_const(-.6) 
# my_grid.eta2[:] = 0
my_grid.plot_grid() # plot level one grid
plt.savefig('_output/figure1')


print my_grid.refine()
my_grid.plot_grid() # plot level two grid
plt.savefig('_output/figure2')


print my_grid.refine(None, 2);
my_grid.plot_grid() # plot level three grid
plt.savefig('_output/figure3')


# system("open figure*")
