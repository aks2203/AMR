################################
#
# Avi Schwarzschild
# Fall 2015 - reaserch - v7
# main2.py
# 
###############################


from Grid import Grid
from BathymetryClass import Bathymetry
import random
from scipy.interpolate import lagrange
import numpy as np
import matplotlib.pyplot as plt
from os import system
import pdb

#################################
#
# model a gausian hump of water with an island
# on the line from 0-20
# initial bath data for every integer
# point. First level grid is uniform
# with cell width of 2.
#
#################################

#inputs:
#bathymetry_data = [(0,-20), (1,-18), (2, -16), (3, -16), (4, -12), (5, -14), (6, -8)]
bathymetry_data = np.array([(0,-20), (1,-18), (2, -16), (3, -16), (4, -12), (5, -14), (6, -8)])

grid_edges = np.array([0, 1, 2, 3, 4, 5, 6])

# initialize the grid object:
my_grid = Grid(grid_edges, bathymetry_data)
my_grid.fill_const(-15.5) 
my_grid.eta2[:] = 0
my_grid.plot_grid() # plot level one grid
plt.savefig('figure1')
print type(my_grid.grid)

print my_grid.refine(np.array([0, 2, 4, 6]))
my_grid.plot_grid() # plot level two grid
plt.savefig('figure2')
print type(my_grid.grid)

print my_grid.refine(None, 2);
my_grid.plot_grid() # plot level three grid
plt.savefig('figure3')

print type(my_grid.grid)
print my_grid.grid

# system("open figure*")
