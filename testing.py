################################
#
# Avi Schwarzschild
# Fall 2015 - reaserch
# testing.py for testing...
# 
###############################

from Grid import Grid
from BathymetryClass import Bathymetry
import random
from scipy.interpolate import lagrange
import numpy as np
import matplotlib.pyplot as plt


eta = [] # a list of water hights with respect to z=0
b =  [] # a list of bathymetry values with repects to z=0
# h = b+eta

#################################
#
# model a gausian hump of water
# on the line from 0-20
# initial bath data for every integer
# point. First level grid is uniform
# with cell width of 2.
#
#################################

#inputs:
bathymetry_data = [(i, -10-int(5*random.random() + 1)) for i in xrange(21)]
grid_edges = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

#initialize bathymetry object:
my_bath = Bathymetry(bathymetry_data)

#initialize the grid object:
my_grid = Grid(grid_edges, bathymetry_data)
my_grid.fill_gaussian()
my_grid.plot_grid()

my_grid.refine()
my_grid.plot_grid()