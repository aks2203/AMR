################################
#
# Avi Schwarzschild
# Fall 2015 - reaserch
# gird functions
# 
###############################

import numpy as np
from misc import cell_values
from BathymetryClass import Bathymetry
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import gridspec
from random import random
import pdb


class Grid:
    
    def __init__(self, grid_edges, bath_data, momentum_data = []):
        self.grid = np.array(grid_edges)      # a list of the grid edges
        self.bath_data = np.array(bath_data)
        self.bathymetry = Bathymetry(bath_data)   
        self.bathvals = self.bathymetry.bath_cell_values(grid_edges)     
        self.num_cells = len(self.grid) - 1 # number of cells in computational grid
        self.eta = np.zeros(self.num_cells)
        self.eta2 = np.zeros(self.num_cells)
        if momentum_data == []: 
            # momentum_data = np.array([(i, 2-(4*random())) for i in xrange(len(grid_edges))])
            momentum_data = np.zeros((len(grid_edges),2))
            momentum_data[:,0] = [i for i in xrange(len(grid_edges))]
        # pdb.set_trace()
        self.mu = Bathymetry(momentum_data)
        self.mu_vals = self.mu.bath_cell_values(grid_edges)
        self.mu2 = Bathymetry(momentum_data)
        self.mu2_vals = self.mu2.bath_cell_values(grid_edges)
        self.interpolants = []
        print self.grid

    def fill_cells(self, depths, momenta):
        self.eta = np.array(depths)
        self.mu = np.array(momenta)

    def fill_const(self, depth, depth2=None):
        if depth2 == None: depth2 = depth
        self.eta[0:self.num_cells/2] = depth
        self.eta[self.num_cells/2:] = depth2
        # self.mu[:] = momentum

    def fill_gaussian(self):
        def gauss(x):
            return np.exp(-(x*x))
        for i in xrange(self.num_cells):
            self.eta[i] = (round(7*gauss(i-self.num_cells/2), 2) - 3)
        self.eta2.fill(10)

    def mass_check(self):
        mass = 0; mass2 = 0;
        mu = np.sum(self.mu_vals)
        mu2 = np.sum(self.mu2_vals)
        for i in xrange(len(self.grid)-1):
            mass += (self.eta[i] - self.bathvals[i])*(self.grid[i+1]-self.grid[i])
            mass2 += (self.eta2[i] - self.eta[i])*(self.grid[i+1]-self.grid[i])

        
        return mass, mass2, mu, mu2

    
    def plot_grid(self):
        # gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        plt.figure()
        # plt.subplot(gs[0])

        pbath = [] 
        for i in self.bathvals: pbath.append(i); pbath.append(i)
        
        pmu = []
        for i in self.mu_vals: pmu.append(i); pmu.append(i)

        peta = []
        for i in self.eta: peta.append(i); peta.append(i)

        peta2 = []
        for i in self.eta2: peta2.append(i); peta2.append(i)

        for j in xrange(len(pbath)):
            if pbath[j] > peta[j]: peta[j] = pbath[j]

        pbins = []  
        for i in self.grid:
            pbins.append(i)
            if self.grid[0] < i < self.grid[-1]: pbins.append(i)
        
        plt.fill_between(pbins, min(self.bathvals) - 2, pbath, facecolor='burlywood')
        plt.fill_between(pbins, pbath, peta, facecolor='lightskyblue')
        # plt.fill_between(pbins, peta, peta2, facecolor='blue')
        # plt.plot(pbins, peta, '-ro')
        # plt.plot([x[0] for x in self.bath_data], [x[1] for x in self.bath_data], '-g')

        # if len(self.grid) < 30: 
        for i in self.grid: 
            plt.plot([i, i], [min(self.bathvals) - 2, max(self.eta)+5], color='black', ls='--')


        tan_patch = mpatches.Patch(color='burlywood', label='Bathymetry')
        blue_patch = mpatches.Patch(color='lightskyblue', label='water')
        green_patch = mpatches.Patch(color='green', label='Bathymetry interpolants')
        #plt.legend(handles=[tan_patch, blue_patch, green_patch])

        # plt.ylim(min(self.bathvals) - 2, max(self.eta)+10)
        plt.xlabel('Horizontal dimenssion (grey dotted lines on grid boundries)')

        midpoints = []
        for i in xrange(len(self.grid)-1):
            midpoints.append(self.grid[i] + (self.grid[i+1]-self.grid[i])/2.0)

        # plt.subplot(gs[1])
        # plt.plot(midpoints, self.mu_vals,'ko')
        # if len(self.grid) < 30: 
        #     for i in self.grid: plt.axvline(i, min(self.bathvals) - 2, max(self.eta)+5, color='grey', ls='--')

    def grid_cell_values(self, grid, interpolants, integrals):
        return cell_values(grid, interpolants, integrals)    

    def get_new_grid(self, factor):
        self.new_grid = []
        for i in xrange(len(self.grid)):
            if self.grid[i] == self.grid[-1]: self.new_grid.append(self.grid[i])
            else:
                for j in xrange(factor):
                    self.new_grid.append(self.grid[i]+j*(self.grid[i+1]-self.grid[i])/(factor+0.0))
        self.new_grid = np.array(self.new_grid)

    def refine_surface(self, eta, factor):
        slope = 0
        interpolants = []
        for i in xrange(len(eta)):
            if i == 0 : 
                slope = (eta[i+1] - eta[i])/(self.grid[i+1] - self.grid[i]+0.0)
            elif i == len(eta) - 1: 
                slope = (eta[i] - eta[i-1])/(self.grid[i] - self.grid[i-1])
            else:
                s1 = (eta[i] - eta[i-1])/(self.grid[i] - self.grid[i-1])
                s2 = (eta[i+1] - eta[i])/(self.grid[i+1] - self.grid[i])
                if s1*s2 <= 0: slope = 0
                else:
                    slope = ((abs(s1) <= abs(s2))*s1 + (abs(s1) > abs(s2))*s2)
            p = np.poly1d([slope, eta[i] - slope*(self.grid[i] + (self.grid[i+1]-self.grid[i])/2.0)])
            interpolants.append((self.grid[i], self.grid[i+1], p))

        integrals = []
        for i in interpolants:
            integrals.append(np.polyint(i[2]))

        new_eta = self.grid_cell_values(self.new_grid, interpolants, integrals)

        return new_eta

    def refine_mu(self, mu, factor):
        slope = 0
        interpolants = []
        for i in xrange(len(mu)):
            if i == 0 : 
                slope = (mu[i+1] - mu[i])/(self.grid[i+1] - self.grid[i]+0.0)
            elif i == len(mu) - 1: 
                slope = (mu[i] - mu[i-1])/(self.grid[i] - self.grid[i-1])
            else:
                s1 = (mu[i] - mu[i-1])/(self.grid[i] - self.grid[i-1])
                s2 = (mu[i+1] - mu[i])/(self.grid[i+1] - self.grid[i])
                if s1*s2 <= 0: slope = 0
                else:
                    slope = ((abs(s1) <= abs(s2))*s1 + (abs(s1) > abs(s2))*s2)
            p = np.poly1d([slope, mu[i] - slope*(self.grid[i] + (self.grid[i+1]-self.grid[i])/2.0)])
            interpolants.append((self.grid[i], self.grid[i+1], p))

        integrals = []
        for i in interpolants:
            integrals.append(np.polyint(i[2]))

        new_eta = self.grid_cell_values(self.new_grid, interpolants, integrals)

        return [i/(factor+0.0) for i in new_eta]

    def refine(self, grid = None, factor=2):
        # Mass check
        mass_a1, mass_b1, mu_a1, mu_b1 = self.mass_check()

        # Get new grid
        if grid != None: self.new_grid = grid
        else: self.get_new_grid(factor)
    
        # Refine bathymetry
        self.bathvals = self.bathymetry.bath_cell_values(self.new_grid)

        # Refine surface one
        self.eta = self.refine_surface(self.eta, factor)
        self.mu_vals = self.refine_mu(self.mu_vals, factor)

        # Refine surface two
        self.eta2 = self.refine_surface(self.eta2, factor)
        self.mu2_vals = self.refine_mu(self.mu2_vals, factor)

        self.grid = self.new_grid

        # Mass check
        mass_a2, mass_b2, mu_a2, mu_b2 = self.mass_check()

        return np.allclose(mass_a1, mass_a2), np.allclose(mu_a1, mu_a2), np.allclose(mass_b1, mass_b2), np.allclose(mu_b1, mu_b2)

     


