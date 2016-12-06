################################
#
# Avi Schwarzschild
# Fall 2016 - reaserch
# grid functions
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

dry_tolerance = 1e-3

class Grid:
    
    def __init__(self, grid_edges, bath_data, momentum_data = []):
        self.grid = np.array(grid_edges)      # a list of the grid edges
        self.bath_data = np.array(bath_data)
        self.bathymetry = Bathymetry(bath_data)   
        self.bathvals = self.bathymetry.bath_cell_values(grid_edges)     
        self.num_cells = len(self.grid) - 1 # number of cells in computational grid
        self.eta = np.zeros(self.num_cells)
        self.eta2 = np.zeros(self.num_cells)
        self.eta_init = []
        self.h = np.zeros(self.num_cells)
        self.h2 = np.zeros(self.num_cells)
        if momentum_data == []: 
            # momentum_data = np.array([(i, 2-(4*random())) for i in xrange(len(grid_edges))])
            momentum_data = np.zeros((len(grid_edges),2))
            momentum_data[:,0] = [i for i in grid_edges]
        self.mu = Bathymetry(momentum_data)
        self.mu_vals = self.mu.bath_cell_values(grid_edges)
        self.mu2 = Bathymetry(momentum_data)
        self.mu2_vals = self.mu2.bath_cell_values(grid_edges)
        self.interpolants = []

    def fill_cells(self, surface_heights, momenta):
        self.eta = np.array(surface_heights)
        self.mu = np.array(momenta)
        self.h = self.eta - self.bathvals

    def fill_const(self, surface_height, surface_height2=None):
        if surface_height2 == None: surface_height2 = surface_height
        self.eta[0:self.num_cells/2] = surface_height
        self.eta[self.num_cells/2:] = surface_height2
        self.eta2[:] = 0
        self.eta_init = [surface_height, 0]
        self.h = (self.eta - self.bathvals)*((self.eta-self.bathvals) > 0) + 0*((self.eta-self.bathvals) <= 0)
        h2 = self.eta2 - self.h - self.bathvals
        self.h2 = h2*(h2>0) + 0 *(h2<=0)
 

    def fill_gaussian(self):
        self.fill_const(self.eta_init[0], self.eta_init[1])
        def gauss(x):
            return np.exp(-(x*x))
        for i in xrange(self.num_cells):
            self.eta[i] = (round(7*gauss(i-self.num_cells/2), 2) - 3)
        self.eta2.fill(10)
        self.h = self.eta - self.eta2
        self.h2 = self.eta2 - self.bathvals

    def mass_check(self):
        mass = 0; mass2 = 0;
        mu = np.sum(self.mu_vals)
        mu2 = np.sum(self.mu2_vals)
        # for i in xrange(len(self.grid)-1):
        #     mass += (self.h[i])*(self.grid[i+1]-self.grid[i])
        #     mass2 += (self.h2[i])*(self.grid[i+1]-self.grid[i])
        mass = np.sum(self.h)*self.grid[1]-self.grid[0]
        mass2 = np.sum(self.h2)*self.grid[1]-self.grid[0]


        return mass, mass2, mu, mu2

    
    def plot_grid(self):
        # gs = gridspec.GridSpec(2, 1, height_ratios=[4,1])
        plt.figure()
        # plt.subplot(gs[0])

        pbath = np.repeat(self.bathvals, 2)
        peta = pbath + np.repeat(self.h, 2)
        peta2 = peta + np.repeat(self.h2, 2)

        for j in xrange(len(pbath)):
            if pbath[j] > peta[j]: peta[j] = pbath[j]

        pbins = []  
        for i in self.grid:
            pbins.append(i)
            if self.grid[0] < i < self.grid[-1]: pbins.append(i)
        
        plt.fill_between(pbins, min(self.bathvals) - 2, pbath, facecolor='burlywood')
        plt.fill_between(pbins, pbath, peta, facecolor='blue')
        plt.fill_between(pbins, peta, peta2, facecolor='lightskyblue')
        # plt.plot(pbins, peta, '-ro')
        plt.plot([x[0] for x in self.bath_data], [x[1] for x in self.bath_data], '-g', linewidth=2)

        if len(self.grid) < 30: 
            for i in self.grid: 
                plt.plot([i, i], [min(self.bathvals) - 2, max(self.eta)+5], color='black', ls='--')


        tan_patch = mpatches.Patch(color='burlywood', label='Bathymetry')
        blue_patch = mpatches.Patch(color='lightskyblue', label='water')
        green_patch = mpatches.Patch(color='green', label='Bathymetry interpolants')
        # plt.legend(handles=[tan_patch, blue_patch, green_patch])

        plt.ylim(min(self.bathvals) - 2, max(self.eta)+4)
        plt.xlabel('Horizontal dimenssion (grey dotted lines on grid boundries)')

        midpoints = []
        for i in xrange(len(self.grid)-1):
            midpoints.append(self.grid[i] + (self.grid[i+1]-self.grid[i])/2.0)

        # plt.subplot(gs[1])
        # plt.plot(midpoints, self.mu_vals,'ko')
        if len(self.grid) < 30: 
            for i in self.grid: plt.axvline(i, min(self.bathvals) - 2, max(self.eta)+5, color='grey', ls='--')

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

    def find_etas(self):
        #This would take into account dry cells
        eta = (self.h > dry_tolerance)*(self.h + self.bathvals) + (self.h <= dry_tolerance)*self.eta_init[0]
        eta2 = (self.h2 > dry_tolerance)*(self.h2 + self.h + self.bathvals) + (self.h2 <= dry_tolerance)*self.eta_init[1]
        return eta, eta2


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
        
        self.eta, self.eta2 = self.find_etas()

        # Refine bathymetry
        self.bathvals = self.bathymetry.bath_cell_values(self.new_grid)

        # Refine surface one
        eta = self.refine_surface(self.eta, factor)
        self.mu_vals = self.refine_mu(self.mu_vals, factor)


        # Refine surface two
        eta2 = self.refine_surface(self.eta2, factor)
        self.mu2_vals = self.refine_mu(self.mu2_vals, factor)

        self.grid = self.new_grid

        #eta to h
        # self.h = self.eta - self.bathvals
        # self.h2 = self.eta2 - self.eta
        self.h = (eta - self.bathvals)*((eta-self.bathvals) > dry_tolerance) + 0*((eta-self.bathvals) <= dry_tolerance)
        h2 = eta2 - self.h - self.bathvals
        self.h2 = h2*(h2>dry_tolerance) + dry_tolerance *(h2<=dry_tolerance)

        # Mass check
        mass_a2, mass_b2, mu_a2, mu_b2 = self.mass_check()

        # print self.bathvals

        return np.allclose(mass_a1+mass_b1, mass_a2+mass_b2), np.allclose(mass_b1, mass_b2)



