#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Contributors: 
# Patrick Sadnoval
# Kelvin Leong


################################## HEADER ####################################
# This python script holds the code, pseudo-code for the simulations for the 
# 2-dimensional molecular simulation for Q2 and all its subsections
##############################################################################
# Pseudocode:
# Import the necessary libraries
# Set N number of particles simulating, x-separation and y-separation between
# each particles (evenly spaced)
# Set the initial positions of all particles as an array, COM located at (0,0)
# Set the initial velocity as zero for all particles as an array
# Define timestep dt=0.01, end time T=1000*dt, and generate time array
# Modify and define the Verlet method for multibody in a function 
# The Verlet function:
#   Initialize the physical constants given by the problem
#   Initialize x,y, full-step vx,vy, half-step vx,vy as 2D numpy array, with 
#   N particles on row axis and time t on column axis
#   Initialize the kinetic and potential energy as a numpy array wrt time
#   Compute the first iteration of the half-step velocity components at dt/2 (Eq.8)
#   for each particle, as well as the energy at t = 0
#   Start a time for-loop from 1 to len(t):
#       Compute the the positions from prev. half-step velocity for all particles
#       first (Eq.9)
#       Calculate acceleration and potential energy of particle i due to all 
#       other j particles
#       Compute k vector from this acceleration (Eq.10)
#       Compute full-step velocity from this k vector (Eq.11)
#       Compute half-step velocity from this k vector (Eq.12)
#       Calculate kinetic energy using the current full-step velocity
#   Sum kinetic and potential
#   Return all relavant arrays from used in Verlet
# Plot the simulation on xy-plane
# Plot the energy as a function of time

import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import Verlet_multibody

# Number of particles simulating
N = 16

# Set initial conditions for the N particles
Lx = 4.0
Ly = 4.0

dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)

x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)

xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)

x_initial = xx_grid.flatten() - Lx/2
y_initial = yy_grid.flatten() - Ly/2


# Assign the initial velocities of all N particles to zero
vx_initial = np.zeros(N)
vy_initial = np.zeros(N)

# Define time step 
dt = 0.01
T = 1000*dt
t = np.arange(0,T,dt)

# Group-up initial conditions
r_initial = np.array([x_initial, y_initial])
v_initial = np.array([vx_initial, vy_initial])


x_all_t, y_all_t, vx_all_t, vy_all_t, vxhalf_all_t, vyhalf_all_t, V_t, K_t = Verlet_multibody(t,dt,r_initial,v_initial,N)

# Calculate total energy of the system
E_t = V_t + K_t

#%%
# Plotting
plt.figure(figsize=(10,10))
for i in range(N):
    plt.plot(x_all_t[i,:],y_all_t[i,:],'.')
plt.title("N=16 particles Simulation",fontsize=16)
plt.xlabel("X-Position",fontsize=16)
plt.ylabel("Y-Position",fontsize=16)
plt.grid('on')
plt.tight_layout()
plt.savefig("Q2a.png")

# Plot energy, median of energy and the 10% interval from median
avg_E = np.median(E_t)
plt.figure(figsize=(8,8))
plt.plot(t,E_t, 'r.')
plt.plot(t, np.ones(len(t))*avg_E, 'k--', label='Median energy', markersize=5)
plt.plot(t, np.ones(len(t))*avg_E + 0.1*avg_E, 'm--', label='+10% from median', markersize=5)
plt.plot(t, np.ones(len(t))*avg_E - 0.1*avg_E, 'b--', label='-10% from median', markersize=5)
plt.xlabel('Time [s]',fontsize=16)
plt.ylabel('Energy',fontsize=16)
plt.title('N=16 particles system energy', fontsize=16)
plt.grid('on')
plt.legend()
plt.tight_layout()
plt.savefig("Q2b.png")
plt.show() 


    
                


