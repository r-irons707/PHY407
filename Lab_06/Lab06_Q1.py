# Contributors: 
# Patrick Sadnoval
# Kelvin Leong


####################################### HEADER ###########################################
# This python script holds the code, pseudo-code for the simulations for the 
# 2-dimensional molecular simulation for Q1 and all its subsections
######################################## Q1.a ############################################
# Pseudocode
# Import needed libraries
# Define timestep from problem and generate time array
# Define the inital positions for each simulation
# Define the Verlet method in a function where we only need
# inital position, and the time vectors to perform the simulation
# The Verlet function must
#   Initialize the physical constants given by the problem
#   Set the inital conditions on the position arrays
#   Compute the first iteration of the velocity components at dt/2
#   Start a for loop over len(t)-1
#   Compute the the next position instance from prev. mid velocity calc
#   Compute k vector from new position
#   Compute mid velocity from new k vector
#   Return the positions vectors from simulation
# Plot the result
# Showing x,y plot and the evolution of eah component

import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import Verlet

# Define time step 
dt = 0.01
t = np.arange(0,1000*dt,dt)

r11 = [4,4]
r21 = [5.2,4]

r12 = [4.5,4]
r22 = [5.2,4]

r13 = [2,3]
r23 = [3.5,4.4]

r1,r2 = Verlet(t,dt,r11,r21)
r3,r4 = Verlet(t,dt,r12,r22)
r5,r6 = Verlet(t,dt,r13,r23)

fig, ((a0,a1,a2),(a3,a4,a5),(a6,a7,a8)) = plt.subplots(figsize=(12,6),ncols=3,nrows=3,\
                                        gridspec_kw={'height_ratios': [3, 1, 1]})
a0.plot(r1[0],r1[1],'.',label="Particle 1",c='k')
a0.plot(r2[0],r2[1],'.',label="Particle 2",c='r')
a0.set_title("First Simulation")
a0.set_xlabel("X-Position")
a0.set_ylabel("Y-Position")
a0.legend()

a1.plot(r3[0],r3[1],'.',label="Particle 1",c='k')
a1.plot(r4[0],r4[1],'.',label="Particle 2",c='r')
a1.set_title("Second Simulation")
a1.set_xlabel("X-Position")
a1.set_ylabel("Y-Position")
a1.legend()

a2.plot(r5[0],r5[1],'.',label="Particle 1",c='k')
a2.plot(r6[0],r6[1],'.',label="Particle 2",c='r')
a2.set_title("Third Simulation")
a2.set_xlabel("X-Position")
a2.set_ylabel("Y-Position")
a2.legend()

a3.plot(t,r1[0],c='k')
a3.plot(t,r2[0],c='r')
a3.set_ylabel("X-Component")
a3.set_xlabel("Time [s]")
a4.plot(t,r3[0],c='k')
a4.plot(t,r4[0],c='r')
a4.set_ylabel("X-Component")
a4.set_xlabel("Time [s]")
a5.plot(t,r5[0],c='k')
a5.plot(t,r6[0],c='r')
a5.set_ylabel("X-Component")
a5.set_xlabel("Time [s]")

a6.plot(t,r1[1],c='k')
a6.plot(t,r2[1],c='r')
a6.set_ylabel("Y-Component")
a6.set_xlabel("Time [s]")
a7.plot(t,r3[1],c='k')
a7.plot(t,r4[1],c='r')
a7.set_ylabel("Y-Component")
a7.set_xlabel("Time [s]")
a8.plot(t,r5[1],c='k')
a8.plot(t,r6[1],c='r')
a8.set_ylabel("Y-Component")
a8.set_xlabel("Time [s]")

plt.tight_layout()
plt.savefig("Q1a.pdf")
plt.show()
