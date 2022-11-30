import numpy as np
import matplotlib.pyplot as plt

# part a
# simulate a particle choosing 1 of 4 directions on a lattice of shape LxL, let L=50, and choosing a separate direction if
# the particle chooses to move outside the lattice
# defining conditions
step = 5000 # number of time steps, increase to 5000 later
time = np.arange(1,step+1,1)
L = 101 # originally 101
position_x = np.zeros(step,dtype=int)
position_x[0] = (L-1)//2
position_y = np.zeros(step,dtype=int)
position_y[0] = (L-1)//2

for i in range(step-1):
    # choose a random integer [1,4] from a uniform distribution (equal probability for each) to determine step direction
    choice = np.random.randint(1, high=5, size=None, dtype=int)
    if position_x[i] >= L:
        choice = 3
    if position_y[i] >= L:
        choice = 2
    if position_x[i] <= 0:
        choice = 4
    if position_y[i] <= 0:
        choice = 1
    #########################
    if choice == 1: # symbols 'up' direction
        position_x[i+1] = position_x[i]
        position_y[i+1] = position_y[i] + 1
    if choice == 2: # symbols 'down' direction
        position_x[i+1] = position_x[i]
        position_y[i+1] = position_y[i] - 1
    if choice == 3: # symbols 'left' direction
        position_x[i+1] = position_x[i] - 1
        position_y[i+1] = position_y[i]
    if choice == 4: # symbols 'right' direction
        position_x[i+1] = position_x[i] + 1
        position_y[i+1] = position_y[i]

fig, (ax0,ax1,ax2) = plt.subplots(figsize=(14,4),ncols=3,nrows=1)

ax0.plot(position_x,position_y)
ax0.set_title('x vs y position of particle in a 100x100 box')
ax0.set_xlabel('x position')
ax0.set_ylabel('y position')
ax1.plot(time,position_x)
ax1.set_title('position of particle in x direction over time')
ax2.plot(time,position_y)
ax2.set_title('position of particle in y direction over time')