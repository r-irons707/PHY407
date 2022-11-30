"""
This program simulates Brownian motion in the presence of walls
Note that the physical behaviour would be to stick to walls,
which is the purpose of Q1a.
Author: Nico Grisouard, University of Toronto
"""

import numpy as np
import matplotlib.pyplot as plt

def nextmove(x, y):
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    direction = np.random.randint(0, high=4, size=None, dtype=int)

    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y


Lp = 101  # size of domain
Nt = 50  # number of time steps
time = np.arange(1,Nt+1,1)
# arrays to record the trajectory of the particle
pos_x = np.zeros(Nt,dtype=int)
pos_x[0] = (Lp-1)//2
pos_y = np.zeros(Nt,dtype=int)
pos_y[0] = (Lp-1)//2

centre_point = (Lp-1)//2  # middle point of domain
xp = centre_point
yp = centre_point

for i in range(Nt-1):
    xpp, ypp = nextmove(xp, yp)
    print(xpp)
    pos_x[i+1] = xpp
    pos_y[i+1] = ypp

print(pos_x)

fig, (ax0,ax1,ax2) = plt.subplots(figsize=(14,4),ncols=3,nrows=1)

ax0.plot(pos_x,pos_y)
ax0.set_title('x vs y position of particle in a 50x50 box')
ax1.plot(time,pos_x)
ax1.set_title('position of particle in x direction over time')
ax2.plot(time,pos_y)
ax2.set_title('position of particle in y direction over time')
plt.show()