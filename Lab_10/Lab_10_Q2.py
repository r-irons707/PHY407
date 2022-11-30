from random import random
import numpy as np

# function we are integrating
def f(r):
    if np.sum(r**2) <= 1:
        return(1)
    else:
        return(0)

a = -1 # define bounds
b = 1
N = 1e6 # sample size
dim = 10 # number of dimensions
tot_points = 0 # sum variable

# define volume of hypercube enclosing the hypersphere
L = 2**dim

for i in range(int(N)):
    # Generate dof random numbers within region (-1,1)
    xi = np.array([(b-a)*random() + a for _ in range(dim)]) # applying a shift on the sampling region
    # Compute function for random numbers
    y = f(xi)
    tot_points += y # update the total number of points in the hypersphere

volume = L/N*tot_points
print(f"Volume of sphere in {dim} dimensions:",volume)

