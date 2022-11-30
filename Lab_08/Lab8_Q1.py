from numpy import empty,zeros,max 
from pylab import imshow,gray,show
import matplotlib.pyplot as plt 
# Constants 
M = 100 
V = 1.0 
target = 1e-6 
omega = 0.85

phi = zeros([M+1,M+1] ,float) 
phi[20:-20,20] = V
phi[20:-20,80] = -V 
#phiprime = empty([M+1,M+1] ,float) 

# Main loop 
delta = 1.0 
while delta>target: 
    # Calculate new values of 
    old = phi.copy()
    for i in range(M+1): 
        for j in range(M+1): 
            if i==0 or i==M or j==0 or j==M: 
                phi[i,j] = phi[i,j] 
            elif 80>=i>=20 and j==20:
                phi[i,j] = phi[i,j]
            elif 80>=i>=20 and j==80:
                phi[i,j] = phi[i,j]
            else: 
                phi[i,j] = (1+omega)*(phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4 - phi[i,j]*omega
    # Calculate maximum difference from old values 
    delta = max(abs(old-phi))
    print(delta)

plt.contourf(phi) 
gray() 
show() 

