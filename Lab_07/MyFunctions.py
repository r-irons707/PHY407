# Import libraries
import numpy as np

# constants
m = 9.1094e-31
hbar = 1.0546e-34
e = 1.6022e-19
epsilon_0 = 8.854e-12 # F/m
a = 5e-11 # m (Bohr radius)
h = 0.0001*a # step size
rinf = 20*a


# potential function
def V(r):
    return -e**2/4/np.pi/epsilon_0/r

def f(v,r,E,l):
    R = v[0]
    S = v[1]
    fR = S/r**2
    fS = ((2*m*r**2/hbar**2)*(V(r)-E)+l*(l+1))*R
    return np.array([fR,fS] ,float)

# wavefunction for a particular energy
def solve_En(E,l):
    psi = 0.0
    phi = 1.0
    v = np.array([psi,phi] ,float)
    for r in np.arange(h,rinf,h):
        k1 = h*f(v,r,E,l)
        k2 = h*f(v+0.5*k1,r+0.5*h,E,l)
        k3 = h*f(v+0.5*k2,r+0.5*h,E,l)
        k4 = h*f(v+k3,r+h,E,l)
        v += (k1+2*k2+2*k3+k4)/6
    return v[0]

def solve_Rn(E,l):
    R = []
    S = 1.0
    v = np.array([0.0,S] ,float)
    for r in np.arange(h,rinf,h):
        R.append(v[0])
        k1 = h*f(v,r,E,l)
        k2 = h*f(v+0.5*k1,r+0.5*h,E,l)
        k3 = h*f(v+0.5*k2,r+0.5*h,E,l)
        k4 = h*f(v+k3,r+h,E,l)
        v += (k1+2*k2+2*k3+k4)/6
    return(R)
def En(n):
    return(-15*e/n**2,-13*e/n**2)