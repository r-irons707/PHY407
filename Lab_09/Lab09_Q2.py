########################## 
# Contributors           # 
# Patrick Sandoval       #
# Kelvin Leong           #
##########################

################################# HEADER ##########################################
# This python script holds the code and plots for Q2 of Lab09 where we use        #  
# the functions defined in dcst.py and dcst_for_q2.py to perform the fourier      #
# decompositon of the functions                                                   #
#                                                                                 #
###################################################################################

####################################### Q2a #######################################
# Import needed libraries
import numpy as np
import matplotlib.pyplot as plt

# Import functions from scripts
from dcst_for_q2 import *

def dXXt2(f,xs,ys):
    """ Takes DXT along x, then DXT along y (X = C/S)
    IN: f, the input 2D numpy array
    OUT: b, the 2D transformed array
    xs and ys are the symmetry of the input
    function along that specific axis """
    M = f.shape[0] # Number of rows
    N = f.shape[1] # Number of columns
    a = np.zeros((M, N)) # Intermediate array
    b = np.zeros((M, N)) # Final array

    # Take transform along x
    for j in range(N):
        # Check if function is even along this direction
        if xs=="E":
            # Perform DCT
            a[:,j] = dct(f[:,j])
        # Check if function is odd
        elif xs=="O":
            # Perform DST
            a[:,j] = dst(f[:,j])
        else:
            raise "Function is not symmetric nor anti-symmetric along x"

    # Take transform along y
    for i in range(M):
        # Check if function is even along this direction
        if ys=="E":
            # Perform DCT
            b[i,:] = dct(a[i,:])
        # Check if function is odd along this direction
        elif ys=='O':
            # Perform DST
            b[i,:] = dst(a[i,:])
        else:
            raise "Function is not symmetric nor anti-symmetric along y"
    return b

def idXXt2(b,xs,ys):
    """ Takes iDXT along y, then iDXT along x (X = C/S)
    IN: b, the input 2D numpy array
    OUT: f, the 2D inverse-transformed array
    xs and ys are the symmetry of the input
    function along that specific axis """
    M = b.shape[0] # Number of rows
    N = b.shape[1] # Number of columns
    a = np.zeros((M, N)) # Intermediate array
    f = np.zeros((M, N)) # Final array
    # Take inverse transform along y
    for i in range(M):
        if ys=="E":
            a[i:] = idct(b[i,:])
        elif ys=="O":
            a[i:] = idst(b[i,:])
        else:
            raise "Function is not symmetric or anti-symmetric along x"
        # iDXT b[i,:] and set as a[i,:]

    # Take inverse transform along x
    for j in range(N):
        if xs=="E":
            f[:,j] = idct(a[:,j])
        elif xs=='O':
            f[:,j] = idst(a[:,j])
        else:
            raise "Function is not symmetric or anti-symmetric along y"
        # iDXT a[:,j] and set as f[:,j]
    return f

# To test this functions we will generate a 2D gaussian 
x = np.linspace(-5,5,100)
y = np.linspace(-5,5,100)
x, y = np.meshgrid(x, y)

def gauss2D(x=0, y=0, mx=0, my=0, sx=1, sy=1):
    return 1. / (2. * np.pi * sx * sy) * np.exp(-((x - mx)**2. / (2. * sx**2.) + (y - my)**2. / (2. * sy**2.)))
z = gauss2D(x,y)

# Take 2D FFT and iFFT 
f = dXXt2(z,"E","E")
zprime = idXXt2(f,"E","E")

# Plot result
fig, (a0,a1,a2) = plt.subplots(figsize=(12,4),ncols=3)
a0.contourf(x,y,z,cmap='inferno')
a0.set_title("Original Signal")
a0.set_xlabel("X")
a0.set_ylabel("Y")
a1.imshow(f,origin='lower',cmap='inferno')
a1.set_xlim([0,10])
a1.set_ylim([0,10])
a1.set_title("2D Fourier Transform")
a1.set_xlabel("Spectral X")
a1.set_ylabel("Spectral Y")
a2.contourf(x,y,zprime,cmap='inferno')
a2.set_title("Recovered Signal")
a2.set_xlabel("X")
a2.set_ylabel("Y")
plt.tight_layout()
#plt.savefig("Lab9/Q2aPlot.pdf")
plt.show()

# Print out result of test
print("#"*15+"Test"+"#"*15)
print("Residuals:",abs(z-zprime))

####################################### Q2b #######################################
# Define physical constant and system parameters
Lx, Ly, Jo, m, n, c = 1, 1, 1, 1, 1, 1
P = 32
dx = Lx/P
dy = Ly/P
T = 20
omega = 3.75
dt = 0.01

# Define solution arrays given initial condtions
Ez = np.zeros((P+1,P+1))
Hx = np.zeros((P+1,P+1))
Hy = np.zeros((P+1,P+1))

# To satisfy boundary condiions the solutions must follow
# Ez odd in (x,y)
# Hx odd in x, even in y
# Hy even in x, odd in x
# Jz odd in (x,y)

# Creat discrete arrays
t = np.arange(0,T,dt)
x_domain = np.arange(0,Lx+dx,dx)
y_domain = np.arange(0,Ly+dx,dy)
x, y = np.meshgrid(x_domain,y_domain)
kx = np.arange(dx/2 +1)*2*np.pi/Lx
ky = np.arange(dy/2 +1)*2*np.pi/Ly

def J(x,y,t):
    '''Returns a N X N X M matrix for the current density
    defined by the problem'''
    Jz = []
    for tval in t:
        tmp = Jo*np.sin(m*np.pi*x/Lx)*np.sin(n*np.pi*y/Ly)*np.sin(omega*tval)
        Jz.append(tmp)
    return np.stack(Jz,axis=0)

# Define current density and its time evolution
# as a cube where every layer is J at some time t + i*dt      
Jz = J(x,y,t)

# Take fourier transform of Jz at all its consequent times
Jzhat_tmp = []
for mat in Jz:
    Jzhat_tmp.append(dXXt2(mat,xs="O",ys="O"))
Jzhat = np.stack(Jzhat_tmp,axis=0)

# Take fourier transform of Hx, Hy and Ez at t=0
Hxhat = dXXt2(Hx,xs='O',ys='E')
Hyhat = dXXt2(Hy,xs='E',ys='O')
Ezhat = dXXt2(Ez,xs='O',ys='O')

# Write arrays as cubes
EzhatS = np.zeros((len(t),P+1,P+1))
HxhatS = np.zeros((len(t),P+1,P+1))
HyhatS = np.zeros((len(t),P+1,P+1))

EzhatS[0] = Ezhat
HxhatS[0] = Hxhat
HyhatS[0] = Hyhat

# Evolve Fourier coefficients using CN-Method
# Define constants
Dx = np.pi*c*dt/2/Lx
Dy = np.pi*c*dt/2/Ly

Eztmp = []
Hxtmp = []
Hytmp = []
for i in range(len(t)):
    Ezt2 = []
    Hxt2 = []
    Hyt2 = []
    for k in range(len(y_domain)):
        Ezt = []
        Hxt = []
        Hyt = []
        for j in range(len(x_domain)):
            Ei = ((1-x_domain[j]**2*Dx**2 - y_domain[k]**2*Dy**2)*EzhatS[i,k,j] \
                + 2*y_domain[k]*Dy*HxhatS[i,k,j] + 2*x_domain[j]*Dx*HyhatS[i,k,j] + dt*Jzhat[i,k,j])\
                    /(1+x_domain[j]**2*Dx**2 + y_domain[k]*Dy**2)
            Hxi = HxhatS[i,k,j] - y_domain[k]*Dy*(Ei + EzhatS[i,k,j])
            Hyi = HyhatS[i,k,j] - x_domain[j]*Dx*(Ei + EzhatS[i,k,j])
            Ezt.append(Ei)
            Hxt.append(Hxi)
            Hyt.append(Hyi)
        Ezt2.append(Ezt)
        Hxt2.append(Hxt)
        Hyt2.append(Hyt)
    Eztmp.append(Ezt2)
    Hxtmp.append(Hxt2)
    Hytmp.append(Hyt2)
EzHatStack = np.stack(Eztmp,axis=0)
HxHatStack = np.stack(Hxtmp,axis=0)
HyHatStack = np.stack(Hytmp,axis=0)

# Now we take the inveser fourier to recover solution
EzSoltmp = []
HxSoltmp = []
HySoltmp = []
for mat in EzHatStack:
    sol = idXXt2(mat,xs='O',ys='O')
    EzSoltmp.append(sol)
for mat in HxHatStack:
    sol = idXXt2(mat,xs="O",ys="E")
    HxSoltmp.append(sol)
for mat in HyHatStack:
    sol = idXXt2(mat,xs='E',ys='O')
    HySoltmp.append(sol)
EzSol = np.stack(EzSoltmp,axis=0)
HxSol = np.stack(HxSoltmp,axis=0)
HySol = np.stack(HySoltmp,axis=0)

# Now I need to find what index along x or y corresponds to the
# value the lab wants me to plot 
x1_tar = int(0.5/dx)
y1_tar = int(0.0/dx)

x2_tar = int(0.0/dx)
y2_tar = int(0.5/dx)

x3_tar = int(0.5/dx)
y3_tar = int(0.5/dx)

# Collect solutions
HxPlotSol = HxSol[:,y2_tar,x2_tar]
HyPlotSol = HySol[:,y1_tar,x1_tar]
EzPlotSol = EzSol[:,y3_tar,x3_tar]

fig, (a0,a1,a2) = plt.subplots(figsize=(12,3),ncols=3)
a0.plot(t,HxPlotSol,c="b")
a0.set_xlabel("Time [s]")
a0.set_ylabel(r"H$_x$(x=0.5,y=0)")
a1.plot(t,HyPlotSol,c="g")
a1.set_xlabel("Time [s]")
a1.set_ylabel(r"H$_y$(x=0,y=0.5)")
a2.plot(t,EzPlotSol,c="r")
a2.set_xlabel("Time [s]")
a2.set_ylabel(r"E$_z$(x=0.5,y=0.5)")
plt.tight_layout()
plt.savefig("Lab9/Q2bPlots.pdf")
plt.show()


# The following code animates the EM cavity its really cool
from pylab import clf, plot, xlim, ylim, show, pause, draw
for sol in EzSol:
    clf()
    plt.contourf(x_domain,y_domain,sol,cmap='inferno',vmin=np.min(EzSol),vmax=np.max(EzSol))
    plt.colorbar()
    xlim([0,Lx])
    ylim([0,Ly])
    plt.xlabel("X",fontsize=16)
    plt.ylabel("Y",fontsize=16)
    plt.title(r"Oscillating E$_z$ Field",fontsize=18)
    plt.tight_layout()
    draw()
    pause(0.001)
