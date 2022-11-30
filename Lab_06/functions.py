import numpy as np

# Define acceleratino components from Q1a
def ax(x,y,a,b,sigma,epsilon,m):
    '''Returns the x-component of the acceleration
    for a particle subject to a Lennard-Jones
    potential. Where a and b are the x,y coordinates
    of the opposite particle'''
    term1 = (48*epsilon*sigma**12)/((x-a)**2 + (y-b)**2)**13
    term2 = (24*epsilon*sigma**6)/((x-a)**2 + (y-b)**2)**7
    return ((x-a)/m)*(term1 - term2)


def ay(x,y,a,b,sigma,epsilon,m):
    '''Returns the y-component of the acceleration
    for a particle subject to a Lennard-Jones
    potential.Where a and b are the x,y coordinates
    of the opposite particle'''
    term1 = (48*epsilon*sigma**12)/((x-a)**2 + (y-b)**2)**13
    term2 = (24*epsilon*sigma**6)/((x-a)**2 + (y-b)**2)**7
    return ((y-b)/m)*(term1 - term2)

def Verlet(t,dt,x,y):
    '''Performs the Verlet method to simulate
    the interaction of 2 particles under Lennard
    Jones potential. r01 and r02 are 1D arrays
    of length 2 where first entry is x0 and second
    entry is y0
    '''
    # Define constants for simulation
    epsilon = 1
    sigma = 1
    m = 1
    
    # Initialize arrays and set initial conditions
    vx1 = [0]
    vy1 = [0]

    vx2 = [0]
    vy2 = [0]

    x1 = [r01[0]]
    y1 = [r01[1]]

    x2 = [r02[0]]
    y2 = [r02[1]]

    # Compute first term
    v0x1 = vx1[0] + dt/2*ax(x1[0],y1[0],x2[0],y2[0],sigma,epsilon,m)
    v0y1 = vy1[0] + dt/2*ay(x1[0],y1[0],x2[0],y2[0],sigma,epsilon,m)

    v0x2 = vx2[0] + dt/2*ax(x2[0],y2[0],x1[0],y1[0],sigma,epsilon,m) 
    v0y2 = vy2[0] + dt/2*ay(x2[0],y2[0],x1[0],y1[0],sigma,epsilon,m)

    for i in range(len(t)-1):
        x1.append(x1[i] + dt*v0x1)
        y1.append(y1[i] + dt*v0y1)

        x2.append(x2[i] + dt*v0x2)
        y2.append(y2[i] + dt*v0y2)

        k1x = dt * ax(x1[i+1],y1[i+1],x2[i+1],y2[i+1],sigma,epsilon,m)
        k1y = dt * ay(x1[i+1],y1[i+1],x2[i+1],y2[i+1],sigma,epsilon,m)

        k2x = dt * ax(x2[i+1],y2[i+1],x1[i+1],y1[i+1],sigma,epsilon,m)
        k2y = dt * ay(x2[i+1],y2[i+1],x1[i+1],y1[i+1],sigma,epsilon,m)

        v0x1 += k1x
        v0y1 += k1y

        v0x2 += k2x
        v0y2 += k2y

    r1 = np.array([x1,y1])
    r2 = np.array([x2,y2])
    return(r1,r2)

def distance(x,y ,a,b):
    '''Compute distance between two particles given position'''
    return (x-a)**2 + (y-b)**2

def Verlet_N(t,dt,N,initial_positions):
    '''Performs the Verlet method to simulate
    the interaction of 2 particles under Lennard
    Jones potential. r01 and r02 are 1D arrays
    of length 2 where first entry is x0 and second
    entry is y0 '''
    
    # Define constants for simulation
    epsilon = 1
    sigma = 1
    m = 1
    
    # Initialize arrays and set initial conditions
    V = np.zeros(shape=(N,2)) # initial velocities set to 0 for each particle
    vx1 = [0]
    vy1 = [0]

    vx2 = [0]
    vy2 = [0]
    
    # compute acceleration (analogous to force here) of each particle wrt all others for 
    for range(N):
        

    x1 = [r01[0]]
    y1 = [r01[1]]

    x2 = [r02[0]]
    y2 = [r02[1]]

    # Compute first term
    v0x1 = vx1[0] + dt/2*ax(x1[0],y1[0],x2[0],y2[0],sigma,epsilon,m)
    v0y1 = vy1[0] + dt/2*ay(x1[0],y1[0],x2[0],y2[0],sigma,epsilon,m)

    v0x2 = vx2[0] + dt/2*ax(x2[0],y2[0],x1[0],y1[0],sigma,epsilon,m) 
    v0y2 = vy2[0] + dt/2*ay(x2[0],y2[0],x1[0],y1[0],sigma,epsilon,m)

    for i in range(len(t)-1):
        x1.append(x1[i] + dt*v0x1)
        y1.append(y1[i] + dt*v0y1)

        x2.append(x2[i] + dt*v0x2)
        y2.append(y2[i] + dt*v0y2)

        k1x = dt * ax(x1[i+1],y1[i+1],x2[i+1],y2[i+1],sigma,epsilon,m)
        k1y = dt * ay(x1[i+1],y1[i+1],x2[i+1],y2[i+1],sigma,epsilon,m)

        k2x = dt * ax(x2[i+1],y2[i+1],x1[i+1],y1[i+1],sigma,epsilon,m)
        k2y = dt * ay(x2[i+1],y2[i+1],x1[i+1],y1[i+1],sigma,epsilon,m)

        v0x1 += k1x
        v0y1 += k1y

        v0x2 += k2x
        v0y2 += k2y

    r1 = np.array([x1,y1])
    r2 = np.array([x2,y2])
    return(r1,r2)