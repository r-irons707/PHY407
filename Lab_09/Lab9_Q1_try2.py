# Imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc

# Constants
L = 1e-8
m = 9.109e-31
sigma = L / 25
k = 500 / L
tau = 1e-18  # time step
x0 = L / 5
N = 3000  # Iteration of time steps
T = N * tau
P = 1024  # Number of segments
V_0 = 6e-17
x1 = L / 4
interval = L / P

#hbar = 6.626e-34

def trapezoidal_method(a, b, N, f):
    h = (b - a) / N
    result = 1 / 2 * (f(a) + f(b))

    for k in range(1, N):
        result += f(a + k * h)
    return h * result

def integrand(x):
    return (np.exp(-(x - x0)**2 / (2 * sigma**2)))**2

integral = trapezoidal_method(-L / 2, L / 2, 1000, integrand)
psi_0 = np.sqrt(1 / integral)
psi_0 = 1/((2*np.pi*(sigma**2))**(1/4))

def Psi_initial(x):
    return psi_0 * np.exp(-((x-x0)**2)/(4*sigma**2) + 1j*k*x)

x = np.linspace(-L/2, L/2, P-1)
plt.plot(x,np.real(np.conj(Psi_initial(x))*Psi_initial(x)))
plt.legend(['T=0'])
#plt.ylim(0,4e9)

# Potential
def V(x):
    return 0



def hamiltonian_matrix(V, N):
    H = np.zeros((N - 1, N - 1))
    A = -sc.hbar**2/(2*m*(L/P)**2)
    
    for i in range(N-1):
        H[i-1][i-1] = V(i*L/P - L/2) - 2*A
        H[i-1][i-2] = A
        H[i-1][i] = A
    return H

H_matrix = hamiltonian_matrix(V, P)
print(H_matrix)

def hamiltonian_matrix(V, N):
    '''
    H = np.zeros((N - 1, N - 1))
    A = -sc.hbar**2/(2*m*(L/P)**2)
    
    for i in range(N-1):
        H[i-1][i-1] = V(i*L/P - L/2) - 2*A
        H[i-1][i-2] = A
        H[i-1][i] = A
    ''' 
    A = -(sc.hbar**2)/(2*m*(interval**2))
    B_p = np.zeros(P-1)
    for p in range(P-1):
        B_p[p] = V(p*interval - L/2) - (2*A)
    H_D = np.zeros((P-1,P-1))
    for a in range(P-1):
        for b in range(P-1):
            if a == b:
                H_D[a,b] = B_p[a]
            elif b == a - 1 or b == a + 1:
                H_D[a,b] = A
    return H_D

H_matrix = hamiltonian_matrix(V, P)
print(H_matrix)



# Initialize L matrix
# L_matrix = 1j*(tau/(2*sc.hbar))*H_matrix + np.identity(P-1)

# Initialize R matrix
# R_matrix = -1j*(tau/(2*sc.hbar))*H_matrix + np.identity(P-1)

I = np.diag(np.full(P-1, 1.0), 0)
i = complex(0,1) 
L_matrix = I + (i*tau/(2*sc.hbar))*H_matrix
R_matrix = I - (i*tau/(2*sc.hbar))*H_matrix

print(L_matrix, R_matrix)

psi = np.zeros((P-1,N), dtype=complex)
psi[:,0] = Psi_initial(x)
t_array = np.linspace(0,T,N)

print(psi[:,0])


for i in range(N-1):
    # Here we solve the system of equations
    # R*Phi^(n+1) = L*Phi^(n)
    # as we iterate over time.
    # L dot psi = v 
    v = np.dot(R_matrix, psi[:,i])
    psi[:,i+1] = np.linalg.solve(L_matrix, v)
    #Phi[:,a+1] = Phi[:,a+1]/normalize(Phi[:,a+1])
    
    print("{} %".format((i/(N-1))*100))


def pdf(Psi):
    return np.real(np.conj(Psi)*Psi)

plt.figure(figsize=(10,7))
plt.plot(x, pdf(psi[:,0]), label = "$t = 0$")
#plt.plot(x, pdf(psi[:,350]), label = "$t = 0->T/4$")
plt.plot(x, pdf(psi[:,750]), label = "$t = T/4$")
#plt.plot(x, pdf(psi[:,900]), label = "$t = T/4->T/2$")
plt.plot(x, pdf(psi[:,1500]), label = "$t = T/2$")
plt.plot(x, pdf(psi[:,N-1]), label = "$t = T$")
plt.xlabel("x (m)")
plt.ylabel("$|\psi(x,t)|^2$")
plt.title("Probability Density")
plt.grid()
plt.legend()


