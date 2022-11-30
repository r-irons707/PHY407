# Contributors 
# Patrick Sandoval
# Kelvin Leong

######################################## HEADER ########################################
# This python script holds the code, pseudo-code and plots for Q2 of Lab07 for 
# finding the energy eigenvalues and radial eigenfunctions of the Hydrogen atom.
######################################### Q1b/c ###########################################
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from MyFunctions import solve_En, solve_Rn, En
from scipy.integrate import simps

# Constant
e = 1.6022e-19
a = 5e-11 # m (Bohr radius)
h = 0.0001*a # step size
rinf = 20*a

# Main program to find the energy using the secant method
# Define n for energy level and l for degree
n = [1,2]
l = [0,1]

# Set figure for plotting and plotting settings
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
fig, a0 = plt.subplots(figsize=(8,6))
color = ['k','r','gray']
i = 0

# Iterate over energy levels
for nval in n:
    # Iterate over degrees
    for lval in l:
        # Avoid equal degree and energy level
        if lval == nval:
            continue
        # Print statement to identify case
        print(f"Energy Eigvanlues for l={lval} and n={nval}")
        # Set initial enegry eigenvalues
        E1, E2 = En(nval)
        # Solve for R for initial guess of E
        psi2 = solve_En(E1,lval)
        # Set target accuracy of result
        target = e/1000
        # Set condition for while loop
        while abs(E1-E2)>target:
            # Compute next instance of R given next E
            psi1,psi2 = psi2,solve_En(E2,lval)
            # Compute next E
            E1,E2 = E2,E2-psi2*(E2-E1)/(psi2-psi1)
            # Print estimated E
            print("E =" ,E2/e, "eV")

        # Once we know eigenvalue we compute eigenfunction
        psi2 = solve_Rn(E2,lval)
        # Generate evenly space position array
        r_domain = np.arange(h,rinf,h)
        # Normalize eigenfunction
        C = simps(np.power(psi2,2),r_domain)
        psi2 /= np.sqrt(C)
        print("|Psi|^2 = ",simps(np.power(psi2,2),r_domain))
        a0.plot(r_domain,psi2,c=color[i],label=f'R(r) for l={lval}, n={nval}')
        a0.legend()
        a0.set_xlabel("Radial Distance (m)",fontsize=16)
        a0.set_ylabel("R(r)",fontsize=16)
        a0.set_title("Hydrogen Atom Radial Solution",fontsize=18)
        i += 1
        print()
a0.xaxis.set_minor_locator(MultipleLocator(2e-11))
a0.yaxis.set_minor_locator(MultipleLocator(10e3))
a0.yaxis.set_ticks_position('both') 
a0.xaxis.set_ticks_position('both')
plt.tight_layout()
plt.savefig("Q2Plot.pdf")
