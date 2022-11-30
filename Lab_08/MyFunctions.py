import numpy as np
def Fmu(mu,eta):
    g = 9.8
    return(0.5*mu**2 + g*eta)

def Feta(mu,eta):
    return(mu*eta)