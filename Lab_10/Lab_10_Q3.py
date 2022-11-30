from random import random 
import matplotlib.pyplot as plt
import numpy as np

def Mean_value(func, sample):
    '''uses the mean value method to evaluate an integral on bounds a=0...b=1
    input:
        func: integrand to be evaluated
        sample: number of sample points on the interval
    output:
        integral: resultant integral'''
    a = 0; b = 1
    f_avg = 0

    for i in range(sample): 
        x = (b-a)*random()
        f_avg += func(x) 

    Integral = (b-a)*f_avg/sample
    return Integral

def Imp_Sampling(func,weight,sample):
    '''calculates an integral using importance sampling on interval a = 0...b = 1
    imput:
        func: integrand
        weight: weighting function, removes the singularity
        sample: number of sample points
    output:
        integral: resultant integral'''
    a = 0
    b = 1

    def p(x):
        return 1/(2*(x**(1/2)))
    Sigma = 0

    for i in range(sample): 
        x = ((b-a)*random())**2
        Sigma += 2*func(x) / weight(x)

    Integral = Sigma/sample

    return Integral

# define integrands
def f(x):
    return (x**(-1/2))/(1+np.exp(x))
def w(x):
    return (x**(-1/2))

# define number of samples
N_mean = 10000
N_impor = 100

# calculate integral using mean value method
mean_integral = Mean_value(f, N_mean)
impor_integral = Imp_Sampling(f,w,N_impor)
print(mean_integral)
print(impor_integral)

impor_sample = []
mean_val = []
for i in range(100): # running the algorithm multiple times to look for occurence for the histograms
    impor_sample.append(impor_integral)
    mean_val.append(mean_integral)

# plotting our histograms for both methods
fig,(ax0,ax1) = plt.subplots(figsize=(10,4),ncols=2,nrows=1)

ax0.hist(impor_sample,density=True,bins=10,range=(0.80,0.88))
ax0.set_title('Importance sampling histogram')
ax0.set_xlabel('Values (10 bins)')
ax0.set_ylabel('Occurence')
ax1.hist(mean_val,density=True,bins=10,range=(0.80,0.88))
ax1.set_title('Mean value method histogram')
ax1.set_xlabel('Values (10 bins)')
ax1.set_ylabel('Occurence')
plt.show()

# wolfram result evaluates to 0.838932960...