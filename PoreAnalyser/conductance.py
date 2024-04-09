import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

e = 1.6022*1e-19 # C

T = 300 # Kelvin
J2kcal = 1/4184 # 1 kcalth = 4184 J
kT = 1.380*1e-23 *T #* J2kcal# J/K

# diffusion (Balme 25deg)
D_Na = 1.334*1e-9 #m^2/s
D_K = 1.8*1e-9 #m^2/s
D_Cl = 2.032*1e-9 #m^2/s
# mobility
mu_Na = D_K/kT
mu_Cl = D_Cl/kT

# concentration
l2m3 = 0.001
mol = 6.022*1e23 # particles
c = 0.15 *1/l2m3 *mol # mol/ m^3
bulk_conductivity = e*e*(mu_Na+mu_Cl)*c
#print('bulk_conductivity', bulk_conductivity)

A2m = 1e-10
pS = 1e-12

popt = [1.40674664, 1.25040698] # trained with all FF and both charged and uncharged
f_size = 14

def bullk_conduct(z,a,b, conduct=bulk_conductivity):
    '''
    Calculate the conductance of a pore using bulk conductivity 
    '''
    z = np.array(z)
    a = np.array(a)
    b = np.array(b)
    R = 0
    R_vec = []
    for i in range(1, len(z)):
        R_slice = (z[i]-z[i-1])/(bulk_conductivity*np.pi*a[i]*b[i]*A2m) 
        if R_slice<0:
            print(i, 'R_slice', R_slice, 'dz', z[i]-z[i-1], 'a[i]', a[i], 'b[i]', b[i] )
        R += R_slice
        R_vec.append(R_slice)
    
    conduct = 1/R # * bulk_conductivity
    return conduct, R_vec

def no_bulk_conduct(z,a,b, popt, conduct=bulk_conductivity, plot=False, shift_bond=1.85):
    '''
    Calculate the conductance of a pore using a conductivity model 
    '''
    z = np.array(z)
    a = np.array(a)
    b = np.array(b)

    x1 = popt[0]
    x2 = popt[0]
    shift = popt[1]

    def sigmoid_2d(x, x1, x2, shift):
        x,y = x[0], x[1]
        s1 = 1/(1+np.exp(-x*x1+shift))
        s2 = 1/(1+np.exp(-y*x2+shift))
        return s1*s2
    R = 0
    R_vec = []
    facs = []
    for i in range(1, len(z)):
        convert_unit = 10
        fac = sigmoid_2d([(b[i]+shift_bond)/convert_unit, (a[i]+shift_bond)/convert_unit], x1, x2, shift) # model was trained in NM not in Angstrom!!!
        facs.append(fac)
        R_slice = (z[i]-z[i-1])/(bulk_conductivity*fac* np.pi*a[i]*b[i]*A2m)
        R += R_slice
        R_vec.append(R_slice)
    if plot:    
        fig = plt.figure()#figsize=(4,4)
        ax = fig.add_subplot(111)
        plt.plot(z[1:], facs)
        ax.set_xlabel("z ($\AA$)", fontsize=f_size)
        ax.set_ylabel(r"rel. factor bulk conductivity", fontsize=f_size)
        ax.tick_params(axis='both', which='major', labelsize=f_size)
        fig.tight_layout()
        plt.show()
    return 1/R, R_vec,facs