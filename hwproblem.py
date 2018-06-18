# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 12:23:37 2018

@author: Veronica Porubsky
"""
import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import time
#%%

r = te.loada("""
# This model simulates simple chemical networks
    J0: X+B -> X+X ; k_1*X*B        #reaction 1
    J1: Y+B -> Y+Y ; k_2*Y*B      #reaction 2
    J2: X+Y -> B+B ; k_3*X*Y   #reaction 3
  
    
    # Parameters
    k_1 = 1
    k_2 = 1
    k_3 = 1
    
    # Initial values (every variable needs I.V.)
    X = 0.2
    Y = 0.4
    B = 0.4
""")

#clears out any previous simulation results if it exists
r.reset() 

#the simulation is run, and is saved in "r.model", 
#but results are also stored in "result."

result = r.simulate(0,10,1000)     # (start, end, timepoints)

r.plot()     # quick plot of model results
#You can do something like the following to alter the quick plot
#r.plot(ylim=[0,300])

### Plotting "result" using matplotlib (more customizable than model.plot())
#plot mRMA
plt.figure(1)
plt.plot(result[:,0], result[:,1],  label = "c0_x(t)", marker = None, linestyle='-', color='blue')
plt.plot(result[:,0], result[:,2],  color = 'green', label = "c0_b(t)")
plt.plot(result[:,0], result[:,3],  color = 'red', label = "c0_y(t)")
plt.title('Analyte Concentration vs. Time')
plt.ylim([0,1.1])
plt.ylabel('Concentration (M)', fontsize='12') 
plt.xlabel('Time (s)', fontsize='12')
plt.legend(loc= 'upper left')
plt.grid(axis='both')
plt.tight_layout()
#%%
r = te.loada("""
# This model simulates simple chemical networks
    J0: A+B -> 2A ; k_1*A*B        #reaction 1
    J1: C+A -> 2C ; k_2*C*A      #reaction 2
    J2: B+C -> 2B ; k_3*B*C   #reaction 3
  
    
    # Parameters
    k_1 = 1
    k_2 = 1
    k_3 = 1
    
    # Initial values (every variable needs I.V.)
    A = 0.2
    B = 0.3
    C = 0.5
""")
#clears out any previous simulation results if it exists
r.reset() 

#the simulation is run, and is saved in "r.model", 
#but results are also stored in "result."
result = r.simulate(0,10,1000)     # (start, end, timepoints)

r.plot()     # quick plot of model results
#You can do something like the following to alter the quick plot
#r.plot(ylim=[0,300])

### Plotting "result" using matplotlib (more customizable than model.plot())
#plot mRMA
plt.figure(1)
plt.plot(result[:,0], result[:,1], label = "c0_a(t)", marker = None, linestyle='-', color='blue')
plt.plot(result[:,0], result[:,2], color = 'green', label = "c0_b(t)")
plt.plot(result[:,0], result[:,3], color = 'red', label = "c0_c(t)")
plt.title('Analyte Concentration vs. Time')
plt.ylim([0.15,0.55])
plt.ylabel('Concentration (M)', fontsize='12') 
plt.xlabel('Time (s)', fontsize='12')
plt.legend(loc= 'upper right')
plt.grid(axis='both')
plt.tight_layout()

#%% Doing a parametric sweep of mRNA production rate using a for loop

# Investigate the input output relationship between 
# steady state Protein levels (P_ss)
# and mRNA production rate (a_m)
# or mRNA decay rate (d_m)
d_m = np.linspace(1,50,50)
a_m = np.array([0.0001,0.001,0.01,0.1,1,10,100])
P_ss = np.zeros(len(a_m))
T_r = np.zeros(len(a_m))
count = 0
for i in a_m:
    r = te.loada("""
        J0: -> M ; a_m
        J1: M -> ; d_m*M
        J2: M -> P ; a_p*M
        J3: P -> ; d_p*P
        
        # Parameters
        d_m = 10;
        a_p = 500; d_p = 0.05;
        
        # Initial values
        M = 0
        P = 0
        """
        # Iterating parameters
        "a_m = " + str(i) + ";"
    )
    result = r.simulate(0,0.01,100)
    
    # Find the steady-state value
    P_ss[count] = result[-1,2]
    
#    print 'a_m =',a_m[count] #verbose print of current a_m value
    
    plt.figure(2)
    plt.semilogy(result[:,0],result[:,2],label=a_m[count])
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("Protein")
    
    plt.figure(3)
    plt.plot(result[:,0],result[:,1],label=a_m[count])
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("mRNA")
    
    time.sleep(.25)
    
### or you can print them one at a time.
#    if a_m[count] < 1:
#        r.plot(ylim=[0,10000000])
#    if a_m[count] >= 1:
#        r.plot(ylim=[0,100000])
#    time.sleep(.25)

    count += 1
    
#compare steady-state levels
plt.figure(4)
plt.hist(np.log10(P_ss))
plt.xlabel("Protein Value at Steady-state on Log Base 10 Scale")

#%% Doing a parametric sweep of protein production rate using a for loop

# Investigate the input output relationship between 
# steady state Protein levels (P_ss)
# and mRNA production rate (a_m)
# or mRNA decay rate (d_m)

a_m = np.linspace(1,50,50)   
a_p = np.array([0.0001,0.001,0.01,0.1,1,10,100])
M_ss = np.zeros(len(a_p))
T_r = np.zeros(len(a_p))
count = 0
for i in a_p:
    r = te.loada("""
        J0: -> M ; a_m
        J1: M -> ; d_m*M
        J2: M -> P ; a_p*M
        J3: P -> ; d_p*P
        
        # Parameters
        a_m = 10; d_m = 10;
        d_p = 0.05;
        
        # Initial values
        M = 0
        P = 0
        """
        # Iterating parameters
        "a_p = " + str(i) + ";"
    )
    result = r.simulate(0,0.6,100)
    
    # Find the mRNA steady-state value
    M_ss[count] = result[-1,1]
    
#    print 'a_m =',a_m[count] #verbose print of current a_p value
    
    plt.figure(5)
    plt.semilogy(result[:,0],result[:,2],label=a_p[count])
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("Protein")
    
    plt.figure(6)
    plt.plot(result[:,0],result[:,1],label=a_p[count])
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("mRNA")
    
    time.sleep(.25)
    
### or you can print them one at a time.
#    if a_p[count] < 1:
#        r.plot(ylim=[0,10000000])
#    if a_p[count] >= 1:
#        r.plot(ylim=[0,100000])
#    time.sleep(.25)

    count += 1
    
#compare steady-state levels
plt.figure(7)
plt.hist((M_ss))
plt.xlabel("mRNA Value at Steady-state")