# -*- coding: utf-8 -*-
"""
Created on Mon May 07 13:42:47 2018

@author: Yoshi
"""

# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:36:48 2018

@author: Yoshi
"""

import numpy as np
import matplotlib.pyplot as plt
import tellurium as te
import time

r = te.loada("""
# This model describes simple central dogma
    J0: -> M ; a_m        #production of mRNA
    J1: M -> ; d_m*M      #degradation of mRNA 
    J2: M -> P ; a_p*M    #production of protein
    J3: P -> ; d_p*P      #degradation of protein
    
    # Parameters
    a_m = 10; d_m = 1     #prod and deg of mRNA
    a_p = 500; d_p = 0.05 #prod and deg of protein
    
    # Initial values (every variable needs I.V.)
    M = 0
    P = 0
""")

#clears out any previous simulation results if it exists
r.reset() 

#the simulation is run, and is saved in "r.model", 
#but results are also stored in "result."
result = r.simulate(0,200,1000) #(start, end, timepoints)

# Find the steady-state value of Protein
P_ss = result[-1,2]

r.plot() # quick plot of model results
#You can do something like the following to alter the quick plot
#r.plot(ylim=[0,300])

### Plotting "result" using matplotlib (more customizable)
#plot mRMA
plt.figure(1)
plt.subplot(211)
plt.plot(result[:,0], result[:,1], 'r', label = "mRNA")
plt.ylabel('Copies', fontsize='12') 
plt.xlabel('Time (hr)', fontsize='12')
plt.legend(loc=4)
#plot protein
plt.subplot(212)
plt.plot(result[:,0], result[:,2], 'b', label = "Protein")
plt.plot([0,200], [P_ss,P_ss], 'k--', label = "Steady-State")
plt.ylim([0,120000])
plt.ylabel('Copies', fontsize='12') 
plt.xlabel('Time (hr)', fontsize='12')
plt.legend(loc=4)
plt.tight_layout()

#%% Doing a parametric sweep using a for loop

# Investigate the input output relationship between 
# steady state Protein levels (P_ss)
# and mRNA production rate (a_m)
# or mRNA decay rate (d_m)

a_m = np.linspace(1,50,50)  
d_m = np.array([0.0001,0.001,0.01,0.1,1,10,100])
P_ss = np.zeros(len(d_m))
T_r = np.zeros(len(d_m))
count = 0
for i in d_m:
    r = te.loada("""
        J0: -> M ; a_m
        J1: M -> ; d_m*M
        J2: M -> P ; a_p*M
        J3: P -> ; d_p*P
        
        # Parameters
        a_m = 10;
        a_p = 500; d_p = 0.05
        
        # Initial values
        M = 0
        P = 0
        """
        # Iterating parameters
        "d_m = " + str(i) + ";"
    )
    result = r.simulate(0,200,1000)
    
    # Find the steady-state value
    P_ss[count] = result[-1,2]
    
#    print 'd_m =',d_m[count] #verbose print of current d_m value
    
    plt.figure(3)
    plt.semilogy(result[:,0],result[:,2],label=d_m[count])
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("Protein")
    
    plt.figure(4)
    plt.plot(result[:,0],result[:,1],label=d_m[count])
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("mRNA")
    
    time.sleep(.25)
    
### or you can print them one at a time.
#    if d_m[count] < 1:
#        r.plot(ylim=[0,10000000])
#    if d_m[count] >= 1:
#        r.plot(ylim=[0,100000])
#    time.sleep(.25)

    count += 1
    
#compare steady-state levels
plt.figure(5)
plt.hist(P_ss)
plt.xlabel("Protein value at Steady-state")

#%%%%%%%%%% Metabolic networks

#### Assume there is no external degradation unless specified.
#### Assume all production rates are at unity. 
#### Assume boundary species are set to 10 counts.
#### Assume initial values of floating species are at 0.

# Models 1 + 2 : https://i.imgur.com/KpeDxJe.jpg
# Model 3 : https://i.imgur.com/kWZZ3PY.jpg


# Make a model of the first linear pathway

r = te.loada("""
    J0: $X0 -> S1 ; v_1 * X0    #you can put a $ sign before a species to set it to a boundary species (concentration doesn't change)
    J1: S1 -> S2 ; v_2 * S1
    J2: S2 -> S3 ; v_3 * S2 
    J3: S3 -> S4 ; v_4 * S3
    J4: S4 -> $X1 ; v_5*S4 
    
    # Parameters
    v_1 = 1 
    v_2 = 3
    v_3 = 5
    v_4 = 2
    v_5 = 1
    
    # Initial values
    X0 = 10
    X1 = 10
    S1 = 0
    S2 = 0
    S3 = 0
    S4 = 0
"""
)
r.reset()
result = r.simulate(0,20,200)
r.plot()
r.getSteadyStateValues()

#%% 

# Make a model of the second pathway with negative feedback.

a = te.loada("""
    J0: $X0 -> S1 ; v_1 * X0 * (K/(K+S3))
    J1: S1 -> S2 ; v_2 * S1
    J2: S2 -> S3 ; v_3 * S2 
    J3: S3 -> $X1 ; v_4 * S3 
    # Parameters
    K = 10
    v_1 = 1
    v_2 = 1
    v_3 = 1
    v_4 = 1
    # Initial values
    X0 = 5
    X1 = 10
    S1 = 0
    S2 = 0
    S3 = 0
"""
)
a.reset()
a.draw()
result = a.simulate(0,30,200)
a.plot()
    
    
#%% 

# Make a model of this third complex pathway.

b = te.loada("""
    J0: -> S1 ; v1
    J1: S1 -> S2 ; v2*S1
    J2: S2 -> ; v3*S2
    J3: S3 -> S1 ; v4*S3
    J4: S3 -> S2 ; v5*S3
    J5: -> S3 ; v6
    
    # Parameters
    v1 = 1
    v2 = 1
    v3 = 1
    v4 = 1
    v5 = 1
    v6 = 1
    
    # Initial values
    S1 = 0
    S2 = 0
    S3 = 0
    
"""
)
b.reset()
result = b.simulate(0,30,200)
b.plot()

#%% Use Libstructural to analyze model b

import structural as st
ls = st.LibStructural()
ls.loadSBMLFromString(b.getSBML())
print(ls.getSummary())
print(ls.getTestDetails())
print ls.getStoichiometryMatrix()
print ls.getElementaryModes()